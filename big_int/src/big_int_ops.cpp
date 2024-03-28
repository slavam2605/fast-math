//
// Created by Vyacheslav.Moklev on 09/03/2024.
//

#include "big_int/big_int_ops.h"

#include "big_int/burnikel_ziegler.h"
#include "big_int/karatsuba.h"
#include "big_int/toom_cook.h"

/**
 * Fast implementation for `a = a - b * c`
 */
void sub_mul_abs_uint64(bint_t& a, const bint_t& b, const uint64_t c) {
    if (c == 0) return; // a - 0 == a

    const __uint128_t factor = c;
    uint64_t carry = 0;
    for (int i = 0; i < a.data.size() || i < b.data.size(); i++) {
        __uint128_t _mc = a.data[i];
        _mc -= carry;
        if (i < b.data.size()) _mc -= factor * b.data[i];
        carry = -static_cast<uint64_t>(_mc >> 64);
        a.data[i] = static_cast<uint64_t>(_mc);
    }

    if (carry > 0) throw std::runtime_error("sub_mul_abs_uint64: carry > 0");
    big_int_impl::normalize(a);
}

bint_t schoolbook_mul_abs(const bint_t& a, const bint_t& b, int a_limit, int b_limit) {
    bint_t result(0ull);
    for (int i = 0; i < b_limit; i++) {
        const __uint128_t factor = b.data[i];
        if (factor == 0)
            continue;

        uint64_t carry = 0;
        for (int j = 0; j < a_limit + i; j++) {
            __uint128_t c = carry;
            if (j - i >= 0 && j - i < a_limit) c += factor * a.data[j - i];
            if (j < result.data.size()) c += result.data[j];
            carry = c >> 64;
            if (j < result.data.size()) {
                result.data[j] = static_cast<uint64_t>(c);
            } else {
                result.data.push_back(static_cast<uint64_t>(c));
            }
        }
        if (carry > 0) {
            if (result.data.size() != a_limit + i) throw std::runtime_error("slow_mul_abs: unexpected result.data size");
            result.data.push_back(carry);
        }
    }
    return result;
}

bint_t mul_uint64(const bint_t& a, const uint64_t b, int a_limit, const bool new_sign) {
    if (b == 0) return bint_t(0ll);
    if (a_limit < 0) a_limit = a.data.size();

    bint_t result;
    result.data.reserve(a_limit + 1);
    const __uint128_t factor = b;
    uint64_t carry = 0;
    for (int i = 0; i < a_limit; i++) {
        const auto& element = a.data[i];
        const __uint128_t c = factor * element + carry;
        carry = c >> 64;
        result.data.push_back(static_cast<uint64_t>(c));
    }
    if (carry > 0) result.data.push_back(carry);
    result.sign = new_sign;
    return result;
}

void divide_knuth_abs_inner(bint_t& a, const bint_t& b, bint_t& rem) {
    if (big_int_impl::compare_abs(a, b) < 0) {
        rem = a;
        a = bint_t(0ll);
        return;
    }

    if (a.data.back() == 0) throw std::runtime_error("divide_knuth_abs_inner: leading zeroes in divident");
    if (b.data.back() == 0) throw std::runtime_error("divide_knuth_abs_inner: leading zeroes in divisor");
    bint_t current;
    for (int i = a.data.size() - 1; i >= 0; i--) {
        current.data.insert(current.data.begin(), a.data[i]);
        if (big_int_impl::compare_abs(current, b) < 0) {
            a.data[i] = 0ull;
            continue;
        }

        // Guess the bounds
        uint64_t left;
        uint64_t right;
        if (b.data.size() < current.data.size()) {
            const __uint128_t current_part = (static_cast<__uint128_t>(current.data.back()) << 64) + current.data[current.data.size() - 2];
            left = current_part / (b.data.back() + 1);
            right = (current_part + 1) / b.data.back();
        } else {
            left = current.data.back() / (b.data.back() + 1);
            right = (current.data.back() + 1) / b.data.back();
        }

        if (right - left > 2) throw std::runtime_error("right - left > 2: " + std::to_string(right - left));
        sub_mul_abs_uint64(current, b, left);
        a.data[i] = left;
        for (int diff = 0; diff <= right - left; diff++) {
            if (diff == right - left || big_int_impl::compare_abs(current, b) < 0) break;
            big_int_impl::sub_abs_inplace(current, b);
            a.data[i]++;
        }
    }
    big_int_impl::normalize(a);
    big_int_impl::normalize(current);
    rem = current;
}

bint_t add_abs(const bint_t& a, const bint_t& b) {
    bint_t copy = a;
    big_int_impl::add_abs_inplace(copy, b);
    return copy;
}

bint_t sub_abs(const bint_t& a, const bint_t& b) { // NOLINT(*-no-recursion)
    if (big_int_impl::compare_abs(a, b) < 0) {
        auto result = sub_abs(b, a);
        result.sign = true;
        return result;
    }

    bint_t copy = a;
    copy.sign = false;
    big_int_impl::sub_abs_inplace(copy, b);
    return copy;
}

void shift_left_inplace_blocks(bint_t& a, const int shift_blocks) {
    if (shift_blocks == 0) return;
    a.data.resize(a.data.size() + shift_blocks);
    for (int i = a.data.size() - 1; i >= shift_blocks; i--) {
        a.data[i] = a.data[i - shift_blocks];
    }
    std::fill_n(a.data.begin(), shift_blocks, 0);
}

void shift_right_inplace_blocks(bint_t& a, const int shift_blocks) {
    if (shift_blocks == 0) return;
    if (shift_blocks >= a.data.size()) {
        a = bint_t(0);
        return;
    }

    for (int i = 0; i < a.data.size() - shift_blocks; i++) {
        a.data[i] = a.data[i + shift_blocks];
    }
    a.data.resize(a.data.size() - shift_blocks);
}

void shift_left_inplace_in_block(bint_t& a, const int shift) {
    if (shift == 0) return;
    const uint64_t last = a.data.back() >> (64 - shift);
    for (int i = a.data.size() - 1; i >= 0; i--) {
        const uint64_t self = a.data[i] << shift;
        const uint64_t next = i == 0 ? 0 : a.data[i - 1] >> (64 - shift);
        a.data[i] = self | next;
    }
    if (last != 0) {
        a.data.push_back(last);
    }
}

void shift_right_inplace_in_block(bint_t& a, const int shift) {
    if (shift == 0) return;
    for (int i = 0; i < a.data.size(); i++) {
        const uint64_t self = a.data[i] >> shift;
        const uint64_t prev = i == a.data.size() - 1 ? 0 : a.data[i + 1] << (64 - shift);
        a.data[i] = self | prev;
    }
    if (a.data.size() > 1 && a.data.back() == 0) {
        a.data.pop_back();
    }
}

// ================ namespace big_int_impl ================

std::strong_ordering big_int_impl::compare_abs(const bint_t&a, const bint_t& b, const int shift) {
    if (a.data.size() < b.data.size() + shift) return std::strong_ordering::less;
    if (a.data.size() > b.data.size() + shift) return std::strong_ordering::greater;
    for (int i = a.data.size() - 1; i >= 0; i--) {
        auto b_value = i >= shift ? b.data[i - shift] : 0;
        if (a.data[i] < b_value) return std::strong_ordering::less;
        if (a.data[i] > b_value) return std::strong_ordering::greater;
    }
    return std::strong_ordering::equal;
}

void big_int_impl::add_abs_inplace(bint_t& a, const bint_t& b, const int b_start, int b_limit, const int b_shift) {
    if (b_limit < 0) b_limit = b.data.size();
    uint64_t carry = 0;
    for (int i = 0; i < a.data.size() || i - b_shift < b_limit - b_start; i++) {
        uint64_t c = carry;
        carry = 0;
        if (i < a.data.size()) carry |= __builtin_uaddll_overflow(a.data[i], c, &c);
        if (i >= b_shift && i - b_shift < b_limit - b_start) carry |= __builtin_uaddll_overflow(b.data[i - b_shift + b_start], c, &c);
        if (i < a.data.size()) {
            a.data[i] = c;
        } else {
            a.data.push_back(c);
        }
    }
    if (carry > 0) {
        a.data.push_back(carry);
    }
}

void big_int_impl::sub_abs_inplace(bint_t& a, const bint_t& b, int b_limit) {
    if (b_limit < 0) b_limit = b.data.size();
    if (b_limit > a.data.size()) throw std::runtime_error("sub_abs_inplace: b_limit > a.data.size()");
    uint64_t carry = 0;
    for (int i = 0; i < a.data.size() || i < b_limit; i++) {
        __uint128_t c = a.data[i];
        c -= carry;
        if (i < b_limit) c -= b.data[i];
        carry = -static_cast<uint64_t>(c >> 64);
        a.data[i] = c;
    }
    if (carry > 0) throw std::runtime_error("sub_abs_inplace: carry > 0");
    normalize(a);
}

/**
 * Optimized version for division by uint64_t
 */
void big_int_impl::div_abs_inplace(bint_t& a, const uint64_t b, uint64_t& rem) {
    __uint128_t current = 0;
    for (int i = a.data.size() - 1; i >= 0; i--) {
        current = (current << 64) + a.data[i];
        a.data[i] = current / b;
        current %= b;
    }
    while (a.data.size() > 1 && a.data.back() == 0) {
        a.data.pop_back();
    }
    rem = current;
}

bint_t big_int_impl::schoolbook_multiply(const bint_t& a, const bint_t& b, int a_limit, int b_limit) {
    if (a_limit < 0) a_limit = a.data.size();
    if (b_limit < 0) b_limit = b.data.size();
    auto result = schoolbook_mul_abs(a, b, a_limit, b_limit);
    result.sign = a.sign ^ b.sign;
    normalize(result);
    return result;
}

void big_int_impl::normalize(bint_t& a) {
    while (a.data.size() > 1 && a.data.back() == 0) {
        a.data.pop_back();
    }
    if (a.data.empty() || a.data.size() == 1 && a.data[0] == 0) {
        a.sign = false;
    }
}

// ================ namespace big_int ================

bint_t big_int::add(const bint_t& a, const bint_t& b) {
    if (!a.sign && b.sign) return sub_abs(a, b);
    if (a.sign && !b.sign) return sub_abs(b, a);
    if (a.sign) { // => b.sign
        auto result = add_abs(a, b);
        result.sign = true;
        big_int_impl::normalize(result);
        return result;
    }
    return add_abs(a, b);
}

bint_t big_int::sub(const bint_t& a, const bint_t& b) {
    if (!a.sign && b.sign) return add_abs(a, b);
    if (a.sign && !b.sign) {
        auto result = add_abs(a, b);
        result.sign = true;
        big_int_impl::normalize(result);
        return result;
    }
    if (a.sign) return sub_abs(b, a);
    return sub_abs(a, b);
}

bint_t big_int::multiply(const bint_t& a, const bint_t& b, int a_limit, int b_limit) {
    if (a_limit < 0) a_limit = a.data.size();
    if (b_limit < 0) b_limit = b.data.size();

    if (a_limit < KARATSUBA_THRESHOLD || b_limit < KARATSUBA_THRESHOLD) {
        if (a_limit == 1)
            return mul_uint64(b, a.data[0], b_limit, a.sign ^ b.sign);
        if (b_limit == 1)
            return mul_uint64(a, b.data[0], a_limit, a.sign ^ b.sign);

        return big_int_impl::schoolbook_multiply(a, b, a_limit, b_limit);
    }

    if (a_limit < TOOM_COOK_THRESHOLD && b_limit < TOOM_COOK_THRESHOLD)
        return karatsuba(a, b, a_limit, b_limit);

    return toom3(a, b, a_limit, b_limit);
}

void big_int::divide_abs(bint_t& a, const bint_t& b, bint_t& rem) {
    if (a.data.size() < BURNIKEL_ZIEGLER_THRESHOLD && b.data.size() < BURNIKEL_ZIEGLER_THRESHOLD) {
        return divide_knuth_abs(a, b, rem);
    }

    return divide_burnikel_ziegler(a, b, rem);
}

/**
 * Knuth's idea to shift left until highest digit of divisor is greater than 2^63
 */
void big_int::divide_knuth_abs(bint_t& a, const bint_t& b, bint_t& rem) {
    if (b.data.back() > 1ull << 63) {
        return divide_knuth_abs_inner(a, b, rem);
    }

    bint_t new_b = b;
    const int shift = std::countl_zero(new_b.data.back());
    new_b <<= shift;
    a <<= shift;
    divide_knuth_abs_inner(a, new_b, rem);
    rem >>= shift;
}

std::strong_ordering big_int::compare(const bint_t& a, const bint_t& b) {
    if (a.sign && !b.sign) return std::strong_ordering::less;
    if (!a.sign && b.sign) return std::strong_ordering::greater;
    const auto result = big_int_impl::compare_abs(a, b);
    return a.sign ? 0 <=> result : result;
}

/**
 * Fast power algorithm, uses O(log n) multiplications
 */
void big_int::fast_pow_inplace(bint_t& a, uint64_t n) {
    if (n == 0) {
        a = bint_t(1ll);
        return;
    }

    // n is power of 2
    if ((n & (n - 1)) == 0) {
        while (n > 1) {
            a = multiply(a, a);
            n >>= 1;
        }
        return;
    }

    n--;
    bint_t c = a;
    while (n > 0) {
        if (n % 2 == 1) a = multiply(a, c);
        n /= 2;
        c = multiply(c, c);
    }
}

void big_int::shift_left_inplace(bint_t& a, const int64_t shift) { // NOLINT(*-no-recursion)
    if (shift == 0) return;
    if (shift < 0) {
        return shift_right_inplace(a, -shift);
    }

    const int shift_blocks = shift / 64;
    const int shift_in_block = shift % 64;
    shift_left_inplace_blocks(a, shift_blocks);
    shift_left_inplace_in_block(a, shift_in_block);
}

void big_int::shift_right_inplace(bint_t& a, const int64_t shift) { // NOLINT(*-no-recursion)
    if (shift == 0) return;
    if (shift < 0) {
        return shift_left_inplace(a, -shift);
    }

    const int shift_blocks = shift / 64;
    const int shift_in_block = shift % 64;
    shift_right_inplace_blocks(a, shift_blocks);
    shift_right_inplace_in_block(a, shift_in_block);
}