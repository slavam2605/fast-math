//
// Created by Vyacheslav.Moklev on 09/03/2024.
//

#include "big_int/big_int_ops.h"

#include "big_int/burnikel_ziegler.h"
#include "big_int/karatsuba.h"
#include "big_int/toom_cook.h"

/**
 * Fast implementation for `a = a - (b * c << offset * 64)`, optimized for Knuth's division.
 * Does not normalize output intentionally.
 */
uint64_t sub_mul_abs_uint64(bint_t& a, const bint_t& b, const uint64_t c, const int offset) {
    if (offset > a.data.size()) throw std::runtime_error("sub_mul_abs_uint64_special: offset > a.data.size()");
    if (c == 0) return 0; // a - 0 == a

    const __uint128_t factor = c;
    uint64_t carry = 0;
    for (int i = offset; i < a.data.size() || i - offset < b.data.size(); i++) {
        __uint128_t _mc = a.data[i];
        _mc -= carry;
        if (i - offset < b.data.size()) _mc -= factor * b.data[i - offset];
        carry = -static_cast<uint64_t>(_mc >> 64);
        a.data[i] = static_cast<uint64_t>(_mc);
    }

    // Don't normalize `a` -- this is intentional
    return carry;
}

/**
 * Fast implementation for `a = a + ((b[b_start..b.size-1] * c) << (64 * left_shift))`
 */
void add_mul_abs_uint64(bint_t& a, const bint_t& b, const uint64_t c, const int b_start, const int left_shift) {
    if (c == 0) return; // a + 0 == a
    if (left_shift > a.data.size()) throw std::runtime_error("add_mul_abs_uint64: left_shift > a.data.size()");

    const __uint128_t factor = c;
    uint64_t carry = 0;
    for (int i = left_shift; i < a.data.size() || i < b.data.size() + left_shift - b_start; i++) {
        __uint128_t _mc = a.data[i];
        _mc += carry;
        if (i - left_shift + b_start < b.data.size()) {
            _mc += factor * b.data[i - left_shift + b_start];
        }
        carry = static_cast<uint64_t>(_mc >> 64);
        if (i < a.data.size()) {
            a.data[i] = static_cast<uint64_t>(_mc);
        } else {
            a.data.push_back(static_cast<uint64_t>(_mc));
        }
    }

    if (carry > 0) {
        a.data.push_back(carry);
    }
}

bint_t schoolbook_mul_abs(const bint_t& a, const bint_t& b, const int a_limit, const int b_limit) {
    bint_t result;
    result.data.resize(a_limit + b_limit);
    for (int i = 0; i < b_limit; i++) {
        const __uint128_t factor = b.data[i];
        if (factor == 0)
            continue;

        uint64_t carry = 0;
        for (int j = i; j < a_limit + i; j++) {
            __uint128_t c = carry;
            c += factor * a.data[j - i];
            c += result.data[j];
            carry = c >> 64;
            result.data[j] = static_cast<uint64_t>(c);
        }
        result.data[a_limit + i] = carry;
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

/**
 * Knuth's algorithm D for the long division.
 * Implementation is inspired by https://skanthak.hier-im-netz.de/division.html
 */
void divide_knuth_abs_inner(bint_t& a, const bint_t& b, bint_t& rem) {
    if (a.data.back() == 0) throw std::runtime_error("divide_knuth_abs_inner: leading zeroes in divident");
    if (b.data.back() == 0) throw std::runtime_error("divide_knuth_abs_inner: leading zeroes in divisor");

    const int m = a.data.size();
    const int n = b.data.size();

    // Keep `a` unnormalized and normalize only in the end, this is important
    a.data.push_back(0);
    bint_t result;

    for (int i = m - n; i >= 0; i--) {
        const __uint128_t a_part = (static_cast<__uint128_t>(a.data[i + n]) << 64) + a.data[i + n - 1];
        __uint128_t qhat = a_part / b.data[n - 1];
        __uint128_t rhat = a_part % b.data[n - 1];

        while (qhat >> 64 != 0 || static_cast<uint64_t>(qhat) * static_cast<__uint128_t>(b.data[n - 2]) > (rhat << 64) + a.data[i + n - 2]) {
            qhat--;
            rhat += b.data[n - 1];
            if (rhat >> 64 != 0) break;
        }

        if (qhat >> 64 != 0) throw std::runtime_error("divide_knuth_abs_inner: qhat >= 2^64");
        if (sub_mul_abs_uint64(a, b, qhat, i) != 0) {
            qhat--;
            big_int_impl::add_abs_inplace(a, b, 0, -1, i, true);
        }
        result.data.push_back(qhat);

        // Remove last digit of `a`: it is always 0
        if (a.data.back() != 0) throw std::runtime_error("divide_knuth_abs_inner: a.data.back() != 0 after division step");
        a.data.pop_back();
    }

    rem = a;
    big_int_impl::normalize(rem);

    a.data.resize(result.data.size());
    std::copy(result.data.rbegin(), result.data.rend(), a.data.begin());
    big_int_impl::normalize(a);
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

bint_t big_int_impl::schoolbook_square(const bint_t& a, int a_limit) {
    if (a_limit < 0) a_limit = a.data.size();

    bint_t result(0);
    result.data.resize(2 * a_limit);

    uint64_t last_low_word = 0;
    for (int i = a_limit - 1; i >= 0; i--) {
        __uint128_t c = a.data[i];
        c *= c;
        result.data[2 * i] = static_cast<uint64_t>(c >> 1);
        result.data[2 * i + 1] = last_low_word << 63 | static_cast<uint64_t>(c >> 65);
        last_low_word = static_cast<uint64_t>(c);
    }

    for (int i = 1; i < a_limit; i++) {
        add_mul_abs_uint64(result, a, a.data[a_limit - i - 1], a_limit - i, (a_limit - i) * 2 - 1);
    }

    result <<= 1;
    result.data[0] |= a.data[0] & 1;
    normalize(result);
    return result;
}

std::strong_ordering big_int_impl::compare_abs(const bint_t&a, const bint_t& b, const int shift) {
    if (!is_normalized(a)) throw std::runtime_error("compare_abs: a is not normalized");
    if (!is_normalized(b)) throw std::runtime_error("compare_abs: b is not normalized");

    if (a.data.size() < b.data.size() + shift) return std::strong_ordering::less;
    if (a.data.size() > b.data.size() + shift) return std::strong_ordering::greater;

    for (int i = a.data.size() - 1; i >= 0; i--) {
        const auto b_value = i >= shift ? b.data[i - shift] : 0;
        if (a.data[i] < b_value) return std::strong_ordering::less;
        if (a.data[i] > b_value) return std::strong_ordering::greater;
    }

    return std::strong_ordering::equal;
}

void big_int_impl::add_abs_inplace(bint_t& a, const bint_t& b, const int b_start, int b_limit, const int b_shift, const bool let_overflow) {
    if (b_limit < 0) b_limit = b.data.size();
    uint64_t carry = 0;
    for (int i = std::min<int>(b_shift, a.data.size()); i < a.data.size() || i - b_shift < b_limit - b_start; i++) {
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
    if (!let_overflow && carry > 0) {
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

uint64_t big_int_impl::count_bits(const bint_t& a) {
    if (a.data.size() > 1 && a.data.back() == 0) {
        throw std::runtime_error("count_bits: value is not normalized");
    }

    return a.data.size() * 64 - std::countl_zero(a.data.back());
}

void big_int_impl::normalize(bint_t& a) {
    while (a.data.size() > 1 && a.data.back() == 0) {
        a.data.pop_back();
    }
    if (a.data.empty() || a.data.size() == 1 && a.data[0] == 0) {
        a.sign = false;
    }
}

bool big_int_impl::is_normalized(const bint_t& a) {
    return a.data.size() == 1 || !a.data.empty() && a.data.back() > 0;
}

void big_int_impl::square(bint_t& a, int a_limit) {
    if (a_limit < 0) a_limit = a.data.size();

    if (a_limit < big_int::KARATSUBA_SQUARE_THRESHOLD) {
        a = schoolbook_square(a, a_limit);
        return;
    }

    if (a_limit < big_int::TOOM_COOK_SQUARE_THRESHOLD) {
        karatsuba_square(a, a_limit);
        return;
    }

    toom3_square(a, a_limit);
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

    // Check if possible to calculate a square instead
    if (&a == &b && a_limit == b_limit) {
        auto result = a;
        big_int_impl::square(result, a_limit);
        return result;
    }

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
 * Preparation for the Knuth's algorithm D: shift left until the highest digit of the divisor is at least 2^63
 */
void big_int::divide_knuth_abs(bint_t& a, const bint_t& b, bint_t& rem) {
    if (big_int_impl::compare_abs(a, b) < 0) {
        rem = a;
        a = bint_t(0ll);
        return;
    }

    if (b.data.back() >= 1ull << 63) {
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
            a *= a;
            n >>= 1;
        }
        return;
    }

    n--;
    bint_t c = a;
    while (n > 0) {
        if (n % 2 == 1) a *= c;
        n /= 2;
        c *= c;
    }
}

void big_int::shift_left_inplace(bint_t& a, const int64_t shift) { // NOLINT(*-no-recursion)
    if (shift == 0) return;
    if (shift < 0) return shift_right_inplace(a, -shift);
    if (a.data.size() == 1 && a.data.back() == 0) return; // 0 << shift == 0

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