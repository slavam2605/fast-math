//
// Created by Vyacheslav.Moklev on 09/03/2024.
//

#include "big_int_ops.h"

#include <iostream>

#include "karatsuba.h"
#include "toom_cook.h"

std::strong_ordering compare_abs(const bint_t&a, const bint_t& b) {
    if (a.data.size() < b.data.size()) return std::strong_ordering::less;
    if (a.data.size() > b.data.size()) return std::strong_ordering::greater;
    for (int i = a.data.size() - 1; i >= 0; i--) {
        if (a.data[i] < b.data[i]) return std::strong_ordering::less;
        if (a.data[i] > b.data[i]) return std::strong_ordering::greater;
    }
    return std::strong_ordering::equal;
}

void add_abs_inplace(bint_t& a, const bint_t& b, int b_limit, const int b_shift) {
    if (b_limit < 0) b_limit = b.data.size();
    uint64_t carry = 0;
    for (int i = 0; i < a.data.size() || i - b_shift < b_limit; i++) {
        uint64_t c = carry;
        carry = 0;
        if (i < a.data.size()) carry |= __builtin_uaddll_overflow(a.data[i], c, &c);
        if (i >= b_shift && i - b_shift < b_limit) carry |= __builtin_uaddll_overflow(b.data[i - b_shift], c, &c);
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

bint_t add_abs(const bint_t& a, const bint_t& b) {
    bint_t copy = a;
    add_abs_inplace(copy, b);
    return copy;
}

void sub_abs_inplace(bint_t& a, const bint_t& b, int b_limit) {
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
    a.normalize();
}

// a = a - b * c
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
    a.normalize();
}

bint_t sub_abs(const bint_t& a, const bint_t& b) {
    if (compare_abs(a, b) < 0) {
        auto result = sub_abs(b, a);
        result.sign = true;
        return result;
    }

    bint_t copy = a;
    copy.sign = false;
    sub_abs_inplace(copy, b);
    return copy;
}

bint_t slow_mul_abs(const bint_t& a, const bint_t& b, int a_limit, int b_limit) {
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

bint_t mul_abs_uint64(const bint_t& a, const uint64_t b) {
    if (b == 0) return bint_t(0ll);

    bint_t result;
    result.data.reserve(a.data.size() + 1);
    const __uint128_t factor = b;
    uint64_t carry = 0;
    for (const auto& element: a.data) {
        const __uint128_t c = factor * element + carry;
        carry = c >> 64;
        result.data.push_back(static_cast<uint64_t>(c));
    }
    if (carry > 0) result.data.push_back(carry);
    return result;
}

void div_abs_inplace_inner(bint_t& a, const bint_t& b, bint_t& rem) {
    if (compare_abs(a, b) < 0) {
        rem = a;
        a = bint_t(0ll);
        return;
    }

    if (b.data.back() == 0) throw std::runtime_error("div_abs_inplace_inner: leading zeroes in divisor");
    bint_t current;
    for (int i = a.data.size() - 1; i >= 0; i--) {
        current.data.insert(current.data.begin(), a.data[i]);
        if (compare_abs(current, b) < 0) {
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
            if (diff == right - left || compare_abs(current, b) < 0) break;
            sub_abs_inplace(current, b);
            a.data[i]++;
        }
    }
    a.normalize();
    current.normalize();
    rem = current;
}

void shift_left_inplace(bint_t& a, const int shift) {
    if (shift >= 64) throw std::runtime_error("shift_left_inplace: shift >= 64: " + std::to_string(shift));
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

void shift_right_inplace(bint_t& a, const int shift) {
    if (shift >= 64) throw std::runtime_error("shift_right_inplace: shift >= 64: " + std::to_string(shift));
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

// Knuth's idea to shift left until highest digit of divisor is greater than 2^63
void div_abs_inplace(bint_t& a, const bint_t& b, bint_t& rem) {
    if (b.data.back() > 1ull << 63) {
        return div_abs_inplace_inner(a, b, rem);
    }

    bint_t new_b = b;
    const int shift = std::countl_zero(new_b.data.back());
    shift_left_inplace(new_b, shift);
    shift_left_inplace(a, shift);
    div_abs_inplace_inner(a, new_b, rem);
    shift_right_inplace(rem, shift);
}

// Optimized version for division by uint64_t
void div_abs_inplace(bint_t& a, const uint64_t b, uint64_t& rem) {
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

void fast_pow_inplace(bint_t& a, uint64_t n) {
    if (n == 0) {
        a = bint_t(1ll);
        return;
    }

    // n is power of 2
    if ((n & (n - 1)) == 0) {
        while (n > 1) {
            a = karatsuba(a, a);
            n >>= 1;
        }
        return;
    }

    n--;
    bint_t c = a;
    while (n > 0) {
        if (n % 2 == 1) a = karatsuba(a, c);
        n /= 2;
        c = karatsuba(c, c);
    }
}

bint_t add(const bint_t& a, const bint_t& b) {
    if (!a.sign && b.sign) return sub_abs(a, b);
    if (a.sign && !b.sign) return sub_abs(b, a);
    if (a.sign) { // => b.sign
        auto result = add_abs(a, b);
        result.sign = true;
        return result;
    }
    return add_abs(a, b);
}

bint_t sub(const bint_t& a, const bint_t& b) {
    if (!a.sign && b.sign) return add_abs(a, b);
    if (a.sign && !b.sign) {
        auto result = add_abs(a, b);
        result.sign = true;
        return result;
    }
    if (a.sign) return sub_abs(b, a);
    return sub_abs(a, b);
}

bint_t multiply(const bint_t& a, const bint_t& b) {
    if (a.data.size() < KARATSUBA_THRESHOLD || b.data.size() < KARATSUBA_THRESHOLD)
        return slow_mul(a, b); // TODO check a.size == 1 || b.size == 1

    if (a.data.size() < TOOM_COOK_THRESHOLD && b.data.size() < TOOM_COOK_THRESHOLD)
        return karatsuba(a, b);

    return toom3(a, b);
}

bint_t slow_mul(const bint_t& a, const bint_t& b, int a_limit, int b_limit) {
    if (a_limit < 0) a_limit = a.data.size();
    if (b_limit < 0) b_limit = b.data.size();
    auto result = slow_mul_abs(a, b, a_limit, b_limit);
    result.sign = a.sign ^ b.sign;
    result.normalize();
    return result;
}

std::strong_ordering compare(const bint_t& a, const bint_t& b) {
    if (a.sign && !b.sign) return std::strong_ordering::less;
    if (!a.sign && b.sign) return std::strong_ordering::greater;
    const auto result = compare_abs(a, b);
    return a.sign ? 0 <=> result : result;
}

void debug_print(const bint_t& a) {
    std::cout << (a.sign ? "-" : "+") << "[";
    for (int i = a.data.size() - 1; i >= 0; i--) {
        std::cout << a.data[i];
        if (i > 0) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
}