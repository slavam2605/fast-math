//
// Created by Vyacheslav.Moklev on 09/03/2024.
//

#include "big_int_ops.h"

#include <iostream>

#include "karatsuba.h"

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
    if (b_limit > a.data.size()) throw std::exception();
    uint64_t carry = 0;
    for (int i = 0; i < a.data.size() || i < b_limit; i++) {
        uint64_t c = carry;
        carry = 0;
        carry |= __builtin_usubll_overflow(a.data[i], c, &c);
        if (i < b_limit) carry |= __builtin_usubll_overflow(c, b.data[i], &c);
        a.data[i] = c;
    }
    while (a.data.size() > 1 && a.data.back() == 0) {
        a.data.pop_back();
    }
    if (carry > 0) throw std::exception();
}

bint_t sub_abs(const bint_t& a, const bint_t& b) {
    if (compare_abs(a, b) < 0) {
        auto result = sub_abs(b, a);
        result.sign = true;
        return result;
    }

    bint_t copy = a;
    sub_abs_inplace(copy, b);
    return copy;
}

bint_t mul_abs(const bint_t& a, const bint_t& b, int a_limit, int b_limit) {
    bint_t result(0ull);
    bint_t temp_result;
    for (int i = 0; i < b_limit; i++) {
        const __uint128_t factor = b.data[i];
        if (factor == 0)
            continue;

        uint64_t carry = 0;
        for (int j = 0; j < a_limit; j++) {
            const __uint128_t c = factor * a.data[j] + carry;
            carry = c >> 64;
            temp_result.data.push_back(static_cast<uint64_t>(c));
        }
        if (carry > 0) temp_result.data.push_back(carry);
        add_abs_inplace(result, temp_result, -1, i);
        temp_result.data.clear();
    }
    return result;
}

void div_abs_inplace(bint_t& a, const bint_t& b, bint_t& rem) {
    if (compare_abs(a, b) < 0) {
        rem = a;
        a = bint_t(0ll);
        return;
    }

    if (b.data.back() == 0) throw std::exception(); // no leading zeroes in divisor
    bint_t current;
    for (int i = a.data.size() - 1; i >= 0; i--) {
        current.data.insert(current.data.begin(), a.data[i]);
        if (compare_abs(current, b) < 0) {
            a.data[i] = 0ull;
            continue;
        }

        // Guess the bounds
        uint64_t left = b.data.size() < current.data.size()
            ? ((static_cast<__uint128_t>(current.data.back()) << 64) + current.data[current.data.size() - 2]) / (b.data.back() + 1)
            : current.data.back() / (b.data.back() + 1);

        uint64_t right = b.data.size() < current.data.size()
             ? ((static_cast<__uint128_t>(current.data.back()) << 64) + current.data[current.data.size() - 2] + 1) / b.data.back()
             : (current.data.back() + 1) / b.data.back();

        const uint64_t initial_right = right;
        while (right - left > 1) {
            const uint64_t middle = left + (right - left) / 2;
            bint_t guess = mul(b, bint_t(middle));
            if (compare_abs(guess, current) <= 0) {
                left = middle;
            } else {
                right = middle;
            }
        }
        if (right > left && right == initial_right) {
            // Check if right bound fits
            bint_t guess = mul(b, bint_t(right));
            if (compare_abs(guess, current) <= 0) {
                sub_abs_inplace(current, guess);
                a.data[i] = right;
                continue;
            }
        }
        bint_t guess = mul(b, bint_t(left));
        sub_abs_inplace(current, guess);
        a.data[i] = left;
    }
    while (a.data.size() > 1 && a.data.back() == 0) {
        a.data.pop_back();
    }
    while (current.data.size() > 1 && current.data.back() == 0) {
        current.data.pop_back();
    }
    rem = current;
}

// Optimized version for division by uint64_t
void div_abs_inplace(bint_t& a, const uint64_t b, uint64_t& rem) {
    if (a.data.size() == 1 && a.data[0] < b) {
        rem = a.data[0];
        a = bint_t(0ll);
        return;
    }

    __uint128_t current = 0;
    for (int i = a.data.size() - 1; i >= 0; i--) {
        current = (current << 64) + a.data[i];
        if (current < b) {
            a.data[i] = 0ull;
            continue;
        }

        // Guess the bounds
        uint64_t left = current / (b + 1);
        uint64_t right = (current + 1) / b;

        const uint64_t initial_right = right;
        while (right - left > 1) {
            const uint64_t middle = left + (right - left) / 2;
            const __uint128_t guess = static_cast<__uint128_t>(b) * middle;
            if (guess <= current) {
                left = middle;
            } else {
                right = middle;
            }
        }
        if (right > left && right == initial_right) {
            // Check if right bound fits
            const __uint128_t guess = static_cast<__uint128_t>(b) * right;
            if (guess <= current) {
                current -= guess;
                a.data[i] = right;
                continue;
            }
        }
        const __uint128_t guess = static_cast<__uint128_t>(b) * left;
        current -= guess;
        a.data[i] = left;
    }
    while (a.data.size() > 1 && a.data.back() == 0) {
        a.data.pop_back();
    }
    rem = current;
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

bint_t mul(const bint_t& a, const bint_t& b, int a_limit, int b_limit) {
    if (a_limit < 0) a_limit = a.data.size();
    if (b_limit < 0) b_limit = b.data.size();
    auto result = mul_abs(a, b, a_limit, b_limit);
    result.sign = a.sign ^ b.sign;
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