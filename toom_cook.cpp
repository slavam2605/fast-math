//
// Created by Vyacheslav.Moklev on 24/03/2024.
//

#include "toom_cook.h"

#include <iostream>

#include "big_int_ops.h"

std::tuple<bint_t, bint_t, bint_t> split_abs(const bint_t& a, int m) {
    if (a.data.size() <= m) {
        auto copy = a;
        copy.sign = false;
        return std::make_tuple(bint_t(0ll), bint_t(0ll), copy);
    }

    bint_t low;
    low.data.resize(m);
    std::copy_n(a.data.begin(), m, low.data.begin());

    bint_t mid;
    mid.data.resize(std::min(m, static_cast<int>(a.data.size()) - m));
    std::copy_n(a.data.begin() + m, mid.data.size(), mid.data.begin());

    if (a.data.size() <= 2 * m)
        return std::make_tuple(bint_t(0ll), mid, low);

    bint_t high;
    high.data.resize(std::min(m, static_cast<int>(a.data.size()) - 2 * m));
    std::copy_n(a.data.begin() + 2 * m, high.data.size(), high.data.begin());

    return std::make_tuple(high, mid, low);
}

void div_abs_inplace_3_assert_rem(bint_t& a) {
    __uint128_t current = 0;
    for (int i = a.data.size() - 1; i >= 0; i--) {
        current = (current << 64) + a.data[i];
        a.data[i] = current / 3;
        current %= 3;
    }
    while (a.data.size() > 1 && a.data.back() == 0) {
        a.data.pop_back();
    }
    if (current != 0) throw std::runtime_error("div_abs_inplace_3_assert_rem: rem was " +
        std::to_string(static_cast<uint64_t>(current)) + ", a = " + a.to_string());
}

bint_t toom3(const bint_t& a, const bint_t& b) {
    int m = (std::max(a.data.size(), b.data.size()) + 2) / 3;
    auto [a2, a1, a0] = split_abs(a, m);
    auto [b2, b1, b0] = split_abs(b, m);

    auto p0 = a0 + a2;
    auto p_1 = p0 + a1;
    auto p_m1 = p0 - a1;
    auto p_m2 = p_m1 + a2;
    shift_left_inplace(p_m2, 1);
    p_m2 = p_m2 - a0;

    auto q0 = b0 + b2;
    auto q_1 = q0 + b1;
    auto q_m1 = q0 - b1;
    auto q_m2 = q_m1 + b2;
    shift_left_inplace(q_m2, 1);
    q_m2 = q_m2 - b0;

    auto r_0 = multiply(a0, b0);
    auto r_1 = multiply(p_1, q_1);
    auto r_m1 = multiply(p_m1, q_m1);
    auto r_m2 = multiply(p_m2, q_m2);
    auto r_inf = multiply(a2, b2);

    auto r4 = r_inf;
    auto r3 = r_m2 - r_1;
    div_abs_inplace_3_assert_rem(r3);
    auto r1 = r_1 - r_m1;
    shift_right_inplace(r1, 1);
    auto r2 = r_m1 - r_0;
    auto diff = r2 - r3;
    r3 = r2 - r3;
    shift_right_inplace(r3, 1);
    shift_left_inplace(r_inf, 1);
    r3 = r3 + r_inf;
    r2 = r2 + r1 - r4;
    r1 = r1 - r3;

    add_abs_inplace(r_0, r1, -1, m);
    add_abs_inplace(r_0, r2, -1, 2 * m);
    add_abs_inplace(r_0, r3, -1, 3 * m);
    add_abs_inplace(r_0, r4, -1, 4 * m);
    r_0.sign = a.sign ^ b.sign;
    r_0.normalize();
    return r_0;
}