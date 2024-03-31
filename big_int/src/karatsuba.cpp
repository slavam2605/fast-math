//
// Created by Vyacheslav.Moklev on 09/03/2024.
//

#include "big_int/karatsuba.h"
#include "big_int/big_int_ops.h"

using namespace big_int_impl;

std::pair<bint_t, int> split_abs(const bint_t& a, int m, int a_limit) {
    if (a_limit <= m)
        return std::make_pair(bint_t(0ll), a_limit);

    bint_t high;
    high.data.resize(a_limit - m);
    std::copy_n(a.data.begin() + m, a_limit - m, high.data.begin());
    return std::make_pair(high, m);
}

bint_t split_abs_inplace(bint_t& a, const int m, int a_limit) {
    if (a_limit <= m) {
        throw std::runtime_error("split_abs_inplace: a_limit [" + std::to_string(a_limit) +
            "] <= m [" + std::to_string(2 * m) + "]");
    }

    bint_t high;
    high.data.resize(a_limit - m);
    std::copy_n(a.data.begin() + m, a_limit - m, high.data.begin());
    normalize(high);

    a.sign = false;
    a.data.resize(m);
    normalize(a);
    return high;
}

bint_t karatsuba(const bint_t& a, const bint_t& b, int a_limit, int b_limit) { // NOLINT(*-no-recursion)
    if (a_limit < 0) a_limit = a.data.size();
    if (b_limit < 0) b_limit = b.data.size();

    const int m = (std::max(a_limit, b_limit) + 1) / 2;
    auto [a1, a0_limit] = split_abs(a, m, a_limit);
    auto [b1, b0_limit] = split_abs(b, m, b_limit);

    auto z2 = big_int::multiply(a1, b1);
    auto z0 = big_int::multiply(a, b, a0_limit, b0_limit);
    add_abs_inplace(a1, a, 0, a0_limit);
    add_abs_inplace(b1, b, 0, b0_limit);
    auto z1 = big_int::multiply(a1, b1);
    sub_abs_inplace(z1, z2);
    sub_abs_inplace(z1, z0);

    add_abs_inplace(z0, z1, 0, -1, m);
    add_abs_inplace(z0, z2, 0, -1, 2 * m);
    z0.sign = a.sign ^ b.sign;
    normalize(z0);

    return z0;
}

void karatsuba_square(bint_t& a, int a_limit) { // NOLINT(*-no-recursion)
    if (a_limit < 0) a_limit = a.data.size();

    const int m = (a_limit + 1) / 2;
    auto a1 = split_abs_inplace(a, m, a_limit);

    auto z2 = a1;
    square(z2);
    add_abs_inplace(a1, a);
    square(a);
    square(a1);
    sub_abs_inplace(a1, z2);
    sub_abs_inplace(a1, a);

    add_abs_inplace(a, a1, 0, -1, m);
    add_abs_inplace(a, z2, 0, -1, 2 * m);
    normalize(a);
}
