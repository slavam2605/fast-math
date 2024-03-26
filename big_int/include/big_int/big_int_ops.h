//
// Created by Vyacheslav.Moklev on 09/03/2024.
//

#ifndef BIG_INT_OPS_H
#define BIG_INT_OPS_H

#include "big_int.h"

namespace big_int {
    constexpr int KARATSUBA_THRESHOLD = 30;
    constexpr int TOOM_COOK_THRESHOLD = 100;

    bint_t add(const bint_t& a, const bint_t& b);
    bint_t sub(const bint_t& a, const bint_t& b);
    bint_t multiply(const bint_t& a, const bint_t& b);

    void div_abs_inplace(bint_t& a, const bint_t& b, bint_t& rem);
    void fast_pow_inplace(bint_t& a, uint64_t n);
    void shift_left_inplace(bint_t& a, int shift);
    void shift_right_inplace(bint_t& a, int shift);

    std::strong_ordering compare(const bint_t& a, const bint_t& b);
}

namespace big_int_impl {
    bint_t schoolbook_multiply(const bint_t& a, const bint_t& b, int a_limit = -1, int b_limit = -1);
    void add_abs_inplace(bint_t& a, const bint_t& b, int b_limit = -1, int b_shift = 0);
    void sub_abs_inplace(bint_t& a, const bint_t& b, int b_limit = -1);
    void div_abs_inplace(bint_t& a, uint64_t b, uint64_t& rem);
    void normalize(bint_t& a);
}

#endif //BIG_INT_OPS_H
