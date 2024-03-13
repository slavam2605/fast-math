//
// Created by Vyacheslav.Moklev on 09/03/2024.
//

#ifndef BIG_INT_OPS_H
#define BIG_INT_OPS_H

#include "big_int.h"

bint_t add(const bint_t& a, const bint_t& b);
bint_t sub(const bint_t& a, const bint_t& b);
bint_t slow_mul(const bint_t& a, const bint_t& b, int a_limit = -1, int b_limit = -1);

std::strong_ordering compare(const bint_t& a, const bint_t& b);

void shift_left_inplace(bint_t& a, int shift);
void shift_right_inplace(bint_t& a, int shift);
void fast_pow_inplace(bint_t& a, uint64_t n);
void add_abs_inplace(bint_t& a, const bint_t& b, int b_limit = -1, int b_shift = 0);
void sub_abs_inplace(bint_t& a, const bint_t& b, int b_limit = -1);
bint_t mul_abs_uint64(const bint_t& a, uint64_t b);
void div_abs_inplace(bint_t& a, const bint_t& b, bint_t& rem);
void div_abs_inplace(bint_t& a, uint64_t b, uint64_t& rem);
void debug_print(const bint_t& a);

#endif //BIG_INT_OPS_H
