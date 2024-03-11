//
// Created by Vyacheslav.Moklev on 09/03/2024.
//

#ifndef KATATSUBA_H
#define KATATSUBA_H

#include "big_int.h"

// extern int SPLIT_LIMIT;

bint_t karatsuba(const bint_t& a, const bint_t& b, int a_limit = -1, int b_limit = -1);

#endif //KATATSUBA_H
