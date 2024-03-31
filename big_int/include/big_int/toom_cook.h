//
// Created by Vyacheslav.Moklev on 24/03/2024.
//

#ifndef TOOM_COOK_H
#define TOOM_COOK_H

#include "big_int.h"

bint_t toom3(const bint_t& a, const bint_t& b, int a_limit = -1, int b_limit = -1);
void toom3_square(bint_t& a, int a_limit = -1);

#endif //TOOM_COOK_H
