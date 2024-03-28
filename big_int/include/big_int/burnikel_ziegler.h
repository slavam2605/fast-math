//
// Created by Vyacheslav.Moklev on 27/03/2024.
//

#ifndef BURNIKEL_ZIEGLER_H
#define BURNIKEL_ZIEGLER_H

#include <big_int/big_int.h>

constexpr int BURNIKEL_ZIEGLER_THRESHOLD = 100;

void divide_burnikel_ziegler(bint_t& a, const bint_t& b, bint_t& rem);

#endif //BURNIKEL_ZIEGLER_H
