//
// Created by Vyacheslav.Moklev on 09/03/2024.
//

#ifndef BIG_INT_H
#define BIG_INT_H

#include <vector>
#include <string>

struct bint_t {
    bool sign;
    std::vector<uint64_t> data;

    bint_t();
    bint_t(const bint_t& other) = default;
    bint_t(bint_t&& other) noexcept;
    explicit bint_t(uint64_t value);
    explicit bint_t(int64_t value);
    explicit bint_t(int value);

    bint_t& pow(int n);
    [[nodiscard]] bint_t pow(int n) const;

    bint_t& operator=(const bint_t& other) = default;
    bint_t& operator=(bint_t&& other) = default;
    std::strong_ordering operator<=>(const bint_t& other) const;
    bool operator==(const bint_t& other) const;
    bint_t operator-() const;

    bint_t operator+(const bint_t& other) const;
    bint_t operator-(const bint_t& other) const;
    bint_t operator*(const bint_t& other) const;
    bint_t operator/(const bint_t& other) const;
    bint_t operator%(const bint_t& other) const;
    bint_t operator<<(int64_t n) const;
    bint_t operator>>(int64_t n) const;

    bint_t& operator+=(const bint_t& other);
    bint_t& operator-=(const bint_t& other);
    bint_t& operator*=(const bint_t& other);
    bint_t& operator/=(const bint_t& other);
    bint_t& operator%=(const bint_t& other);
    bint_t& operator<<=(int64_t n);
    bint_t& operator>>=(int64_t n);

    [[nodiscard]] std::string to_string() const;
};

#endif //BIG_INT_H
