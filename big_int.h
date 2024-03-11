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
    explicit bint_t(uint64_t value);
    explicit bint_t(int64_t value);

    bint_t& operator=(const bint_t& other) = default;
    std::strong_ordering operator<=>(const bint_t& other) const;
    bool operator==(const bint_t& other) const;
    bint_t operator+(const bint_t& other) const;
    bint_t operator-(const bint_t& other) const;
    bint_t operator*(const bint_t& other) const;

    [[nodiscard]] std::string to_string() const;
};

#endif //BIG_INT_H
