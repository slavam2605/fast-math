#include <iostream>
#include "big_int.h"
#include "big_int_ops.h"
#include "karatsuba.h"

using namespace std;

bint_t get_mul(bool kar, int limit) {
    auto start = chrono::steady_clock::now();
    bint_t a(3ll);
    for (int i = 0; i < limit; i++) {
        a = kar ? karatsuba(a, a) : mul(a, a);
    }
    auto duration = chrono::steady_clock::now() - start;
    cout << chrono::duration_cast<chrono::milliseconds>(duration).count() << " ms" << endl;
    return a;
}

int main() {
    auto val = get_mul(true, 15);
    auto start = chrono::steady_clock::now();
    cout << val.to_string() << endl;
    auto duration = chrono::steady_clock::now() - start;
    cout << "to_string time: " << chrono::duration_cast<chrono::milliseconds>(duration).count() << " ms" << endl;

    if (false) {
        cout << "Time test: " << endl;
        auto kar_res = get_mul(true, 22);
        cout << endl;
    }
    if (false) {
        cout << "Equality test: " << endl;
        auto mul_res = get_mul(false, 21);
        auto kar_res = get_mul(true, 21);
        cout << (compare(mul_res, kar_res) == 0 ? "equal" : "different") << endl;
    }

    return 0;
}
