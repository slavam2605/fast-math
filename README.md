# fast-math: Library for big integers

fast-math is a C++ library for fast computing with big integers.
It uses base-2<sup>64</sup> number representation and `__uint128_t` GCC-extension
for faster operations.

## Features

* Basic operations: addition, subtraction, shifting, comparison.
* Fast multiplication with Karatsuba and Toom-Cook (Toom-3) algorithms, O(n<sup>log<sub>3</sub>5</sup>) ≈ O(n<sup>1.46</sup>) time complexity.
* Fast division with [Burnikel-Ziegler algorithm](https://pure.mpg.de/rest/items/item_1819444_4/component/file_2599480/content): O(n<sup>log<sub>3</sub>5</sup> log n) ≈ O(n<sup>1.46</sup> log n) time complexity with Toom-3 multiplication.
* Fast divide-and-conquer `to_string` implementation

## Installation
Clone the repository, build the project and install using CMake:
```bash
git clone https://github.com/slavam2605/fast-math
cd fast-math && mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
sudo make install
```
After that, just add these lines to CMakeLists.txt in your project:
```cmake
find_package(BigIntLib REQUIRED)
target_link_libraries(<YourTargetName> BigIntLib::BigIntLib)
```

## Usage
An example of computing a factorial of 1000:
```c++
#include <iostream>
#include <big_int/big_int.h>

int main() {
    bint_t factorial(1);
    for (int i = 2; i <= 1000; i++) {
        factorial = factorial * bint_t(i);
    }

    // Prints 402387260077093773543702433923003985...(2532 more digits)
    std::cout << factorial.to_string() << std::endl;
    return 0;
}
```