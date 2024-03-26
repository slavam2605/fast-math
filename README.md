# fast-math: Library for big integers

fast-math is a C++ library for fast computing with big integers.
It uses base-2<sup>64</sup> number representation and `__uint128_t` GCC-extension
for faster operations.

## Features

* Basic operations: addition, subtraction, shifting, comparison.
* Fast multiplication with Karatsuba and Toom-Cook (Toom-3) algorithms, O(n<sup>log<sub>3</sub>5</sup>) ≈ O(n<sup>1.46</sup>) time complexity.
* Division with Knuth's algorithm D ideas: still O(n<sup>2</sup>) time complexity, but faster than trivial implementation.
* Fast power, using O(log n) multiplications.
* Fast divide-and-conquer `to_string` implementation.