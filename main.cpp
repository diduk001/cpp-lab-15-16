#include "Matrix.h"
#include <iostream>

int main() {
    Matrix mtx1 = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            mtx1.set(i, j, i * 3 + j);
        }
    }

    Matrix mtx2 = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            mtx2.set(i, j, i * 3 - j);
        }
    }

    const auto mtx3 = mtx1.multiplyMultithread(mtx2);
    const auto mtx4 = mtx1.multiplyAsync(mtx2);
    std::cout << mtx3 << std::endl;
    std::cout << mtx4 << std::endl;
    std::cout << mtx3.isEqualMultithread(mtx4) << std::endl;
    std::cout << mtx3.isEqualAsync(mtx4) << std::endl;

    return 0;
}