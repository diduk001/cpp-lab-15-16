#include "Matrix.h"
#include "Timer.h"
#include <iostream>
#include <random>

#define MATRIX_SIZE 500

int main() {
    std::mt19937 rng;
    rng.seed(1234);
    std::uniform_real_distribution<double> distribution(0, 1);

    Matrix mtx1 = Matrix(MATRIX_SIZE, MATRIX_SIZE);
    Matrix mtx2 = Matrix(MATRIX_SIZE, MATRIX_SIZE);

    for (size_t row_idx = 0; row_idx < mtx1.getRows(); row_idx++) {
        for (size_t col_idx = 0; col_idx < mtx1.getCols(); col_idx++) {
            double value1 = distribution(rng);
            double value2 = distribution(rng);
            mtx1.set(row_idx, col_idx, value1);
            mtx2.set(row_idx, col_idx, value2);
        }
    }

    Timer timer = Timer();
    timer.start();

    Matrix mtx3 = mtx1.multiplyMultithread(mtx2);

    const auto duration = timer.get();
    std::cout << "Matrix Size: " << MATRIX_SIZE << std::endl;
    std::cout << "Threads: " << THREADS << std::endl;
    std::cout << "Duration: " << duration.count() << " ms" << std::endl;

    return 0;
}