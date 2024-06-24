//
// Created by Stepan Didurenko on 23.06.2024.
//

#ifndef CPP_LAB_15_16_MATRIX_H
#define CPP_LAB_15_16_MATRIX_H

#include <iostream>

// For comparing doubles
#define EPS 1e-9
#define THREADS 1024

class Matrix {
private:
    size_t rows;
    size_t cols;
    double **elements{};

public:
    // Constructors and destructor

    Matrix();

    Matrix(size_t newRows, size_t newCols);

    explicit Matrix(size_t rows) : Matrix(rows, rows) {}

    Matrix(size_t newRows, size_t newCols, double **newElements);

    Matrix(const Matrix &other);

    static Matrix identity(size_t rows);

    ~Matrix();

    void updateSize(size_t newRows, size_t newCols);

    // Getters, setters, indexing

    const double &get(size_t row, size_t col) const;

    void set(size_t row, size_t col, double newVal);

    size_t getRows() const {
        return rows;
    }

    size_t getCols() const {
        return cols;
    }

    bool isValidIndex(unsigned long long row, unsigned long long col) const {
        return row < this->getRows() &&
               col < this->getCols();
    }

    // Equality
    bool operator==(const Matrix &other) const;
    bool isEqualMultithread(const Matrix &other) const;
    bool isEqualAsync(const Matrix &other) const;

    bool operator!=(const Matrix &other) const;

    bool operator==(const double &k) const;
    bool isEqualMultithread(const double &k) const;
    bool isEqualAsync(const double &k) const;

    bool operator!=(const double &k) const;

    // Arithmetics

    Matrix operator*(const double &k) const;
    Matrix multiplyMultithread(const double &k) const;
    Matrix multiplyAsync(const double &k) const;

    Matrix operator-() const;

    Matrix operator+(const Matrix &other) const;
    Matrix addMultithread(const Matrix &other) const;
    Matrix addAsync(const Matrix &other) const;

    Matrix operator-(const Matrix &other) const;
    Matrix subtractMultithread(const Matrix &other) const;
    Matrix subtractAsync(const Matrix &other) const;

    Matrix operator*(const Matrix &other) const;
    Matrix multiplyMultithread(const Matrix &other) const;
    Matrix multiplyAsync(const Matrix &other) const;

    // Transpose matrix
    Matrix T() const;

    // Elementary operations
    void swapRows(size_t row1, size_t row2);

    void multiplyRow(size_t row, double k);

    void addRow(size_t row1, size_t row2, double k);

    // Find matrix determinant
    double det() const;

    // Inverse matrix
    Matrix operator!() const;
};

// Input and Output
std::istream &operator>>(std::istream &in, Matrix &mtx);

std::ostream &operator<<(std::ostream &out, const Matrix &mtx);

// LHS multiplication
Matrix operator*(const double &k, const Matrix &mtx);



#endif //CPP_LAB_15_16_MATRIX_H
