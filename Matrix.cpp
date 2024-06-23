//
// Created by Stepan Didurenko on 23.06.2024.
//

#include "Matrix.h"
#include "Matrix.h"
#include <stdexcept>

// Constructors and destructor

Matrix::Matrix() {
    this->rows = 0;
    this->cols = 0;
}

Matrix::Matrix(size_t newRows, size_t newCols) {
    if (newRows <= 0 || newCols <= 0)
        throw std::invalid_argument("Rows and newCols must be positive numbers");
    this->rows = newRows;
    this->cols = newCols;
    this->elements = new double *[newRows];
    for (int rowIdx = 0; rowIdx < newRows; rowIdx++) {
        this->elements[rowIdx] = new double[newCols];
        for (int colIdx = 0; colIdx < newCols; colIdx++)
            this->elements[rowIdx][colIdx] = 0;
    }
}

Matrix::Matrix(size_t newRows, size_t newCols, double **newElements) {
    if (newRows <= 0 || newCols <= 0)
        throw std::invalid_argument("Rows and newCols must be positive numbers");
    this->rows = newRows;
    this->cols = newCols;
    this->elements = new double *[newRows];
    for (int rowIdx = 0; rowIdx < newRows; rowIdx++) {
        this->elements[rowIdx] = new double[newCols];
        for (int colIdx = 0; colIdx < newCols; colIdx++)
            this->elements[rowIdx][colIdx] =
                    std::abs(newElements[rowIdx][colIdx]) < EPS ? 0 : newElements[rowIdx][colIdx];
    }
}

Matrix::Matrix(const Matrix &other) {
    this->rows = other.rows;
    this->cols = other.cols;
    this->elements = new double *[rows];
    for (int rowIdx = 0; rowIdx < rows; rowIdx++) {
        this->elements[rowIdx] = new double[cols];
        for (int colIdx = 0; colIdx < cols; colIdx++)
            this->elements[rowIdx][colIdx] = other.elements[rowIdx][colIdx];
    }
}

Matrix Matrix::identity(size_t rows) {
    Matrix result(rows);
    for (int i = 0; i < rows; i++)
        result.set(i, i, 1);
    return result;
}

Matrix::~Matrix() {
    for (int rowIdx = 0; rowIdx < rows; rowIdx++)
        delete[] elements[rowIdx];
    delete[] elements;
}

void Matrix::updateSize(size_t newRows, size_t newCols) {
    for (int rowIdx = 0; rowIdx < this->getRows(); rowIdx++)
        delete[] this->elements[rowIdx];
    delete[] this->elements;

    this->rows = newRows;
    this->cols = newCols;
    this->elements = new double *[this->rows];
    for (int rowIdx = 0; rowIdx < this->getRows(); rowIdx++) {
        this->elements[rowIdx] = new double[this->getCols()];
        for (int colIdx = 0; colIdx < this->getCols(); colIdx++)
            this->elements[rowIdx][colIdx] = 0;
    }
}

// Getters, setters, indexing

const double &Matrix::get(size_t row, size_t col) const {
    if (!this->isValidIndex(row, col))
        throw std::invalid_argument(
                "Maximal row is " + std::to_string(this->getRows() - 1) + ", " +
                "maximal col is " + std::to_string(this->getCols() - 1)
        );

    return this->elements[row][col];
}

void Matrix::set(size_t row, size_t col, double newVal) {
    if (!this->isValidIndex(row, col))
        throw std::invalid_argument(
                "Maximal row is " + std::to_string(this->getRows() - 1) + ", " +
                "maximal col is " + std::to_string(this->getCols() - 1)
        );

    this->elements[row][col] = std::abs(newVal) < EPS ? 0 : newVal;
}

// Equality

bool Matrix::operator==(const Matrix &other) const {
    if (this->getCols() != other.getCols() || this->getRows() != other.getRows()) return false;
    for (size_t rowIdx = 0; rowIdx < this->getRows(); rowIdx++)
        for (size_t colIdx = 0; colIdx < this->getCols(); colIdx++) {
            double thisValue = this->get(rowIdx, colIdx);
            double otherValue = other.get(rowIdx, colIdx);
            if (std::abs(thisValue - otherValue) > EPS)
                return false;
        }

    return true;
}

bool Matrix::operator!=(const Matrix &other) const {
    return !(*this == other);
}

bool Matrix::operator==(const double &k) const {
    return *this == k * Matrix::identity(this->getRows());
}

bool Matrix::operator!=(const double &k) const {
    return !(*this == k)
}


// Arithmetics

Matrix Matrix::operator*(const double &k) const {
    Matrix result(this->getRows(), this->getCols());
    for (size_t rowIdx = 0; rowIdx < this->getRows(); rowIdx++) {
        for (size_t colIdx = 0; colIdx < this->getCols(); colIdx++) {
            double thisValue = this->get(rowIdx, colIdx);
            result.set(
                    rowIdx,
                    colIdx,
                    thisValue * k
            );
        }
    }
    return result;
}

Matrix operator*(const double &k, const Matrix &mtx) {
    return mtx * k;
}

Matrix Matrix::operator-() const {
    return *this * -1;
}

Matrix Matrix::operator+(const Matrix &other) const {
    if (this->getCols() != other.getCols() || this->getRows() != other.getRows())
        throw std::invalid_argument("Matrices must have same dimensions");

    Matrix result(this->getRows(), this->getCols());
    for (size_t rowIdx = 0; rowIdx < this->getRows(); rowIdx++) {
        for (size_t colIdx = 0; colIdx < this->getCols(); colIdx++) {
            double thisValue = this->get(rowIdx, colIdx);
            double otherValue = other.get(rowIdx, colIdx);
            result.set(
                    rowIdx,
                    colIdx,
                    thisValue + otherValue
            );
        }
    }
    return result;
}

Matrix Matrix::operator-(const Matrix &other) const {
    return *this + (-other);
}

Matrix Matrix::operator*(const Matrix &other) const {
    if (this->getCols() != other.getRows())
        throw std::invalid_argument("Number of columns in first matrix must be equal to number of rows in second");

    Matrix result(this->getRows(), other.getCols());
    for (size_t rowIdx = 0; rowIdx < this->getRows(); rowIdx++) {
        for (size_t colIdx = 0; colIdx < other.getCols(); colIdx++) {
            double valueSum = 0;
            for (size_t i = 0; i < this->getCols(); i++) {
                double thisValue = this->get(rowIdx, i);
                double otherValue = other.get(i, colIdx);
                valueSum += thisValue * otherValue;
            }
            result.set(rowIdx, colIdx, valueSum);
        }
    }
    return result;
}

// Transpose matrix
Matrix Matrix::T() const {
    Matrix result(this->getCols(), this->getRows());
    for (size_t rowIdx = 0; rowIdx < this->getRows(); rowIdx++) {
        for (size_t colIdx = 0; colIdx < this->getCols(); colIdx++) {
            double thisValue = this->get(rowIdx, colIdx);
            result.set(colIdx, rowIdx, thisValue);
        }
    }
    return result;
}

// Elementary operations

void Matrix::swapRows(size_t row1, size_t row2) {
    if (!this->isValidIndex(row1, 0) || !this->isValidIndex(row2, 0))
        throw std::invalid_argument(
                "Maximal row is " + std::to_string(this->getRows() - 1)
        );

    std::swap(this->elements[row1], this->elements[row2]);
}

void Matrix::multiplyRow(size_t row, double k) {
    if (!this->isValidIndex(row, 0))
        throw std::invalid_argument(
                "Maximal row is " + std::to_string(this->getRows() - 1)
        );

    if (k == 0)
        throw std::invalid_argument("k must not be equal to 0");

    for (size_t colIdx = 0; colIdx < this->getCols(); colIdx++) {
        double thisValue = this->get(row, colIdx);
        this->set(row, colIdx, thisValue * k);
    }
}

void Matrix::addRow(size_t row1, size_t row2, double k) {
    if (!this->isValidIndex(row1, 0) || !this->isValidIndex(row2, 0))
        throw std::invalid_argument(
                "Maximal row is " + std::to_string(this->getRows() - 1)
        );

    if (row1 == row2)
        throw std::invalid_argument("Rows must be different");

    for (size_t colIdx = 0; colIdx < this->getCols(); colIdx++) {
        double row1Value = this->get(row1, colIdx);
        double row2Value = this->get(row2, colIdx);
        this->set(row1, colIdx, row1Value + row2Value * k);
    }
}

// Find matrix determinant
// Convert into diagonal matrix (like in Gauss method) and return product of all diagonal elements
// https://cp-algorithms.com/linear_algebra/determinant-gauss.html
double Matrix::det() const {
    if (this->getRows() != this->getCols())
        throw std::invalid_argument("Matrix must be square");

    const size_t size = this->getCols();
    if (size == 2) {
        double a = this->get(0, 0),
                b = this->get(0, 1),
                c = this->get(1, 0),
                d = this->get(1, 1);
        return a * d - b * c;
    }

    Matrix mtx(*this);
    double result = 1;
    for (size_t stepIdx = 0; stepIdx < size; stepIdx++) {
        // Choose row with maximal absolute value of stepIdx_th element
        size_t maxByAbsRowIdx = stepIdx;
        double maxByAbsValue = mtx.get(maxByAbsRowIdx, stepIdx);
        for (size_t maxByAbsRowCandidate = stepIdx + 1; maxByAbsRowCandidate < size; maxByAbsRowCandidate++) {
            double candidateValue = std::abs(mtx.get(maxByAbsRowCandidate, stepIdx));
            if (candidateValue > std::abs(maxByAbsValue)) {
                maxByAbsRowIdx = maxByAbsRowCandidate;
                maxByAbsValue = mtx.get(maxByAbsRowIdx, stepIdx);
            }
        }
        if (std::abs(maxByAbsValue) < EPS)
            return 0;

        // Swap current step row and max abs row
        if (stepIdx != maxByAbsRowIdx)
            result *= -1;
        mtx.swapRows(stepIdx, maxByAbsRowIdx);

        // Divide max abs row by first stepIdx element
        result *= maxByAbsValue;
        mtx.multiplyRow(stepIdx, 1 / maxByAbsValue);

        // Nullify all elements under [stepIdx][stepIdx]
        for (size_t lowerRowIdx = stepIdx + 1; lowerRowIdx < size; lowerRowIdx++) {
            if (std::abs(mtx.get(lowerRowIdx, stepIdx)) < EPS)
                continue;
            double coefficient = -mtx.get(lowerRowIdx, stepIdx);
            mtx.addRow(lowerRowIdx, stepIdx, coefficient);
        }
    }

    return result;
}

// Inverse matrix
// Convert into identity matrix (like in Gauss method) and do same operation on identity matrix
// https://cp-algorithms.com/linear_algebra/linear-system-gauss.html
// Returns 0 matrix if matrix is not invertible
Matrix Matrix::operator!() const {
    if (this->getRows() != this->getCols())
        throw std::invalid_argument("Matrix must be square");

    const size_t size = this->getCols();
    const double det = this->det();
    if (det == 0)
        return {size, size};

    if (size == 2) {
        double a = this->get(0, 0),
                b = this->get(0, 1),
                c = this->get(1, 0),
                d = this->get(1, 1);

        double **resultElements;
        resultElements = new double *[2];
        resultElements[0] = new double[2]{d, -b};
        resultElements[1] = new double[2]{-c, a};

        return (1 / det) * Matrix(2, 2, elements);
    }

    Matrix mtx = Matrix(*this);
    Matrix result = Matrix::identity(size);

    for (size_t stepIdx = 0; stepIdx < size; stepIdx++) {
        // Choose row with maximal absolute value of stepIdx_th element
        size_t maxByAbsRowIdx = stepIdx;
        double maxByAbsValue = mtx.get(maxByAbsRowIdx, stepIdx);
        for (size_t maxByAbsRowCandidate = stepIdx + 1; maxByAbsRowCandidate < size; maxByAbsRowCandidate++) {
            double candidateValue = std::abs(mtx.get(maxByAbsRowCandidate, stepIdx));
            if (candidateValue > std::abs(maxByAbsValue)) {
                maxByAbsRowIdx = maxByAbsRowCandidate;
                maxByAbsValue = mtx.get(maxByAbsRowIdx, stepIdx);
            }
        }

        // Swap current step row and max abs row
        mtx.swapRows(stepIdx, maxByAbsRowIdx);
        result.swapRows(stepIdx, maxByAbsRowIdx);

        // Divide max abs row by first stepIdx element
        mtx.multiplyRow(stepIdx, 1 / maxByAbsValue);
        result.multiplyRow(stepIdx, 1 / maxByAbsValue);

        // Nullify all elements under [stepIdx][stepIdx]
        for (size_t lowerRowIdx = 0; lowerRowIdx < size; lowerRowIdx++) {
            if (lowerRowIdx == stepIdx || std::abs(mtx.get(lowerRowIdx, stepIdx)) < EPS) {
                continue;
            }
            double coefficient = -mtx.get(lowerRowIdx, stepIdx);
            mtx.addRow(lowerRowIdx, stepIdx, coefficient);
            result.addRow(lowerRowIdx, stepIdx, coefficient);
        }
    }

    return result;
}

// Input and Output

std::istream &operator>>(std::istream &in, Matrix &mtx) {
    size_t newRows, newCols;
    in >> newRows >> newCols;
    mtx.updateSize(newRows, newCols);
    for (size_t rowIdx = 0; rowIdx < mtx.getRows(); rowIdx++)
        for (size_t colIdx = 0; colIdx < mtx.getCols(); colIdx++) {
            double value;
            in >> value;
            mtx.set(rowIdx, colIdx, value);
        }
    return in;
}

std::ostream &operator<<(std::ostream &out, const Matrix &mtx) {
    const size_t rows = mtx.getRows();
    const size_t cols = mtx.getCols();

    out << "[" << std::endl;
    for (size_t rowIdx = 0; rowIdx < rows; rowIdx++) {
        out << "\t[";
        for (size_t colIdx = 0; colIdx < cols; colIdx++) {
            out << "\t" << mtx.get(rowIdx, colIdx);
            if (colIdx < cols - 1)
                out << ",";
        }
        out << "]";
        if (rowIdx < rows - 1)
            out << ",";
        out << std::endl;
    }
    out << "]" << std::endl;

    return out;
}