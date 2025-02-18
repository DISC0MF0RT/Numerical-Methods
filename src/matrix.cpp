#include "matrix.h"

Matrix::Matrix(int r, int c) : rows(r), cols(c), data(r, std::vector<double>(c, 0.0)) {}

void Matrix::fill() {
    std::cout << "Enter matrix elements:\n";
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cin >> data[i][j];
        }
    }
}

void Matrix::printMatrix() {
    for (const auto &row : data) {
        for (double val : row) {
            std::cout << val << "\t";
        }
        std::cout << std::endl;
    }
}

void Matrix::set(double element, int r, int c) {
    if (r >= rows || c >= cols || r < 0 || c < 0) {
        throw std::out_of_range("Index out of range");
    }
    data[r][c] = element;
}

void Matrix::printElement(int r, int c) {
    if (r >= rows || c >= cols || r < 0 || c < 0) {
        std::cerr << "Index out of range" << std::endl;
        return;
    }
    std::cout << "Element at (" << r << ", " << c << "): " << data[r][c] << std::endl;
}

std::vector<double> Matrix::gaussJordanSolve(std::vector<double> b) {
    if (b.size() != rows) {
        throw std::invalid_argument("Vector b size must match the number of rows in the matrix");
    }

    std::vector<std::vector<double>> extended = data;
    for (int i = 0; i < rows; ++i) {
        extended[i].push_back(b[i]);
    }

    for (int i = 0; i < rows; ++i) {
        int maxRow = i;
        for (int j = i + 1; j < rows; ++j) {
            if (abs(extended[j][i]) > abs(extended[maxRow][i])) {
                maxRow = j;
            }
        }

        if (maxRow != i) {
            std::swap(extended[i], extended[maxRow]);
        }

        double diag = extended[i][i];
        if (abs(diag) < 1e-12) {
            throw std::runtime_error("Matrix is singular or system has no unique solution");
        }

        for (int j = i; j <= cols; ++j) {
            extended[i][j] /= diag;
        }

        for (int j = 0; j < rows; ++j) {
            if (j != i) {
                double factor = extended[j][i];
                for (int k = i; k <= cols; ++k) {
                    extended[j][k] -= factor * extended[i][k];
                }
            }
        }
    }

    std::vector<double> x(rows);
    for (int i = 0; i < rows; ++i) {
        x[i] = extended[i][cols];
    }

    return x;
}
