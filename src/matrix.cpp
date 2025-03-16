#include "matrix.h"

Matrix::Matrix(int r, int c) : rows(r), cols(c), data(r, std::vector<double>(c, 0.0)) {}

void Matrix::fill() {
    std::cout << "Enter matrix elements:\n";
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << "i = " << i << " j = " << j << ": ";
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

bool Matrix::isTriDiagonal(){
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++) {
            if (std::abs(i - j) > 1 && data[i][j] != 0){
                return false;
            }
        }
    }
    return true;
}

double Matrix::determinant() const {    
    if (rows != cols) {
        throw std::invalid_argument("Determinant is only defined for square matrices");
    }
    return determinantRecursive(data);
}

double Matrix::determinantRecursive(const std::vector<std::vector<double>>& mat) const {
    int n = mat.size();
    if (n == 1) {
        return mat[0][0];
    }
    if (n == 2) {
        return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    }

    double det = 0.0;
    for (int i = 0; i < n; ++i) {
        std::vector<std::vector<double>> subMatrix(n - 1, std::vector<double>(n - 1));
        for (int row = 1; row < n; ++row) {
            int colIdx = 0;
            for (int col = 0; col < n; ++col) {
                if (col == i) continue;
                subMatrix[row - 1][colIdx] = mat[row][col];
                ++colIdx;
            }
        }
        det += (i % 2 == 0 ? 1 : -1) * mat[0][i] * determinantRecursive(subMatrix);
    }
    return det;
}

std::vector<std::vector<double>> Matrix::inverse() const {
    if (rows != cols) {
        throw std::invalid_argument("Matrix must be square to find inverse");
    }
    if (std::abs(determinant()) < 1e-12) {
        throw std::runtime_error("Matrix is singular and cannot be inverted");
    }

    std::vector<std::vector<double>> extended = data;
    for (int i = 0; i < rows; ++i) {
        extended[i].resize(2 * cols, 0.0);
        extended[i][cols + i] = 1.0;
    }

    for (int i = 0; i < rows; ++i) {
        int maxRow = i;
        for (int j = i + 1; j < rows; ++j) {
            if (std::abs(extended[j][i]) > std::abs(extended[maxRow][i])) {
                maxRow = j;
            }
        }

        if (maxRow != i) {
            std::swap(extended[i], extended[maxRow]);
        }

        double diag = extended[i][i];
        if (abs(diag) < 1e-12) {
            throw std::runtime_error("Matrix is singular and cannot be inverted");
        }

        for (int j = 0; j < 2 * cols; ++j) {
            extended[i][j] /= diag;
        }

        for (int j = 0; j < rows; ++j) {
            if (j != i) {
                double factor = extended[j][i];
                for (int k = 0; k < 2 * cols; ++k) {
                    extended[j][k] -= factor * extended[i][k];
                }
            }
        }
    }

    std::vector<std::vector<double>> inverse(rows, std::vector<double>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            inverse[i][j] = extended[i][cols + j];
        }
    }
    return inverse;
}

std::vector<double> Matrix::gaussJordanSolve(std::vector<double> b) {
    if (b.size() != rows) {
        throw std::invalid_argument("Vector b size must match the number of rows in the matrix");
    }
    if (std::abs(determinant()) < 1e-12){
         throw std::runtime_error("Matrix is singular or system has no unique solution");
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

std::vector<double> Matrix::progon(std::vector<double> b){
    if (rows != b.size()){
        throw std::invalid_argument("Vector b size must match the number of rows in the matrix");
    }
    if (!isTriDiagonal()){
        throw std::invalid_argument("Not tridiagonal");
    }
    if(std::abs(determinant()) < 1e-12){
        throw std::runtime_error("Matrix is singular or system has no unique solution");

    }

    int n = rows;
    std::vector<double> alpha(n, 0.0), beta(n, 0.0), x(n, 0.0);

    alpha[0] = -data[0][1] / data[0][0];
    beta[0] = b[0] / data[0][0];

    for (int i = 1; i < n - 1; i++){
        double denom = data[i][i] + data[i][i-1] * alpha[i -1];
        alpha[i] = -data[i][i + 1] / denom;
        beta[i] = (b[i] - data[i][i - 1] * beta[i - 1]) / denom;
    }

    beta[n - 1] = (b[n - 1] - data[n - 1][n - 2] * beta[n - 2]) / 
    (data[n - 1][n - 1] + data[n - 1][n - 2] * alpha[n - 2]);

    x[n - 1] = beta[n - 1];
    for (int i = n - 2; i >= 0; i--){
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }

    return x;
}

bool Matrix::check(std::vector<double> x, std::vector<double> b){
    std::vector<double> calculated_b(rows, 0.0);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            calculated_b[i] += data[i][j] * x[j];
        }
    }

    for (int i = 0; i < rows; ++i) {
        if (std::abs(calculated_b[i] - b[i]) > 1e-9) {
            return false;
        }
    }
    
    return true;
}
