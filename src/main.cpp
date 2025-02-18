#include "matrix.h"
#include <iostream>

int main() {
    try {
        int n;
        std::cout << "Enter the number of equations: ";
        std::cin >> n;

        Matrix A(n, n);
        std::cout << "Enter coefficients of matrix A:\n";
        A.fill();

        std::vector<double> b(n);
        std::cout << "Enter the constants vector b:\n";
        for (int i = 0; i < n; ++i) {
            std::cin >> b[i];
        }

        std::vector<double> solution = A.gaussJordanSolve(b);

        std::cout << "Solution:\n";
        for (int i = 0; i < n; ++i) {
            std::cout << "x" << i + 1 << " = " << solution[i] << std::endl;
        }
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
