#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

class Matrix {
private:
    int rows, cols;
    std::vector<std::vector<double>> data;

public:
    Matrix(int r, int c);

    void fill();

    void printMatrix();

    void set(double element, int r, int c);

    void printElement(int r, int c);

    std::vector<double> gaussJordanSolve(std::vector<double> b);
};

#endif
