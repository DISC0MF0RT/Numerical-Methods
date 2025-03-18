#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
class Matrix {
private:
    int rows, cols;
    std::vector<std::vector<double>> data;

    double determinantRecursive(const std::vector<std::vector<double>>& mat) const;
public:
    Matrix(int r, int c);

    void fill();

    void printMatrix();

    void set(double element, int r, int c);

    void printElement(int r, int c);

    bool isTriDiagonal();

    double determinant() const;

    std::vector<std::vector<double>> inverse() const;

    std::vector<double> gaussJordanSolve(std::vector<double> b);

    std::vector<double> progon(std::vector<double> b);

    bool check(std::vector<double> x, std::vector<double> b);
    bool check(std::vector<double> x, std::vector<double> b, double EPS);

    std::vector<double> methodSimpleIteration(const std::vector<double> &b, double EPS, int maxIteration);

    std::vector<double> methodZeidelua(const std::vector<double> &b, double EPS, int maxIteration);


};

#endif
