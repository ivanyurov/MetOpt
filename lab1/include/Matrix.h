#pragma once
#include <vector>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <cmath>

#include "Utils.h"

using std::vector;

class Matrix
{
public:
    std::vector<std::vector<double>> matrix;
    int n;
    int m;
    std::vector<int> line_indexes;
    std::vector<int> column_indexes;

    Matrix(const std::vector<std::vector<double>> &matrix);
    Matrix(const std::vector<double> &vect);

    Matrix &transpose();
    void print();
    void set_columns(std::vector<int> columns);
    void set_lines(std::vector<int> lines);
    Matrix &allocate_matrix(std::vector<int> lines, std::vector<int> columns);
    Matrix &multiply(const Matrix &a);
    Matrix &get_inverse_matrix();
    double determinant(std::vector<std::vector<double>> matrix);
    std::vector<int> get_addition_to_square_matrix(std::vector<int> &available_indexes);
    Matrix &subtract(Matrix &a);

    bool operator==(const Matrix &rhs);
    bool operator!=(const Matrix &rhs);
    bool operator<(const Matrix &rhs);
    bool operator>(double num);
    bool operator<(double num);
    bool operator<=(double num);
    bool operator>=(double num);
    Matrix &operator*(double num);
    int column_size();

    void gaussian_elimination();
    int compute_rank();
    bool is_full_rank();

private:
    std::vector<std::vector<double>> create_augmented_matrix(size_t size);
    void gauss_jordan_elimination(std::vector<std::vector<double>> &augmented);
};