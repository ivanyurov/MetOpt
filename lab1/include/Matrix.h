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
    // Храним базовую полную матрицу; line_indexes/column_indexes определяют
    // какие строки/столбцы рассматриваются как активная подматрица.
    std::vector<std::vector<double>> matrix;
    size_t n; // количество строк в полной матрице
    size_t m; // количество столбцов в полной матрице
    std::vector<int> line_indexes;
    std::vector<int> column_indexes;

    Matrix() = default;
    Matrix(const std::vector<std::vector<double>> &matrix);
    Matrix(const std::vector<double> &vect);

    // Методы возвращают объект Matrix по значению (без утечек и без "new")
    Matrix transpose() const;
    void print() const;
    void set_columns(const std::vector<int> &columns);
    void set_lines(const std::vector<int> &lines);
    Matrix allocate_matrix(const std::vector<int> &lines, const std::vector<int> &columns) const;
    Matrix multiply(const Matrix &a) const;
    Matrix get_inverse_matrix() const;
    double determinant(const std::vector<std::vector<double>> &mat) const;
    std::vector<int> get_addition_to_square_matrix(const std::vector<int> &available_indexes);
    Matrix subtract(const Matrix &a) const;

    // операции сравнения
    bool operator==(const Matrix &rhs) const;
    bool operator!=(const Matrix &rhs) const;
    bool operator<(const Matrix &rhs) const;

    // сравнения с числами (проверяют для всех элементов подматрицы)
    bool operator>(double num) const;
    bool operator<(double num) const;
    bool operator<=(double num) const;
    bool operator>=(double num) const;

    // умножение на скаляр — возвращает новый объект
    Matrix operator*(double num) const;
    int column_size() const;

    void gaussian_elimination(); // изменяет матрицу (использует line/column indexes)
    int compute_rank() const;
    bool is_full_rank() const;

private:
    void gauss_jordan_elimination(std::vector<std::vector<double>> &augmented) const;
};
