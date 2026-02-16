#include "../include/Matrix.h"

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <cmath>

using namespace std;

// Конструкторы
Matrix::Matrix(const vector<vector<double>> &matrix)
{
    if (matrix.empty()) {
        this->matrix = matrix;
        n = 0;
        m = 0;
    } else {
        this->matrix = matrix;
        n = matrix.size();
        m = matrix[0].size();
    }
    line_indexes.clear();
    column_indexes.clear();
    for (size_t i = 0; i < n; ++i) line_indexes.push_back((int)i);
    for (size_t j = 0; j < m; ++j) column_indexes.push_back((int)j);
}

Matrix::Matrix(const vector<double> &vect)
{
    matrix.clear();
    for (double el : vect)
        matrix.push_back(vector<double>(1, el));

    n = matrix.size();
    m = (n == 0 ? 0 : matrix[0].size());
    line_indexes.clear();
    column_indexes.clear();
    for (size_t i = 0; i < n; ++i) line_indexes.push_back((int)i);
    for (size_t j = 0; j < m; ++j) column_indexes.push_back((int)j);
}

// Транспонирование (возвращаем новую матрицу)
Matrix Matrix::transpose() const
{
    if (n == 0 || m == 0) return Matrix(std::vector<std::vector<double>>{});
    vector<vector<double>> t_matrix(m, vector<double>(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            t_matrix[j][i] = matrix[i][j];

    Matrix ret(t_matrix);
    std::vector<int> new_columns(line_indexes.begin(), line_indexes.end());
    std::vector<int> new_lines(column_indexes.begin(), column_indexes.end());
    ret.set_columns(new_columns);
    ret.set_lines(new_lines);
    return ret;
}

// Вывод матрицы (только активная подматрица)
void Matrix::print() const
{
    cout << endl;
    for (size_t i = 0; i < line_indexes.size(); ++i)
    {
        for (size_t j = 0; j < column_indexes.size(); ++j)
        {
            cout << matrix[line_indexes[i]][column_indexes[j]] << " ";
        }
        cout << std::endl;
    }
    cout << "\n";
}

// Установка индексов столбцов и строк
void Matrix::set_columns(const std::vector<int> &columns)
{
    column_indexes = columns;
}

void Matrix::set_lines(const std::vector<int> &lines)
{
    line_indexes = lines;
}

// Выделение подматрицы (возвращаем новую матрицу)
Matrix Matrix::allocate_matrix(const std::vector<int> &lines, const std::vector<int> &columns) const
{
    if (lines.empty() || columns.empty())
        return Matrix(std::vector<std::vector<double>>{});

    vector<vector<double>> new_matrix(lines.size(), vector<double>(columns.size(), 0.0));
    for (size_t i = 0; i < lines.size(); ++i)
    {
        for (size_t j = 0; j < columns.size(); ++j)
        {
            int li = lines[i];
            int cj = columns[j];
            if (li < 0 || (size_t)li >= (size_t)matrix.size() || cj < 0 || (size_t)cj >= (size_t)matrix[0].size())
                throw std::out_of_range("allocate_matrix: index out of range");
            new_matrix[i][j] = matrix[li][cj];
        }
    }
    Matrix ret(new_matrix);
    return ret;
}

// Умножение активной подматрицы this (rows=line_indexes, cols=column_indexes)
// на активную подматрицу a (rows=a.line_indexes, cols=a.column_indexes).
Matrix Matrix::multiply(const Matrix &a) const
{
    size_t lhs_rows = line_indexes.size();
    size_t lhs_cols = column_indexes.size();
    size_t rhs_rows = a.line_indexes.size();
    size_t rhs_cols = a.column_indexes.size();

    if (lhs_cols != rhs_rows)
        throw std::invalid_argument("multiply: incompatible dimensions");

    vector<vector<double>> result(lhs_rows, vector<double>(rhs_cols, 0.0));

    for (size_t i = 0; i < lhs_rows; ++i)
    {
        for (size_t k = 0; k < rhs_cols; ++k)
        {
            double sum = 0.0;
            for (size_t j = 0; j < lhs_cols; ++j)
            {
                double left_val = matrix[line_indexes[i]][column_indexes[j]];
                double right_val = a.matrix[a.line_indexes[j]][a.column_indexes[k]];
                sum += left_val * right_val;
            }
            result[i][k] = sum;
        }
    }

    Matrix ans(result);
    return ans;
}

Matrix Matrix::get_inverse_matrix() const
{
    size_t sz = line_indexes.size();
    if (sz == 0 || sz != column_indexes.size())
    {
        throw std::invalid_argument("Matrix must be square (active submatrix)");
    }

    vector<vector<double>> augmented(sz, vector<double>(2 * sz, 0.0));
    for (size_t i = 0; i < sz; ++i)
    {
        for (size_t j = 0; j < sz; ++j)
        {
            augmented[i][j] = matrix[line_indexes[i]][column_indexes[j]];
        }
        augmented[i][i + sz] = 1.0;
    }

    for (size_t i = 0; i < sz; ++i)
    {
        // Правильный проход по строкам (а не по column_indexes)
        size_t max_row = i;
        for (size_t k = i + 1; k < sz; ++k)
        {
            if (std::fabs(augmented[k][i]) > std::fabs(augmented[max_row][i]))
            {
                max_row = k;
            }
        }

        if (max_row != i)
        {
            std::swap(augmented[i], augmented[max_row]);
        }

        if (std::fabs(augmented[i][i]) < 1e-14)
        {
            throw std::invalid_argument("Matrix is singular and cannot be inverted");
        }

        double pivot = augmented[i][i];
        for (size_t j = 0; j < 2 * sz; ++j)
        {
            augmented[i][j] = (augmented[i][j] / pivot);
        }

        for (size_t k = 0; k < sz; ++k)
        {
            if (k != i)
            {
                double factor = augmented[k][i];
                for (size_t j = 0; j < 2 * sz; ++j)
                {
                    augmented[k][j] -= factor * augmented[i][j];
                }
            }
        }
    }

    vector<vector<double>> inverse(sz, vector<double>(sz, 0));
    for (size_t i = 0; i < sz; ++i)
    {
        for (size_t j = 0; j < sz; ++j)
        {
            inverse[i][j] = augmented[i][j + sz];
        }
    }

    Matrix inversed_matrix(inverse);
    return inversed_matrix;
}

double Matrix::determinant(const vector<vector<double>> &mat) const
{
    size_t sz = mat.size();
    if (sz == 0 || mat[0].size() != sz)
        throw std::invalid_argument("determinant: matrix must be square and non-empty");

    if (sz == 1) return mat[0][0];
    if (sz == 2) return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];

    double det = 0.0;
    for (size_t col = 0; col < sz; ++col)
    {
        vector<vector<double>> sub(sz - 1, vector<double>(sz - 1));
        for (size_t i = 1; i < sz; ++i)
        {
            size_t subcol = 0;
            for (size_t j = 0; j < sz; ++j)
            {
                if (j == col) continue;
                sub[i - 1][subcol++] = mat[i][j];
            }
        }
        double minor = determinant(sub);
        det += ((col % 2 == 0) ? 1.0 : -1.0) * mat[0][col] * minor;
    }
    return det;
}

bool Matrix::operator==(const Matrix &rhs) const
{
    if (line_indexes.size() != rhs.line_indexes.size() ||
        column_indexes.size() != rhs.column_indexes.size())
        return false;

    for (size_t i = 0; i < line_indexes.size(); ++i)
        for (size_t j = 0; j < column_indexes.size(); ++j)
            if (matrix[line_indexes[i]][column_indexes[j]] !=
                rhs.matrix[rhs.line_indexes[i]][rhs.column_indexes[j]])
                return false;

    return true;
}

bool Matrix::operator!=(const Matrix &rhs) const { return !(*this == rhs); }

bool Matrix::operator<(const Matrix &rhs) const
{
    if (rhs.line_indexes.size() != line_indexes.size() || column_indexes.size() != rhs.column_indexes.size())
        return false;

    for (size_t i = 0; i < line_indexes.size(); ++i)
        for (size_t j = 0; j < column_indexes.size(); ++j)
            if (!(matrix[line_indexes[i]][column_indexes[j]] < rhs.matrix[rhs.line_indexes[i]][rhs.column_indexes[j]]))
                return false;

    return true;
}

Matrix Matrix::subtract(const Matrix &a) const
{
    size_t rows = line_indexes.size();
    size_t cols = column_indexes.size();
    if (rows != a.line_indexes.size() || cols != a.column_indexes.size())
        throw std::invalid_argument("subtract: incompatible sizes");

    vector<vector<double>> ans(rows, vector<double>(cols, 0.0));
    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            double left = matrix[line_indexes[i]][column_indexes[j]];
            double right = a.matrix[a.line_indexes[i]][a.column_indexes[j]];
            ans[i][j] = std::round((left - right) * 10000000.0) / 10000000.0;
        }
    }
    return Matrix(ans);
}

bool Matrix::operator>(double num) const
{
    for (size_t i = 0; i < line_indexes.size(); ++i)
        for (size_t j = 0; j < column_indexes.size(); ++j)
            if (!(matrix[line_indexes[i]][column_indexes[j]] > num))
                return false;
    return true;
}

bool Matrix::operator<(double num) const
{
    for (size_t i = 0; i < line_indexes.size(); ++i)
        for (size_t j = 0; j < column_indexes.size(); ++j)
            if (!(matrix[line_indexes[i]][column_indexes[j]] < num))
                return false;
    return true;
}

Matrix Matrix::operator*(double num) const
{
    vector<vector<double>> out = matrix;
    for (size_t i = 0; i < out.size(); ++i)
        for (size_t j = 0; j < out[0].size(); ++j)
            out[i][j] *= num;
    return Matrix(out);
}

bool Matrix::operator<=(double num) const
{
    for (size_t i = 0; i < line_indexes.size(); ++i)
        for (size_t j = 0; j < column_indexes.size(); ++j)
            if (!(matrix[line_indexes[i]][column_indexes[j]] <= num))
                return false;
    return true;
}

bool Matrix::operator>=(double num) const
{
    for (size_t i = 0; i < line_indexes.size(); ++i)
        for (size_t j = 0; j < column_indexes.size(); ++j)
            if (!(matrix[line_indexes[i]][column_indexes[j]] >= num))
                return false;
    return true;
}

std::vector<int> Matrix::get_addition_to_square_matrix(const std::vector<int> &available_indexes)
{
    int need_to_add = (int)line_indexes.size() - (int)column_indexes.size();
    if (need_to_add <= 0)
    {
        return column_indexes;
    }

    auto combs = combinations(available_indexes, need_to_add);
    for (auto &comb : combs)
    {
        auto candidate_cols = concatenate_vectors<int>(column_indexes, comb);
        Matrix sub = allocate_matrix(line_indexes, candidate_cols);
        double det = determinant(sub.matrix);
        if (std::fabs(det) > 1e-12)
        {
            set_columns(candidate_cols);
            return column_indexes;
        }
    }
    return std::vector<int>{};
}

int Matrix::column_size() const
{
    if (matrix.empty()) return 0;
    return (int)matrix[0].size();
}

void Matrix::gaussian_elimination()
{
    size_t rows = line_indexes.size();
    size_t cols = column_indexes.size();
    size_t r = 0;

    for (size_t c = 0; c < cols && r < rows; ++c)
    {
        size_t pivot_row = r;
        double max_val = std::fabs(matrix[line_indexes[pivot_row]][column_indexes[c]]);
        for (size_t i = r + 1; i < rows; ++i)
        {
            double val = std::fabs(matrix[line_indexes[i]][column_indexes[c]]);
            if (val > max_val)
            {
                max_val = val;
                pivot_row = i;
            }
        }

        if (std::fabs(matrix[line_indexes[pivot_row]][column_indexes[c]]) < 1e-14)
        {
            continue;
        }

        if (pivot_row != r)
            std::swap(line_indexes[pivot_row], line_indexes[r]);

        double pivot = matrix[line_indexes[r]][column_indexes[c]];
        for (size_t j = c; j < cols; ++j)
        {
            matrix[line_indexes[r]][column_indexes[j]] /= pivot;
        }

        for (size_t i = r + 1; i < rows; ++i)
        {
            double factor = matrix[line_indexes[i]][column_indexes[c]];
            for (size_t j = c; j < cols; ++j)
            {
                matrix[line_indexes[i]][column_indexes[j]] -= factor * matrix[line_indexes[r]][column_indexes[j]];
            }
        }

        ++r;
    }
}

int Matrix::compute_rank() const
{
    Matrix temp = *this;
    temp.gaussian_elimination();

    int rank = 0;
    for (size_t i = 0; i < temp.line_indexes.size(); ++i)
    {
        bool nonzero = false;
        for (size_t j = 0; j < temp.column_indexes.size(); ++j)
        {
            if (std::fabs(temp.matrix[temp.line_indexes[i]][temp.column_indexes[j]]) > 1e-12)
            {
                nonzero = true;
                break;
            }
        }
        if (nonzero) ++rank;
    }
    return rank;
}

bool Matrix::is_full_rank() const
{
    int rank = compute_rank();
    int min_dim = (int)std::min(line_indexes.size(), column_indexes.size());
    return rank == min_dim;
}
