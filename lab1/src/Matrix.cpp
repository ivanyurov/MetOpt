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
    this->matrix = matrix;
    n = matrix.size();
    m = matrix[0].size();
    for (int i = 0; i < n; i++)
        line_indexes.push_back(i);
    for (int j = 0; j < m; j++)
        column_indexes.push_back(j);
}

Matrix::Matrix(const vector<double> &vect)
{
    for (auto &el : vect)
    {
        matrix.push_back(vector<double>(1, el));
    }

    n = matrix.size();
    m = matrix[0].size();
    for (int i = 0; i < n; i++)
    {
        line_indexes.push_back(i);
    }
    for (int j = 0; j < m; j++)
    {
        column_indexes.push_back(j);
    }
}

// Операция транспонирования
Matrix &Matrix::transpose()
{
    vector<vector<double>> t_matrix(m, vector<double>(n));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            t_matrix[j][i] = matrix[i][j];
        }
    }
    Matrix *ret = new Matrix(t_matrix);
    ret->set_columns(line_indexes);
    ret->set_lines(column_indexes);
    return *ret;
}

// Вывод матрицы
void Matrix::print()
{
    cout << endl;
    for (auto &line_ind : line_indexes)
    {
        for (auto &col_ind : column_indexes)
        {
            std::cout << matrix[line_ind][col_ind] << " ";
        }
        std::cout << std::endl;
    }
    cout << "\n";
}

// Установка индексов столбцов и строк
void Matrix::set_columns(vector<int> columns)
{
    column_indexes = columns;
}

void Matrix::set_lines(vector<int> lines)
{
    line_indexes = lines;
}

// Выделение подматрицы
Matrix &Matrix::allocate_matrix(vector<int> lines, vector<int> columns)
{
    vector<vector<double>> new_matrix(lines.size(), vector<double>(columns.size()));
    for (size_t i = 0; i < lines.size(); i++)
    {
        for (size_t j = 0; j < columns.size(); j++)
        {
            new_matrix[i][j] = matrix[lines[i]][columns[j]];
        }
    }
    Matrix *ret = new Matrix(new_matrix);
    return *ret;
}

// Умножение матриц
Matrix &Matrix::multiply(const Matrix &a)
{
    vector<vector<double>> result(line_indexes.size(), vector<double>(a.column_indexes.size(), 0));

    for (int i = 0; i < line_indexes.size(); i++)
    {
        vector<double> line = matrix[line_indexes[i]];
        for (int k = 0; k < a.column_indexes.size(); k++)
        {

            for (int j = 0; j < column_indexes.size(); j++)
            {
                result[i][k] += (line[column_indexes[j]] * a.matrix[a.line_indexes[j]][a.column_indexes[k]]);
            }
        }
    }
    Matrix *ans = new Matrix(result);
    return *ans;
}

Matrix &Matrix::get_inverse_matrix()
{
    if (line_indexes.size() == 0 || line_indexes.size() != column_indexes.size())
    {
        throw std::invalid_argument("Matrix must be square");
    }
    vector<vector<double>> augmented(line_indexes.size(), vector<double>(2 * line_indexes.size(), 0));

    for (size_t i = 0; i < line_indexes.size(); ++i)
    {
        for (size_t j = 0; j < column_indexes.size(); ++j)
        {
            augmented[i][j] = matrix[line_indexes[i]][column_indexes[j]];
        }
        augmented[i][i + line_indexes.size()] = 1;
    }

    for (size_t i = 0; i < line_indexes.size(); ++i)
    {
        size_t max_row = i;
        for (size_t k = i + 1; k < column_indexes.size(); ++k)
        {
            if (abs(augmented[k][i]) > abs(augmented[max_row][i]))
            {
                max_row = k;
            }
        }

        if (max_row != i)
        {
            std::swap(augmented[i], augmented[max_row]);
        }

        if (augmented[i][i] == 0)
        {
            throw std::invalid_argument("Matrix is singular and cannot be inverted");
        }

        double pivot = augmented[i][i];
        for (size_t j = 0; j < 2 * line_indexes.size(); ++j)
        {
            augmented[i][j] = (augmented[i][j] / pivot);
        }

        for (size_t k = 0; k < line_indexes.size(); ++k)
        {
            if (k != i)
            {
                double factor = augmented[k][i];
                for (size_t j = 0; j < 2 * line_indexes.size(); ++j)
                {
                    augmented[k][j] -= factor * augmented[i][j];
                }
            }
        }
    }

    vector<vector<double>> inverse(line_indexes.size(), vector<double>(line_indexes.size(), 0));
    for (size_t i = 0; i < line_indexes.size(); ++i)
    {
        for (size_t j = 0; j < column_indexes.size(); ++j)
        {
            inverse[i][j] = augmented[i][j + line_indexes.size()];
        }
    }

    Matrix *inversed_matrix = new Matrix(inverse);

    return *inversed_matrix;
}

// Детерминант (рекурсивная реализация)
double Matrix::determinant(vector<vector<double>> matrix)
{
    if (matrix.size() == 0 || matrix[0].size() != matrix.size())
    {
        throw std::invalid_argument("Matrix must be square");
    }

    if (matrix.size() == 1)
    {
        return matrix[0][0];
    }

    if (matrix.size() == 2)
    {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }

    if (matrix.size() == 3)
    {
        return matrix[0][0] * matrix[1][1] * matrix[2][2] +
               matrix[0][1] * matrix[1][2] * matrix[2][0] +
               matrix[0][2] * matrix[1][0] * matrix[2][1] -
               matrix[0][2] * matrix[1][1] * matrix[2][0] -
               matrix[0][0] * matrix[1][2] * matrix[2][1] -
               matrix[0][1] * matrix[1][0] * matrix[2][2];
    }

    double det = 0;
    for (size_t col = 0; col < matrix.size(); ++col)
    {
        vector<vector<double>> submatrix(matrix.size() - 1, vector<double>(matrix.size() - 1));
        for (size_t i = 1; i < matrix.size(); ++i)
        {
            size_t subcol = 0;
            for (size_t j = 0; j < matrix.size(); ++j)
            {
                if (j == col)
                    continue;
                submatrix[i - 1][subcol] = matrix[i][j];
                subcol++;
            }
        }

        double minor = determinant(submatrix);
        det += (col % 2 == 0 ? 1 : -1) * matrix[0][col] * minor;
    }

    return det;
}

// Операции сравнения
bool Matrix::operator==(const Matrix &rhs)
{
    if (line_indexes.size() != rhs.line_indexes.size() ||
        column_indexes.size() != rhs.column_indexes.size())
        return false;

    for (size_t i = 0; i < line_indexes.size(); ++i)
    {
        for (size_t j = 0; j < column_indexes.size(); ++j)
        {
            if (matrix[line_indexes[i]][column_indexes[j]] !=
                rhs.matrix[rhs.line_indexes[i]][rhs.column_indexes[j]])
            {
                return false;
            }
        }
    }
    return true;
}

bool Matrix::operator!=(const Matrix &rhs)
{
    return !(*this == rhs);
}

// Остальные операторы и методы
Matrix &Matrix::subtract(Matrix &a)
{
    vector<vector<double>> ans(line_indexes.size(), vector<double>(column_indexes.size()));
    for (int i = 0; i < line_indexes.size(); i++)
    {
        for (int j = 0; j < column_indexes.size(); j++)
        {
            ans[i][j] = std::round((matrix[line_indexes[i]][column_indexes[j]] - a.matrix[a.line_indexes[i]][a.column_indexes[j]]) * 10000000) / 10000000;
        }
    }
    Matrix *ret = new Matrix(ans);
    return *ret;
}

// Операторы сравнения с числами
bool Matrix::operator>(double num)
{
    for (int i = 0; i < line_indexes.size(); i++)
    {
        for (int j = 0; j < column_indexes.size(); j++)
        {
            if (matrix[line_indexes[i]][column_indexes[j]] <= num)
            {
                return false;
            }
        }
    }
    return true;
}

bool Matrix::operator<(double num)
{
    for (int i = 0; i < line_indexes.size(); i++)
    {
        for (int j = 0; j < column_indexes.size(); j++)
        {
            if (matrix[line_indexes[i]][column_indexes[j]] >= num)
            {
                return false;
            }
        }
    }
    return true;
}

bool Matrix::operator<(const Matrix &rhs)
{
    if (rhs.line_indexes.size() != line_indexes.size() || column_indexes.size() != rhs.column_indexes.size())
    {
        return false;
    }
    for (int i = 0; i < line_indexes.size(); i++)
    {
        for (int j = 0; j < column_indexes.size(); j++)
        {
            if (matrix[line_indexes[i]][column_indexes[j]] >= rhs.matrix[rhs.line_indexes[i]][rhs.column_indexes[j]])
            {
                return false;
            }
        }
    }
    return true;
}

// Умножение на скаляр
Matrix &Matrix::operator*(double num)
{
    for (int i = 0; i < line_indexes.size(); i++)
    {
        for (int j = 0; j < column_indexes.size(); j++)
        {
            matrix[line_indexes[i]][column_indexes[j]] *= num;
        }
    }
    return *this;
}

bool Matrix::operator<=(double num)
{
    for (int i = 0; i < line_indexes.size(); i++)
    {
        for (int j = 0; j < column_indexes.size(); j++)
        {
            if (matrix[line_indexes[i]][column_indexes[j]] > num)
            {
                return false;
            }
        }
    }
    return true;
}
bool Matrix::operator>=(double num)
{
    for (int i = 0; i < line_indexes.size(); i++)
    {
        for (int j = 0; j < column_indexes.size(); j++)
        {
            if (matrix[line_indexes[i]][column_indexes[j]] < num)
            {
                return false;
            }
        }
    }
    return true;
}

// Дополнительные методы
vector<int> Matrix::get_addition_to_square_matrix(vector<int> &available_indexes)
{

    int need_to_add = line_indexes.size() - column_indexes.size();
    if (need_to_add <= 0)
    {
        return column_indexes;
    }
    vector<vector<int>> combs = combinations(available_indexes, need_to_add);
    for (auto &comb : combs)
    {

        if (determinant(
                allocate_matrix(
                    line_indexes,
                    concatenate_vectors<int>(column_indexes, comb))
                    .matrix) != 0)
        {

            set_columns(concatenate_vectors<int>(column_indexes, comb));

            return column_indexes;
        }
    }
}

int Matrix::column_size()
{
    if (matrix.empty())
        return 0;
    return matrix[0].size();
}

void Matrix::gaussian_elimination()
{
    int rank = 0;
    for (size_t col = 0; col < column_indexes.size() && rank < line_indexes.size(); ++col)
    {
        // Поиск строки с максимальным элементом в текущем столбце
        size_t max_row = rank;
        for (size_t row = rank + 1; row < line_indexes.size(); ++row)
        {
            if (abs(matrix[line_indexes[row]][column_indexes[col]]) >
                abs(matrix[line_indexes[max_row]][column_indexes[col]]))
            {
                max_row = row;
            }
        }

        // Если максимальный элемент равен нулю, пропускаем столбец
        if (matrix[line_indexes[max_row]][column_indexes[col]] == 0)
        {
            continue;
        }

        // Меняем строки местами
        if (max_row != rank)
        {
            std::swap(line_indexes[max_row], line_indexes[rank]);
        }

        // Нормализация текущей строки
        double pivot = matrix[line_indexes[rank]][column_indexes[col]];
        for (size_t j = col; j < column_indexes.size(); ++j)
        {
            matrix[line_indexes[rank]][column_indexes[j]] /= pivot;
        }

        // Обнуление элементов ниже текущего
        for (size_t row = rank + 1; row < line_indexes.size(); ++row)
        {
            double factor = matrix[line_indexes[row]][column_indexes[col]];
            for (size_t j = col; j < column_indexes.size(); ++j)
            {
                matrix[line_indexes[row]][column_indexes[j]] -=
                    factor * matrix[line_indexes[rank]][column_indexes[j]];
            }
        }

        // Увеличиваем ранг
        rank++;
    }
}

int Matrix::compute_rank()
{
    // Создаем копию матрицы, чтобы не изменять оригинальную
    Matrix temp(*this);
    temp.gaussian_elimination();

    // Подсчет ненулевых строк
    int rank = 0;
    for (size_t row = 0; row < temp.line_indexes.size(); ++row)
    {
        bool is_non_zero = false;
        for (size_t col = 0; col < temp.column_indexes.size(); ++col)
        {
            if (temp.matrix[temp.line_indexes[row]][temp.column_indexes[col]] != 0)
            {
                is_non_zero = true;
                break;
            }
        }
        if (is_non_zero)
        {
            rank++;
        }
    }
    return rank;
}

bool Matrix::is_full_rank()
{
    int rank = compute_rank();
    int min_dim = std::min(line_indexes.size(), column_indexes.size());
    return rank == min_dim;
}