#pragma once
#include <vector>
#include <algorithm>
#include <iostream>

using std::vector;

template <typename T>
std::vector<T> subtract_vectors(const std::vector<T> &a, const std::vector<T> &b)
{
    std::vector<T> res;
    res.reserve(a.size());
    for (size_t i = 0; i < a.size(); ++i)
    {
        if (std::find(b.begin(), b.end(), a[i]) == b.end())
        {
            res.push_back(a[i]);
        }
    }
    return res;
}

template <typename T>
std::vector<T> concatenate_vectors(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
    std::vector<T> result;
    result.reserve(vec1.size() + vec2.size());
    for (const auto &el : vec1)
        result.push_back(el);
    for (const auto &el : vec2)
        result.push_back(el);
    // Не сортируем по умолчанию: порядок индексов важен для соответствия колонок.
    return result;
}

// Генерация всех сочетаний (комбинаций) size n из elements
void generate_combinations(const std::vector<int> &elements, int n, size_t start,
                           std::vector<int> &current, std::vector<std::vector<int>> &result);

// Удобная оболочка
std::vector<std::vector<int>> combinations(const std::vector<int> &elements, int n);

// Печать вектора int
void print_vector(const std::vector<int> &a);
