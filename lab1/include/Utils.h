#pragma once
#include <vector>
#include <algorithm>
#include <iostream>

using std::vector;

template <typename T>
std::vector<T> subtract_vectors(std::vector<T> a, std::vector<T> b)
{
    std::vector<T> res;
    for (int i = 0; i < a.size(); i++)
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
    vector<T> result;
    for (auto &el : vec1)
    {
        result.push_back(el);
    }
    for (auto &el : vec2)
    {
        result.push_back(el);
    }
    std::sort(result.begin(), result.end());
    return result;
}

void generate_combinations(const std::vector<int> &elements, int n, size_t start,
                           std::vector<int> &current, std::vector<std::vector<int>> &result);

std::vector<std::vector<int>> combinations(const std::vector<int> &elements, int n);

void print_vector(const std::vector<int> &a);