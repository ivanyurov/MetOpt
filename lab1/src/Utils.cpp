#include "../include/Utils.h"
#include <iostream>

void generate_combinations(const std::vector<int> &elements, int n, size_t start,
                           std::vector<int> &current, std::vector<std::vector<int>> &result)
{
    if ((int)current.size() == n)
    {
        result.push_back(current);
        return;
    }

    for (size_t i = start; i < elements.size(); ++i)
    {
        current.push_back(elements[i]);
        generate_combinations(elements, n, i + 1, current, result);
        current.pop_back();
    }
}

std::vector<std::vector<int>> combinations(const std::vector<int> &elements, int n)
{
    std::vector<std::vector<int>> result;
    if (n <= 0 || (size_t)n > elements.size()) return result;
    std::vector<int> current;
    generate_combinations(elements, n, 0, current, result);
    return result;
}

void print_vector(const std::vector<int> &a)
{
    for (size_t i = 0; i < a.size(); ++i)
    {
        std::cout << a[i];
        if (i + 1 < a.size()) std::cout << " ";
    }
    std::cout << "\n";
}
