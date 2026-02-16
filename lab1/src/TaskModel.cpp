#include <iostream>
#include "../include/TaskModel.h"
#include <stdexcept>
#include <algorithm>
#include <iomanip>

using std::cout;
using std::endl;

LPProblem::LPProblem()
{
    n = 0;
}

LPProblem::LPProblem(int n)
{
    this->n = n;
}

LPProblem::LPProblem(const LPProblem &problem)
{
    this->n = problem.n;
    this->to_max = problem.to_max;
    this->objective_type = problem.objective_type;
    this->objective = problem.objective;
    this->constraints = problem.constraints;
    this->bounds = problem.bounds;
}

LPProblem::LPProblem(vector<double> objective, ObjectiveType objective_type)
{
    this->n = static_cast<int>(objective.size());
    this->objective = objective;
    this->objective_type = objective_type;
}

void LPProblem::set_solution_dim(int n)
{
    this->n = n;
}

void LPProblem::set_objective(const std::vector<double> &coeffs, ObjectiveType type)
{
    objective = coeffs;
    objective_type = type;
    n = static_cast<int>(objective.size());
}

void LPProblem::add_constraint(const Constraint &c)
{
    if (!objective.empty() && static_cast<int>(c.coefficients.size()) != n)
        throw std::runtime_error("add_constraint: coefficients size != n");
    constraints.push_back(c);
}

void LPProblem::add_var_bound(const VariableBound &vb)
{
    if (vb.component < 1)
        throw std::runtime_error("add_var_bound: component must be >= 1");
    // We allow adding bounds even if n not set yet; validation will occur on convert or later
    bounds.push_back(vb);
}

void LPProblem::print_problem() const
{
    cout << "Problem of " << ((objective_type == ObjectiveType::MINIMIZE) ? "Minimization" : "Maximization") << "\n";
    cout << "Objective:\n";

    // Печатаем коэффициент с явным знаком, но индекс без знака:
    for (int i = 0; i < static_cast<int>(objective.size()); ++i)
    {
        // Печатать знак только для коэффициента:
        cout << std::showpos << objective[i] << std::noshowpos;
        cout << "x" << (i + 1) << " ";
    }
    cout << "\n\n";

    cout << "Constraints:\n";
    for (const auto &c : constraints)
    {
        for (int i = 0; i < static_cast<int>(c.coefficients.size()); ++i)
        {
            cout << std::showpos << c.coefficients[i] << std::noshowpos;
            cout << "x" << (i + 1) << " ";
        }
        if (c.type == InequalityType::LESS_EQUAL)
            cout << "<= ";
        else if (c.type == InequalityType::GREATER_EQUAL)
            cout << ">= ";
        else
            cout << "= ";

        cout << c.b << "\n";
    }
    cout << "\nBounds:\n";
    cout << "\nBounds:\n";
    for (const auto &b : bounds)
    {
        cout << "x" << b.component;
        if (b.type == BoundType::NOT_NEGATIVE)
            cout << " >= 0";
        else if (b.type == BoundType::NOT_POSITIVE)
            cout << " <= 0";
        else // BoundType::NO
            cout << " free";
        cout << "\n";
    }
    cout << std::noshowpos << "\n";
}

std::vector<double> LPProblem::get_objective() const
{
    return objective;
}

std::vector<Constraint> LPProblem::get_constraints() const
{
    return constraints;
}

//
// LPProblemGeneral
//

void LPProblemGeneral::convert()
{
    // 1) MAX -> MIN
    if (objective_type == ObjectiveType::MAXIMIZE)
    {
        for (auto &v : objective)
            v *= -1;
        objective_type = ObjectiveType::MINIMIZE;
        to_max = -1;
    }

    // 2) <= -> >=
    for (auto &c : constraints)
    {
        if (c.type == InequalityType::LESS_EQUAL)
        {
            for (auto &coef : c.coefficients)
                coef *= -1;
            c.b *= -1;
            c.type = InequalityType::GREATER_EQUAL;
        }
    }

    // 3) x <= 0 -> -x >= 0
    for (auto &b : bounds)
    {
        if (b.type == BoundType::NOT_POSITIVE)
        {
            int idx = b.component - 1;
            if (idx < 0 || idx >= static_cast<int>(objective.size()))
                continue;

            negated_vars.push_back(idx); // ← ЗАПОМИНАЕМ индекс

            b.type = BoundType::NOT_NEGATIVE;
            for (auto &c : constraints)
            {
                if (idx < static_cast<int>(c.coefficients.size()))
                    c.coefficients[idx] *= -1;
            }
            objective[idx] *= -1;
        }
    }

    // 4) add explicit NO bounds if missing
    for (int i = 1; i <= n; ++i)
    {
        bool found = false;
        for (const auto &b : bounds)
            if (b.component == i)
            {
                found = true;
                break;
            }
        if (!found)
            bounds.push_back({i, BoundType::NO});
    }
}

// Также обнови get_initial_solution в LPProblemGeneral (чтобы работало восстановление после standard form)
std::vector<double> LPProblemGeneral::get_initial_solution(std::vector<double> solution)
{
    std::vector<double> ans(initial_dim, 0.0);
    for (int i = 0; i < initial_dim; ++i)
        ans[i] = solution[i];

    // 1. Восстанавливаем свободные переменные (x = x⁺ - x⁻)
    for (const auto &p : bounds_to_subtract)
    {
        if (p.second < static_cast<int>(solution.size()))
            ans[p.first] -= solution[p.second];
    }

    // 2. Восстанавливаем переменные, бывшие NOT_POSITIVE (x = -x')
    for (int idx : negated_vars)
    {
        if (idx < static_cast<int>(ans.size()))
            ans[idx] *= -1.0; // ← КРИТИЧЕСКОЕ ИСПРАВЛЕНИЕ!
    }

    // 3. Учёт исходного знака (если была максимизация)
    for (auto &v : ans)
        v *= to_max;

    return ans;
}

std::unique_ptr<LPProblem> LPProblemGeneral::dual()
{
    // Работаем с КОПИЕЙ без вызова convert() — сохраняем исходные типы ограничений!
    LPProblemGeneral copy = *this;

    int m = static_cast<int>(copy.constraints.size());
    int nvars = copy.n;

    auto dual = std::make_unique<LPProblemGeneral>(m);

    // objective = b (правые части прямой задачи)
    std::vector<double> dual_obj(m);
    for (int i = 0; i < m; ++i)
        dual_obj[i] = copy.constraints[i].b;
    dual->set_objective(dual_obj, ObjectiveType::MAXIMIZE);

    // bounds for dual variables — КРИТИЧЕСКИ ВАЖНО: учитываем ИСХОДНЫЙ тип ограничения!
    for (int i = 0; i < m; ++i)
    {
        VariableBound vb;
        vb.component = i + 1;

        // Для задачи МИНИМИЗАЦИИ:
        //   >=  → двойственная переменная >= 0
        //   <=  → двойственная переменная <= 0   ← ЭТО БЫЛО УТЕРЯНО!
        //   =   → двойственная переменная свободна
        if (copy.constraints[i].type == InequalityType::GREATER_EQUAL)
            vb.type = BoundType::NOT_NEGATIVE; // y_i >= 0
        else if (copy.constraints[i].type == InequalityType::LESS_EQUAL)
            vb.type = BoundType::NOT_POSITIVE; // y_i <= 0  ← ИСПРАВЛЕНО!
        else                                   // EQUAL
            vb.type = BoundType::NO;           // y_i free

        dual->add_var_bound(vb);
    }

    // constraints of dual — для каждой переменной прямой задачи
    for (int j = 0; j < nvars; ++j)
    {
        Constraint dc;
        dc.coefficients.resize(m);
        for (int i = 0; i < m; ++i)
        {
            if (j < static_cast<int>(copy.constraints[i].coefficients.size()))
                dc.coefficients[i] = copy.constraints[i].coefficients[j];
            else
                dc.coefficients[i] = 0.0;
        }
        if (j < static_cast<int>(copy.objective.size()))
            dc.b = copy.objective[j];
        else
            dc.b = 0.0;

        // Тип ограничения двойственной задачи зависит от границ переменной прямой задачи
        dc.type = InequalityType::EQUAL;
        for (const auto &pb : copy.bounds)
        {
            if (pb.component - 1 == j)
            {
                if (pb.type == BoundType::NOT_NEGATIVE)      // x_j >= 0
                    dc.type = InequalityType::LESS_EQUAL;    // для максимизации: a^T y <= c_j
                else if (pb.type == BoundType::NOT_POSITIVE) // x_j <= 0
                    dc.type = InequalityType::GREATER_EQUAL; // для максимизации: a^T y >= c_j
                else                                         // x_j свободна
                    dc.type = InequalityType::EQUAL;         // a^T y = c_j
                break;
            }
        }

        dual->add_constraint(dc);
    }

    return dual;
}

void LPProblemGeneral::to_standard_form()
{
    // Сохраняем исходную размерность для восстановления решения
    if (initial_dim == 0)
        initial_dim = n;

    // 1. Max → Min
    if (objective_type == ObjectiveType::MAXIMIZE)
    {
        for (auto &coef : objective)
            coef *= -1.0;
        objective_type = ObjectiveType::MINIMIZE;
        to_max = -1;
    }

    // 2. <= → >= (переворачиваем знаки)
    for (auto &c : constraints)
    {
        if (c.type == InequalityType::LESS_EQUAL)
        {
            for (auto &coef : c.coefficients)
                coef *= -1.0;
            c.b *= -1.0;
            c.type = InequalityType::GREATER_EQUAL;
        }
    }

    // 3. x <= 0 → x >= 0 (переворачиваем коэффициенты переменной)
    for (auto &b : bounds)
    {
        if (b.type == BoundType::NOT_POSITIVE)
        {
            int idx = b.component - 1;
            if (idx >= 0 && idx < n)
            {
                b.type = BoundType::NOT_NEGATIVE;
                for (auto &c : constraints)
                {
                    if (idx < static_cast<int>(c.coefficients.size()))
                        c.coefficients[idx] *= -1.0;
                }
                objective[idx] *= -1.0;
            }
        }
    }

    // 4. Добавляем недостающие bounds (все переменные должны иметь bound)
    for (int i = 1; i <= n; ++i)
    {
        bool has_bound = false;
        for (const auto &b : bounds)
        {
            if (b.component == i)
            {
                has_bound = true;
                break;
            }
        }
        if (!has_bound)
            bounds.push_back({i, BoundType::NO});
    }

    // 5. Собираем свободные переменные заранее
    std::vector<int> free_var_indices;
    for (const auto &b : bounds)
    {
        if (b.type == BoundType::NO)
            free_var_indices.push_back(b.component - 1);
    }

    // 6. Расщепляем свободные переменные x = x+ - x-
    for (int orig_idx : free_var_indices)
    {
        if (orig_idx < 0 || orig_idx >= n)
            continue;

        // Оригинальная становится x+ >= 0
        for (auto &b : bounds)
        {
            if (b.component - 1 == orig_idx)
            {
                b.type = BoundType::NOT_NEGATIVE;
                break;
            }
        }

        // Добавляем x- >= 0
        bounds.push_back({n + 1, BoundType::NOT_NEGATIVE});
        bounds_to_subtract.push_back({orig_idx, n});

        // Столбец x- = -столбец оригинальной
        for (auto &c : constraints)
        {
            double orig_coef = (orig_idx < static_cast<int>(c.coefficients.size())) ? c.coefficients[orig_idx] : 0.0;
            c.coefficients.push_back(-orig_coef);
        }

        // Коэффициент в цели для x- = -коэффициент оригинальной
        double orig_obj = (orig_idx < static_cast<int>(objective.size())) ? objective[orig_idx] : 0.0;
        objective.push_back(-orig_obj);

        ++n;
    }

    // 7. Расщепляем равенства = на два >=
    std::vector<Constraint> new_constraints;
    new_constraints.reserve(constraints.size() * 2);

    for (const auto &c : constraints)
    {
        if (c.type == InequalityType::EQUAL)
        {
            // >= b
            Constraint positive = c;
            positive.type = InequalityType::GREATER_EQUAL;
            new_constraints.push_back(positive);

            // >= -b (перевёрнутое)
            Constraint negative = c;
            for (auto &coef : negative.coefficients)
                coef *= -1.0;
            negative.b *= -1.0;
            negative.type = InequalityType::GREATER_EQUAL;
            new_constraints.push_back(negative);
        }
        else
        {
            // Уже >= — оставляем как есть (b может быть < 0)
            new_constraints.push_back(c);
        }
    }

    constraints = std::move(new_constraints);
}

void LPProblemSlack::convert()
{
    // LPProblemGeneral::convert();

    if (initial_dim == 0)
        initial_dim = n;

    // 2. Расщепляем свободные переменные: x = x⁺ - x⁻
    std::vector<int> free_vars;
    for (const auto &b : bounds)
        if (b.type == BoundType::NO)
            free_vars.push_back(b.component - 1);

    for (int idx : free_vars)
    {
        if (idx < 0 || idx >= n)
            continue;

        // Оригинальная переменная становится x⁺ >= 0
        for (auto &b : bounds)
            if (b.component - 1 == idx)
            {
                b.type = BoundType::NOT_NEGATIVE;
                break;
            }

        // Добавляем x⁻ >= 0
        bounds_to_subtract.push_back({idx, n});

        for (auto &c : constraints)
        {
            double coef = (idx < static_cast<int>(c.coefficients.size()))
                              ? c.coefficients[idx]
                              : 0.0;
            c.coefficients.push_back(-coef);
        }

        objective.push_back(-objective[idx]);
        bounds.push_back({n + 1, BoundType::NOT_NEGATIVE});
        ++n;
    }

    // 3. Добавляем slack-переменные для преобразования неравенств в равенства
    // Для >=: aᵀx >= b  →  aᵀx - s = b  (s >= 0)
    // Для <=: aᵀx <= b  →  aᵀx + s = b  (s >= 0)
    int m = static_cast<int>(constraints.size());
    for (int i = 0; i < m; ++i)
    {
        if (constraints[i].type == InequalityType::GREATER_EQUAL ||
            constraints[i].type == InequalityType::LESS_EQUAL)
        {
            // Добавляем столбец для новой slack-переменной
            for (int j = 0; j < m; ++j)
            {
                if (j == i)
                {
                    // Для >=: коэффициент -1; для <=: коэффициент +1
                    double slack_coef = (constraints[i].type == InequalityType::GREATER_EQUAL)
                                            ? -1.0
                                            : 1.0;
                    constraints[j].coefficients.push_back(slack_coef);
                }
                else
                {
                    constraints[j].coefficients.push_back(0.0);
                }
            }

            objective.push_back(0.0);
            bounds.push_back({n + 1, BoundType::NOT_NEGATIVE});
            ++n;

            // Преобразуем в равенство
            constraints[i].type = InequalityType::EQUAL;
        }
    }

    // 4. Нормализуем правые части: обеспечиваем b >= 0
    // Если b < 0, умножаем всё уравнение на -1
    for (auto &c : constraints)
    {
        if (c.b < 0)
        {
            for (auto &coef : c.coefficients)
                coef *= -1.0;
            c.b *= -1.0;
            // Тип остаётся EQUAL, знаки коэффициентов изменены
        }
    }
}

std::unique_ptr<LPProblem> LPProblemSlack::dual()
{
    auto base_dual = LPProblemGeneral::dual();
    auto dual_slack = std::make_unique<LPProblemSlack>(base_dual->get_objective(), base_dual->objective_type);

    for (const auto &c : base_dual->get_constraints())
        dual_slack->add_constraint(c);
    for (const auto &b : base_dual->bounds)
        dual_slack->add_var_bound(b);

    return dual_slack;
}

std::vector<double> LPProblemSlack::get_initial_solution(std::vector<double> solution)
{
    if (static_cast<int>(solution.size()) < n)
        throw std::runtime_error("get_initial_solution: solution size too small");

    std::vector<double> ans(initial_dim, 0.0);
    for (int i = 0; i < initial_dim; ++i)
        ans[i] = solution[i];

    for (const auto &p : bounds_to_subtract)
    {
        if (p.first >= 0 && p.first < static_cast<int>(ans.size()) && p.second >= 0 && p.second < static_cast<int>(solution.size()))
            ans[p.first] -= solution[p.second];
        else
            throw std::runtime_error("get_initial_solution: index out of range");
    }

    return ans;
}
