#pragma once

#include <vector>
#include <memory>

using std::vector;

enum class ObjectiveType
{
    MINIMIZE,
    MAXIMIZE
};

enum class InequalityType
{
    LESS_EQUAL,
    GREATER_EQUAL,
    EQUAL
};

enum class BoundType
{
    NOT_NEGATIVE,
    NOT_POSITIVE,
    NO
};

// Ограничение: a1*x1 + ... + an*xn <=/>=/= b
struct Constraint
{
    vector<double> coefficients;
    double b = 0.0;
    InequalityType type = InequalityType::EQUAL;
};

// Граничное условие на переменную (component — 1-based индекс)
struct VariableBound
{
    int component = 0;
    BoundType type = BoundType::NO;
};

class LPProblem
{
public:
    int n = 0;      // число переменных (текущий размер векторов objective/коэффициентов)
    int to_max = 1; // 1 = original as-is, -1 = original was MAX and was negated

    ObjectiveType objective_type = ObjectiveType::MINIMIZE;
    vector<double> objective;
    vector<Constraint> constraints;
    vector<VariableBound> bounds;

    LPProblem();
    LPProblem(int n);
    LPProblem(const LPProblem &other);
    LPProblem(vector<double> objective, ObjectiveType objective_type = ObjectiveType::MINIMIZE);
    virtual ~LPProblem() = default;

    void set_solution_dim(int n);
    void set_objective(const std::vector<double> &coeffs, ObjectiveType type);
    void add_constraint(const Constraint &c);
    void add_var_bound(const VariableBound &vb);
    void print_problem() const;

    std::vector<double> get_objective() const;
    std::vector<Constraint> get_constraints() const;

    virtual std::vector<double> get_initial_solution(std::vector<double> solution) = 0;
    virtual void convert() = 0;
    virtual std::unique_ptr<LPProblem> dual() = 0;
    virtual void to_standard_form() = 0;
};

// Общая (general) форма
class LPProblemGeneral : public LPProblem
{
public:
    std::vector<std::pair<int, int>> bounds_to_subtract; // (original_index_0based, negative_part_index_0based)
    std::vector<int> negated_vars;
    int initial_dim = 0;
    using LPProblem::LPProblem;

    LPProblemGeneral(const LPProblemGeneral &other)
        : LPProblem(other),
          bounds_to_subtract(other.bounds_to_subtract),
          initial_dim(other.initial_dim)
    {
    }

    void convert() override;
    std::vector<double> get_initial_solution(std::vector<double> solution) override;
    std::unique_ptr<LPProblem> dual() override;
    void to_standard_form() override;
};

// Каноническая (slack) форма
class LPProblemSlack : public LPProblemGeneral
{
public:
    using LPProblemGeneral::bounds_to_subtract;
    using LPProblemGeneral::initial_dim;

    // Конструктор копирования — НЕ копируем локальное поле (его больше нет!)
    LPProblemSlack(const LPProblemSlack &other)
        : LPProblemGeneral(other) // базовый класс скопирует bounds_to_subtract
    {
        // initial_dim уже скопирован через базовый класс
    }

    using LPProblemGeneral::LPProblemGeneral;

    void convert() override;
    std::vector<double> get_initial_solution(std::vector<double> solution) override;
    std::unique_ptr<LPProblem> dual() override;
};
