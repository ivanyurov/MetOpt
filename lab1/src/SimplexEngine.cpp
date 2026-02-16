#include <vector>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include "../include/SimplexEngine.h"
#include <cmath>
#include <iomanip>
#include <string>
#include <ostream>
#include <fstream>

using std::pair;
using std::vector;

void print_vector_file(const std::vector<int> &a, std::ofstream &file)
{
    for (auto &j : a)
        file << j << " ";
    file << "\n\n";
}

void print_matrix_file(Matrix &M, std::ofstream &file)
{
    file << "\n";
    for (auto &line_ind : M.line_indexes)
    {
        for (auto &col_ind : M.column_indexes)
        {
            file << M.matrix[line_ind][col_ind] << " ";
        }
        file << "\n";
    }
    file << "\n";
}

void print_simplex_table(std::ostream &out,
                         Matrix A,
                         Matrix B,
                         Matrix X,
                         Matrix c,
                         Matrix dkt,
                         std::vector<int> n_all,
                         std::vector<int> nk)
{
    out << "\n   SIMPLEX TABLE\n\n";

    // Restore all columns in A
    A.set_columns(n_all);

    // Compute tableau coefficients: B * A
    Matrix tableau_coeff = B.multiply(A);

    // RHS values from X for basic variables
    std::vector<double> rhs;
    for (int i : nk)
    {
        rhs.push_back(X.matrix[i][0]);
    }

    // Determine column widths
    std::vector<size_t> column_widths(tableau_coeff.column_size() + 2, 0);
    // +2 → Basis column + RHS column

    // Basis column width
    column_widths[0] = 7; // enough for "x123"

    // Variable columns width
    for (int j = 0; j < tableau_coeff.column_size(); ++j)
    {
        column_widths[j + 1] = 6; // default minimum width

        for (size_t i = 0; i < tableau_coeff.matrix.size(); ++i)
        {
            size_t len = std::to_string(tableau_coeff.matrix[i][j]).length();
            if (len > column_widths[j + 1])
                column_widths[j + 1] = len;
        }
    }

    // RHS width
    column_widths.back() = 6;
    for (double value : rhs)
    {
        size_t len = std::to_string(value).length();
        if (len > column_widths.back())
            column_widths.back() = len;
    }

    out << std::setw(column_widths[0]) << "Basis" << " | ";
    for (int j = 0; j < tableau_coeff.column_size(); ++j)
    {
        out << std::setw(column_widths[j + 1])
            << ("x" + std::to_string(j + 1))
            << " | ";
    }
    out << std::setw(column_widths.back()) << "RHS" << "\n";

    for (size_t width : column_widths)
    {
        out << std::string(width + 3, '-') << "+";
    }
    out << "\n";

    for (size_t i = 0; i < nk.size(); ++i)
    {
        out << std::setw(column_widths[0])
            << ("x" + std::to_string(nk[i] + 1)) << " | ";

        for (int j = 0; j < tableau_coeff.column_size(); ++j)
        {
            out << std::setw(column_widths[j + 1])
                << std::fixed << std::setprecision(4)
                << tableau_coeff.matrix[i][j] << " | ";
        }

        out << std::setw(column_widths.back())
            << std::fixed << std::setprecision(4)
            << rhs[i] << "\n";
    }

    double obj_value = c.transpose().multiply(X).matrix[0][0];

    out << std::setw(column_widths[0]) << "z" << " | ";
    for (int j = 0; j < dkt.matrix[0].size(); ++j)
    {
        out << std::setw(column_widths[j + 1])
            << std::fixed << std::setprecision(4)
            << dkt.matrix[0][j] << " | ";
    }
    out << "\n";

    out << std::setw(column_widths.back())
        << std::fixed << std::setprecision(6)
        << obj_value << "\n";

    for (size_t width : column_widths)
    {
        out << std::string(width + 3, '-') << "+";
    }
    out << "\n\n";
}

// Find initial basic feasible solution using the artificial basis method
vector<double> SimplexSolver::find_initial_basic_solution(vector<Constraint> &constraints)
{
    int m = constraints.size();
    int n = constraints[0].coefficients.size();

    vector<vector<double>> A(m, vector<double>(n + m, 0));
    vector<double> b(m, 0);
    vector<double> c(n + m, 0);

    vector<Constraint> tmp(constraints);

    vector<double> solution(n + m, 0);
    for (int i = 0; i < m; i++)
    {

        for (int j = 0; j < n; j++)
        {
            A[i][j] = constraints[i].coefficients[j];
        }
        b[i] = constraints[i].b;
        if (b[i] < 0)
        {
            b[i] *= -1;
            for (int j = 0; j < n; j++)
            {
                A[i][j] *= -1;
            }
            A[i][n + i] = -1;
        }
        else
        {
            A[i][n + i] = 1;
        }
        c[n + i] = 1;
        solution[n + i] = b[i];
    }

    LPProblem *problem = new LPProblemSlack(m + n);
    problem->set_objective(c, ObjectiveType::MINIMIZE);
    for (int i = 0; i < m; i++)
    {

        Constraint constraint;
        constraint.coefficients = A[i];
        constraint.b = b[i];
        constraint.type = InequalityType::EQUAL;
        problem->add_constraint(constraint);
    }

    for (int i = 0; i < n + m; i++)
    {
        problem->add_var_bound({i + 1, BoundType::NOT_NEGATIVE});
    }
    problem->print_problem();
    LPProblemSolution artificial_solution = solve(*problem, false, solution);

    bool infeasible = false;
    for (int i = n; i < n + m; ++i)
    {
        if (artificial_solution.solution[i] < 0)
        {
            infeasible = true;
            break;
        }
    }
    if (infeasible)
    {
        throw std::runtime_error("Problem is infeasible");
    }

    vector<double> ans;
    for (int i = 0; i < n; i++)
    {
        ans.push_back(artificial_solution.solution[i]);
    }

    return ans;
}

LPProblemSolution &SimplexSolver::solve(LPProblem &problem, bool logs, vector<double> support)
{
    std::ofstream logfile;

    if (logs)
    {
        logfile.open("log.txt", std::ios::out);
        logfile << "        Simplex method full log\n\n";
        std::cout << "Simplex method started.\n";
    }

    if (support.empty())
    {
        auto constraints = problem.get_constraints();
        support = find_initial_basic_solution(constraints);
    }

    Matrix X(support);

    if (logs)
    {
        logfile << "Initial support vector:\n";
        print_matrix_file(X, logfile);
    }

    // Build full constraint matrix A_ from problem (never mutate this one)
    vector<vector<double>> A_;
    for (auto &constr : problem.get_constraints())
    {
        A_.push_back(constr.coefficients);
    }

    // ORIGINAL, immutable copy of A used to create fresh temporaries
    Matrix A_original(A_);

    if (logs)
    {
        logfile << "Constraint matrix A (original):\n";
        print_matrix_file(A_original, logfile);
    }

    // Basic validations (use original)
    if (A_original.matrix.size() > A_original.matrix[0].size())
        throw std::runtime_error("The number of rows in matrix A must be <= number of columns.");

    if (!(A_original.is_full_rank()))
        throw std::runtime_error("Matrix A must be of full rank.");

    int iter = 1;
    vector<int> nk;
    vector<int> lk;
    vector<pair<int, int>> used_in_basis_change;
    bool basis_changed = false;

    while (true)
    {
        if (!basis_changed)
            used_in_basis_change.clear();

        if (logs)
        {
            std::cout << "\n--- Iteration " << iter << " ---\n";
            logfile << "\n--- Iteration " << iter << " ---\n";
        }

        // Build index sets from current X (X is the vector of basic variable values)
        vector<int> n_plus, n_zero, n_all;
        for (int i = 0; i < static_cast<int>(X.matrix.size()); ++i)
        {
            if (X.matrix[i][0] != 0)
                n_plus.push_back(i);
            else
                n_zero.push_back(i);

            n_all.push_back(i);
        }

        // Work on a copy of A when we need to select columns
        Matrix A_for_nplus = A_original;
        A_for_nplus.set_columns(n_plus);

        // Find basis columns if needed — operate on a fresh copy
        if (nk.empty() || !basis_changed)
        {
            // get_addition_to_square_matrix expects available indices relative to the full A_for_nplus
            nk = A_for_nplus.get_addition_to_square_matrix(n_zero);
        }
        else
        {
            // ensure nk is applied to a fresh matrix (so subsequent uses start from original)
            Matrix A_for_nk = A_original;
            A_for_nk.set_columns(nk);
            // we don't need to store it here, just ensure nk is valid
        }

        lk = subtract_vectors<int>(n_all, nk);

        if (logs)
        {
            std::cout << "Current basis: ";
            print_vector(nk);

            logfile << "n_plus:\n";
            print_vector_file(n_plus, logfile);

            logfile << "n_zero:\n";
            print_vector_file(n_zero, logfile);

            logfile << "nk (basis):\n";
            print_vector_file(nk, logfile);

            logfile << "lk (non-basis):\n";
            print_vector_file(lk, logfile);
        }

        // c is full objective vector (as a column)
        Matrix c(problem.get_objective());

        // Build B = inverse of basis matrix (we must build basis matrix from ORIGINAL)
        Matrix A_basis = A_original;
        A_basis.set_columns(nk);
        Matrix B = A_basis.get_inverse_matrix();

        if (logs)
        {
            logfile << "Basis inverse matrix B:\n";
            print_matrix_file(B, logfile);
        }

        // cnk is c restricted to basis columns: use a fresh allocate on original
        Matrix c_nk = c.allocate_matrix(nk, {0}); // column vector of c for basis

        // For reduced costs and tableau we need A with all columns (full A), but B.multiply(A_full)
        Matrix A_full = A_original;    // fresh copy
        A_full.set_columns(n_all);     // usually n_all is full set; safe to set

        // d_k^T = c^T - c_Nk^T * B * A
        Matrix dkt = c.transpose().subtract(c_nk.transpose().multiply(B.multiply(A_full)));

        if (logs)
        {
            logfile << "Reduced costs d_k:\n";
            print_matrix_file(dkt, logfile);

            print_simplex_table(logfile, A_full, B, X, c, dkt, n_all, nk);
        }

        // take reduced costs for non-basis columns
        Matrix dklkt = dkt.allocate_matrix({0}, lk);

        if (logs)
            std::cout << "Checking optimality...\n";

        if (dklkt >= 0)
        {
            if (logs)
            {
                std::cout << "Optimal solution found.\n";
                logfile << "All reduced costs >= 0.\n";
                logfile << "Optimal solution reached.\n";
                logfile << "\n      End of log\n";
                logfile.close();
            }

            LPProblemSolution *solution =
                new LPProblemSolution(Status::OPTIMAL,
                                      X.transpose().matrix[0],
                                      c.transpose().multiply(X).matrix[0][0]);
            return *solution;
        }

        if (logs)
            std::cout << "Not optimal. Selecting entering variable...\n";

        // Choose first entering variable with negative reduced cost
        vector<int> jk;
        for (int i = 0; i < static_cast<int>(dkt.matrix[0].size()); ++i)
        {
            if (dkt.matrix[0][i] < 0)
            {
                jk.push_back(i);
                break;
            }
        }

        if (logs)
        {
            std::cout << "Entering variable: ";
            print_vector(jk);

            logfile << "Chosen entering variable:\n";
            print_vector_file(jk, logfile);
        }

        // Build A_jk as a fresh matrix (columns = jk) from original
        Matrix A_jk = A_original;
        A_jk.set_columns(jk);

        // BA = B * A_jk (direction for basis variables)
        Matrix BA = B.multiply(A_jk);

        // Build full uk_ vector of size n_all, fill using BA rows mapped to basis indices
        vector<double> uk_(n_all.size(), 0.0);

        for (int i = 0; i < static_cast<int>(BA.matrix.size()); ++i)
            uk_[nk[i]] = BA.matrix[i][0];

        for (auto &j : jk)
            uk_[j] = -1;

        Matrix uk(uk_);

        if (logs)
        {
            logfile << "Direction vector u_k:\n";
            print_matrix_file(uk, logfile);
            std::cout << "Checking unboundedness...\n";
        }

        // Check if uk[Nk] <= 0 (no positive entries among basis positions)
        Matrix uk_nk = uk.allocate_matrix(nk, {0});
        if (uk_nk <= 0)
        {
            if (logs)
            {
                std::cout << "Problem is unbounded.\n";
                logfile << "All u_k[Nk] <= 0. Problem is unbounded.\n";
                logfile << "\n      End of log\n";
                logfile.close();
            }

            LPProblemSolution *solution =
                new LPProblemSolution(Status::UNBOUNDED);
            return *solution;
        }

        // indices of basis rows with positive uk
        vector<int> i_;
        for (auto &i : nk)
            if (uk.matrix[i][0] > 0)
                i_.push_back(i);

        // check degeneracy / choose theta
        if (nk == n_plus ||
            subtract_vectors<int>(nk, n_plus).empty() ||
            uk.allocate_matrix(subtract_vectors<int>(nk, n_plus), {0}) <= 0)
        {
            double theta_k = 1e18;
            for (auto &i : i_)
            {
                double denom = uk.matrix[i][0];
                if (denom > 0)
                {
                    double cand = X.matrix[i][0] / denom;
                    if (cand < theta_k) theta_k = cand;
                }
            }

            if (logs)
            {
                std::cout << "Theta = " << theta_k << "\n";
                logfile << "Theta = " << theta_k << "\n";
            }

            X = X.subtract(uk * theta_k);
            basis_changed = false;
        }
        else
        {
            if (logs)
                std::cout << "Degeneracy detected. Attempting basis change...\n";

            bool ext = false;
            vector<int> choice = subtract_vectors<int>(nk, n_plus);

            for (int idx_l = 0; idx_l < static_cast<int>(lk.size()) && !ext; ++idx_l)
            {
                int l = lk[idx_l];

                for (int idx_choice = 0; idx_choice < static_cast<int>(choice.size()) && !ext; ++idx_choice)
                {
                    int n_candidate = choice[idx_choice];

                    bool skip = false;
                    for (auto &p : used_in_basis_change)
                        if (p.first == n_candidate && p.second == l)
                            skip = true;

                    if (skip) continue;

                    // Candidate basis: replace n_candidate by l
                    vector<int> columns = concatenate_vectors<int>({l}, subtract_vectors<int>(nk, {n_candidate}));

                    // Use a fresh copy of original matrix to test determinant
                    Matrix A_test = A_original;
                    A_test.set_columns(columns);

                    double det = A_test.determinant(A_test.allocate_matrix(A_test.line_indexes, A_test.column_indexes).matrix);
                    if (std::fabs(det) > 1e-12)
                    {
                        nk = columns;
                        // no in-place mutation of A_original; downstream iterations will recreate copies as needed

                        used_in_basis_change.push_back({l, n_candidate});
                        used_in_basis_change.push_back({n_candidate, l});

                        if (logs)
                        {
                            std::cout << "Basis changed: column "
                                      << n_candidate << " replaced with column "
                                      << l << "\n";

                            logfile << "Basis change performed:\n";
                            logfile << "Replaced column "
                                    << n_candidate << " with column "
                                    << l << "\n";

                            logfile << "New basis:\n";
                            print_vector_file(nk, logfile);
                        }

                        ext = true;
                        basis_changed = true;
                        break;
                    }
                }
            }
        }

        ++iter;
    }

    // unreachable normally, but keep cleanup consistent
    if (logs)
    {
        logfile << "\n      End of log\n";
        logfile.close();
    }

    LPProblemSolution *solution =
        new LPProblemSolution(Status::INFEASABLE);

    return *solution;
}

std::vector<std::vector<double>> SimplexSolver::enumerate_vertices(LPProblem &problem, bool logs) const
{
    if (logs)
        std::cout << "Starting enumeration of extreme points...\n";

    // === БЕЗОПАСНОЕ КОПИРОВАНИЕ: создаём новую задачу в канонической форме ===
    LPProblemSlack copy_problem;
    
    // Копируем исходные данные
    copy_problem.n = problem.n;
    copy_problem.objective = problem.get_objective();
    copy_problem.objective_type = problem.objective_type;
    copy_problem.to_max = problem.to_max;
    copy_problem.constraints = problem.get_constraints();
    copy_problem.bounds = problem.bounds;
    copy_problem.initial_dim = problem.n;
    
    // Применяем преобразование к канонической форме
    copy_problem.convert();

    int m = static_cast<int>(copy_problem.constraints.size());
    int n = copy_problem.n;

    if (logs)
    {
        std::cout << "DEBUG: After convert() — m = " << m << " constraints, n = " << n << " variables\n";
        // Для вашей задачи должно быть: m=4, n=10 → C(10,4)=210
        long long expected = 1;
        for (int i = 1; i <= m; ++i)
            expected = expected * (n - m + i) / i;
        std::cout << "DEBUG: Expected combinations C(" << n << "," << m << ") = " << expected << "\n";
    }

    if (m == 0)
        return {std::vector<double>(n, 0.0)};

    // Формируем матрицу A и вектор b
    std::vector<std::vector<double>> A_mat(m, std::vector<double>(n, 0.0));
    std::vector<double> b_vec(m, 0.0);
    for (int i = 0; i < m; ++i)
    {
        b_vec[i] = copy_problem.constraints[i].b;
        const auto &coeffs = copy_problem.constraints[i].coefficients;
        for (size_t j = 0; j < coeffs.size() && j < (size_t)n; ++j)
            A_mat[i][j] = coeffs[j];
    }
    Matrix A_full(A_mat);

    // Генерируем ВСЕ сочетания из n по m
    std::vector<int> all_columns(n);
    for (int i = 0; i < n; ++i) 
        all_columns[i] = i;

    auto basis_combinations = combinations(all_columns, m);

    //if (logs)
        //std::cout << "DEBUG: Generated " << basis_combinations.size() << " basis combinations\n";

    std::vector<std::vector<double>> vertices;
    int feasible_count = 0;

    for (const auto &basis : basis_combinations)
    {
        // Формируем базисную матрицу B (m x m)
        std::vector<int> row_indices(m);
        for (int i = 0; i < m; ++i) 
            row_indices[i] = i;

        Matrix B = A_full.allocate_matrix(row_indices, basis);

        // Проверяем невырожденность
        double det = B.determinant(B.matrix);
        if (std::abs(det) < 1e-10)
            continue;

        // Решаем B * x_B = b
        Matrix B_inv = B.get_inverse_matrix();
        Matrix b_mat(b_vec);
        Matrix x_b = B_inv.multiply(b_mat);

        // Формируем полное решение (небазисные = 0)
        std::vector<double> x(n, 0.0);
        for (int i = 0; i < m; ++i)
        {
            int var_idx = basis[i];
            if (var_idx < n)
                x[var_idx] = x_b.matrix[i][0];
        }

        // Проверяем допустимость: ВСЕ переменные в расширенном пространстве >= 0
        bool feasible = true;
        for (double val : x)
        {
            if (val < -1e-8)
            {
                feasible = false;
                break;
            }
        }

        if (feasible)
        {
            // Восстанавливаем исходные переменные (до расщепления свободных)
            std::vector<double> restored = copy_problem.get_initial_solution(x);
            vertices.push_back(restored);
            ++feasible_count;
        }
    }

    if (logs)
    {
        std::cout << "Checked " << basis_combinations.size() << " bases.\n";
        std::cout << "Found " << feasible_count << " feasible vertices.\n";
        std::cout << "Enumeration completed.\n";
    }

    // Удаляем дубликаты
    std::sort(vertices.begin(), vertices.end());
    auto last = std::unique(vertices.begin(), vertices.end(),
        [](const std::vector<double>& a, const std::vector<double>& b) {
            if (a.size() != b.size()) return false;
            for (size_t i = 0; i < a.size(); ++i)
                if (std::abs(a[i] - b[i]) > 1e-6) 
                    return false;
            return true;
        });
    vertices.erase(last, vertices.end());

    return vertices;
}


void LPProblemSolution::set_solution(vector<double> solution)
{
    this->solution = solution;
}
void LPProblemSolution::set_objective(double objective)
{
    objective_value = objective;
}

void LPProblemSolution::print_sol()
{
    std::cout << "STATUS: " << (status == Status::OPTIMAL ? "Optimal" : status == Status::INFEASABLE ? "Infeasible"
                                                                                                     : "Unbounded")
              << std::endl;
    if (!solution.empty())
    {
        std::cout << "X* = ( ";
        for (auto &el : solution)
        {
            std::cout << el << " ";
        }
        std::cout << ")" << std::endl;
        std::cout << "Objective value is: " << objective_value << std::endl;
    }
}