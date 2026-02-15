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

    vector<vector<double>> A_;
    for (auto &constr : problem.get_constraints())
    {
        A_.push_back(constr.coefficients);
    }

    Matrix A(A_);

    if (logs)
    {
        logfile << "Constraint matrix A:\n";
        print_matrix_file(A, logfile);
    }

    if (A.matrix.size() > A.matrix[0].size())
        throw std::runtime_error("The number of rows in matrix A must be <= number of columns.");

    if (!(A.is_full_rank()))
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

        vector<int> n_plus, n_zero, n_all;
        for (int i = 0; i < X.matrix.size(); i++)
        {
            if (X.matrix[i][0] != 0)
                n_plus.push_back(i);
            else
                n_zero.push_back(i);

            n_all.push_back(i);
        }

        A.set_columns(n_plus);

        if (nk.empty() || !basis_changed)
            nk = A.get_addition_to_square_matrix(n_zero);
        else
            A.set_columns(nk);

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

        Matrix c(problem.get_objective());

        Matrix B = A.get_inverse_matrix();

        if (logs)
        {
            logfile << "Basis inverse matrix B:\n";
            print_matrix_file(B, logfile);
        }

        Matrix cnk = c.allocate_matrix(nk, {0});

        A.set_columns(n_all);

        Matrix dkt = c.transpose().subtract(cnk.transpose().multiply(B.multiply(A)));

        if (logs)
        {
            logfile << "Reduced costs d_k:\n";
            print_matrix_file(dkt, logfile);

            print_simplex_table(logfile, A, B, X, c, dkt, n_all, nk);
        }

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

        vector<int> jk;
        for (int i = 0; i < dkt.matrix[0].size(); i++)
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

        A.set_columns(jk);

        Matrix BA = B.multiply(A);

        vector<double> uk_(n_all.size(), 0);

        for (int i = 0; i < BA.matrix.size(); i++)
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

        if (uk.allocate_matrix(nk, {0}) <= 0)
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

        vector<int> i_;
        for (auto &i : nk)
            if (uk.matrix[i][0] > 0)
                i_.push_back(i);

        if (nk == n_plus ||
            subtract_vectors<int>(nk, n_plus).empty() ||
            uk.allocate_matrix(subtract_vectors<int>(nk, n_plus), {0}) <= 0)
        {
            double theta_k = 1e18;

            for (auto &i : i_)
                if (theta_k > X.matrix[i][0] / uk.matrix[i][0])
                    theta_k = X.matrix[i][0] / uk.matrix[i][0];

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

            for (int i = 0; i < lk.size(); i++)
            {
                int l = lk[i];

                for (int j = 0; j < choice.size(); j++)
                {
                    int n = choice[j];

                    bool skip = false;
                    for (auto &p : used_in_basis_change)
                        if (p.first == n && p.second == l)
                            skip = true;

                    if (skip)
                        continue;

                    vector<int> columns =
                        concatenate_vectors<int>({l},
                                                 subtract_vectors<int>(nk, {n}));

                    A.set_columns(columns);

                    if (A.determinant(
                            A.allocate_matrix(A.line_indexes,
                                              A.column_indexes)
                                .matrix) != 0)
                    {
                        nk = columns;
                        A.set_columns(nk);

                        used_in_basis_change.push_back({l, n});
                        used_in_basis_change.push_back({n, l});

                        if (logs)
                        {
                            std::cout << "Basis changed: column "
                                      << n << " replaced with column "
                                      << l << "\n";

                            logfile << "Basis change performed:\n";
                            logfile << "Replaced column "
                                    << n << " with column "
                                    << l << "\n";

                            logfile << "New basis:\n";
                            print_vector_file(nk, logfile);
                        }

                        ext = true;
                        basis_changed = true;
                        break;
                    }
                }
                if (ext)
                    break;
            }
        }

        iter++;
    }

    if (logs)
    {
        logfile << "\n      End of log\n";
        logfile.close();
    }

    LPProblemSolution *solution =
        new LPProblemSolution(Status::INFEASABLE);

    return *solution;
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