#include "include/TaskLoader.h"
#include <iostream>
#include <memory>

void run_interface()
{
    std::string filename;
    std::cout << "Enter input file name: ";
    std::cin >> filename;

    std::unique_ptr<LPProblem> problem;

    try
    {
        problem = LPProblemParser::parse(filename);
        std::cout << "Parsed problem:" << std::endl;
        problem->print_problem();
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return;
    }

    while (true)
    {
        int opt;
        std::cout << "Choose option:\n"
                     "1 - print problem\n"
                     "2 - convert problem to canonical form\n"
                     "3 - get dual problem (and set as main problem)\n"
                     "4 - get standart problem\n"
                     "5 - solve problem without step-by-step logs\n"
                     "6 - solve with enumerate vertices\n"
                     "7 - exit\n"
                  << std::endl;
        std::cin >> opt;

        switch (opt)
        {
        case 1:
            problem->print_problem();
            break;

        case 2:
            problem->convert();
            std::cout << "Converted to general/canonical (modified) form:\n";
            problem->print_problem();
            break;

        case 3:
            problem = problem->dual();
            std::cout << "Replaced by dual problem:\n";
            problem->print_problem();
            break;

        case 4:
            problem->to_standard_form();
            std::cout << "Replaced by dual problem:\n";
            problem->print_problem();
            break;

        case 5:
        {
            problem->convert();
            SimplexSolver solver;
            bool step_by_step = (opt == 5);
            LPProblemSolution solution = solver.solve(*problem, step_by_step);

            std::cout << "Solution to canonical problem:" << std::endl;
            solution.print_sol();

            if (solution.status == Status::OPTIMAL)
            {
                std::vector<double> init_sol = problem->get_initial_solution(solution.solution);
                solution.set_solution(init_sol);

                if (problem->to_max == -1)
                {
                    solution.set_objective(solution.objective_value * -1);
                }

                std::cout << "\n\nSolution to initial problem:" << std::endl;
                solution.print_sol();
            }
            break;
        }

        case 6:
        {
            SimplexSolver solver;
            auto vertices = solver.enumerate_vertices(*problem, true);
            for (size_t i = 0; i < vertices.size(); ++i)
            {
                std::cout << "Vertex " << (i + 1) << ": ";
                for (double val : vertices[i])
                    std::cout << val << " ";
                std::cout << "\n";
            }

            // Один цикл для поиска оптимума с учётом типа задачи (мин/макс)
            size_t best_idx = 0;
            double best_val = 0.0;
            for (size_t j = 0; j < problem->objective.size(); ++j)
                best_val += problem->objective[j] * vertices[0][j];

            bool is_minimize = (problem->objective_type == ObjectiveType::MINIMIZE);

            for (size_t i = 1; i < vertices.size(); ++i)
            {
                double val = 0.0;
                for (size_t j = 0; j < problem->objective.size(); ++j)
                    val += problem->objective[j] * vertices[i][j];

                // КЛЮЧЕВОЕ ИСПРАВЛЕНИЕ: выбор минимума или максимума в зависимости от типа задачи
                if ((is_minimize && val < best_val) || (!is_minimize && val > best_val))
                {
                    best_val = val;
                    best_idx = i;
                }
            }

            std::cout << "\nOptimal vertex: ";
            for (double val : vertices[best_idx])
                std::cout << val << " ";
            std::cout << "\nObjective value: " << best_val << "\n";
            break;
        }

        case 7:
            return;

        default:
            std::cout << "Invalid option, try again." << std::endl;
            break;
        }
    }
}

int main()
{
    run_interface();
    return 0;
}
