#pragma once

#include <vector>
#include "TaskModel.h"
#include "Matrix.h"
#include <cmath>

using std::vector;

enum Status
{
    OPTIMAL,
    INFEASABLE,
    UNBOUNDED
};

class LPProblemSolution
{
public:
    Status status;
    vector<double> solution;
    double objective_value;
    LPProblemSolution(Status status, vector<double> solution = {}, double objective_value = 0)
    {
        this->status = status;
        if (status == Status::OPTIMAL)
        {
            this->solution = solution;
            this->objective_value = objective_value;
        }
        else if (status == Status::INFEASABLE)
        {
            this->solution = {};
            this->objective_value = 0;
        }
        else if (status == Status::UNBOUNDED)
        {
            this->solution = {};
            this->objective_value = INFINITY;
        }
    }
    void set_solution(vector<double> solution);
    void set_objective(double objective);
    void print_sol();
};

class LPProblemSolver
{
public:
    virtual LPProblemSolution &solve(LPProblem &problem, bool logs = false, vector<double> support = vector<double>(0)) = 0;
};

class SimplexSolver : public LPProblemSolver
{
private:
    vector<double> find_initial_basic_solution(vector<Constraint> &);

public:
    LPProblemSolution &solve(LPProblem &problem, bool logs = false, vector<double> support = vector<double>(0)) override;
};