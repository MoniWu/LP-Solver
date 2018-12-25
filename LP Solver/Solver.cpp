#include "Solver.h"
#include <tuple>

namespace LPSolver
{
    using namespace arma;

    void check_dimension(
        const vec & f,
        const mat & A,
        const vec & b,
        const mat & Aeq,
        const vec & beq,
        const vec & lb,
        const vec & ub)
    {
        int n_variables = A.n_cols;
        int n_equ_constraints = Aeq.n_rows;
        int n_inequ_constraints = A.n_rows;

        if (f.n_rows != n_variables)
            throw std::invalid_argument("Dimension mismatch.");
        if (b.n_rows != n_inequ_constraints)
            throw std::invalid_argument("Dimension mismatch.");
        if (beq.n_rows != n_equ_constraints)
            throw std::invalid_argument("Dimension mismatch.");
        if (lb.n_rows != n_variables)
            throw std::invalid_argument("Dimension mismatch.");
        if (ub.n_rows != n_variables)
            throw std::invalid_argument("Dimension mismatch.");
        if (Aeq.n_cols != n_variables)
            throw std::invalid_argument("Dimension mismatch.");
    }

    /*
    返回极大化问题的
    */
    std::tuple<mat, mat, mat> normal_form(
        bool is_min, 
        const vec & f,
        const mat & A, 
        const vec & b,
        const mat & Aeq, 
        const vec & beq,
        const vec & lb,
        const vec & ub)
    {
        int n_variables = A.n_cols;
        int n_equ_constraints = Aeq.n_rows;
        int n_inequ_constraints = A.n_rows;

        vec normal_f = vec(3 * n_variables + n_equ_constraints + n_inequ_constraints, 0);
        normal_f.head(n_variables) = is_min ? -f : f;


        return std::tuple<mat, mat, mat>();
    }

    Status first_phase(
        const vec & f,
        const mat & A,
        const vec & b,
        vec& x)
    {
        return Status::infinite;
    }

    Status simplex(
        bool is_min,
        const vec & f,
        const mat & A,
        const vec & b,
        const mat & Aeq,
        const vec & beq,
        const vec & lb,
        const vec & ub,
        vec & x)
    {
        check_dimension(f, A, b, Aeq, beq, lb, ub);
        

        return Status::unique;
    }

    /*
    normal form：
    max  (f, 0)'(x, x_bar)
    s.t. (A, I)(x, x_bar)' = b > 0
          x >= 0
    */
    Status SimplexNormalFormWithSlackVariables(
        vec & f,
        mat & A,
        vec & b,
        vec & x,
        double & optimun)
    {
        int n_constraints = A.n_rows;
        int n_variables = A.n_cols;
        int n_original_variables = n_variables - n_constraints;

        /*vec new_f = vec(n_variables + n_constraints);
        new_f.tail(n_constraints).fill(0);
        new_f.head(n_variables) = f;

        uvec pos = b < 0;
        vec new_b = vec(b);
        new_b.elem(pos) *= -1;

        mat new_A = mat(n_constraints, n_variables + n_constraints);
        new_A.submat(0, 0, n_constraints - 1, n_variables - 1) = A;
        new_A.rows(pos) *= -1;
        new_A.submat(0, n_variables, n_constraints - 1, n_variables + n_constraints - 1).eye();*/

        //一阶段目标函数
        vec h = sum(A, 0).t();
        h.tail(n_constraints).fill(0);
        double first_opt = sum(b);
        double second_opt = 0;
        //一阶段
        while (any(h > 0))
        {
            int swap_out = (h > 0).index_max();
            auto positive_index = find(A.col(swap_out) > 0);
            int swap_in = (b.elem(positive_index) / A.col(swap_out)(positive_index));
        }

        return Status::infinite;
    }
}

