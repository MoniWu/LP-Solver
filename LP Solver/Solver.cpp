#include "Solver.h"
#include <tuple>
#include <limits>

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
    ���ؼ��������
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

    Status Simplex(
        bool is_min,
        const vec & f,
        const mat & A,
        const vec & b,
        const mat & Aeq,
        const vec & beq,
        const vec & lb,
        const vec & ub,
        vec & x,
        double & optimun)
    {
        check_dimension(f, A, b, Aeq, beq, lb, ub);
        

        return Status::unique;
    }

    int FindSwapIn(const vec & f)
    {
        for (int i = 0; i < f.n_cols; i++)
        {
            if (f[i] > 0)
            {
                return i;
            }
        }
        return -1;
    }

    int FindSwapOut(const mat & Ab, int swap_in)
    {
        int rows = Ab.n_rows;
        int cols = Ab.n_cols;
        auto b = Ab.col(cols - 1);
        double min_ratio = std::numeric_limits<double>::max();
        int index = -1;
        for (int i = 0; i < rows; i++)
        {
            if (Ab.at(i, swap_in) > 0)
            {
                double ratio = b[i] / Ab.at(i, swap_in);
                if (ratio < min_ratio)
                {
                    min_ratio = ratio;
                    index = i;
                }
            }
        }
        return index;
    }

    void GaussianElimination(int swap_in, int swap_out, vec & f, mat & Ab)
    {
        Ab.row(swap_out) /= Ab.at(swap_out, swap_in);
        for (int i = 0; i < Ab.n_rows; i++)
        {
            if (i == swap_out)continue;
            Ab.row(i) -= Ab.row(i)[swap_in] * Ab.row(swap_out);
        }
        f -= f[swap_in] * Ab.row(swap_out).t();
    }

    void GaussianElimination(int swap_in, int swap_out, vec & f1, vec & f2, mat & Ab)
    {
        GaussianElimination(swap_in, swap_out, f1, Ab);
        f2 -= f2[swap_in] * Ab.row(swap_out).t();
    }

    /*
    normal form��
    max  (f, 0)'(x, x_bar)
    s.t. (A, I)(x, x_bar)' = b > 0
          x >= 0

    input:
    f: [f0, ..., fn+m-1, opt]
    Ab: [ A0,0 , ...,  A0,n+m-1 ,  b0 ]
        [ ...  ,    ,    ...    ,  ...]
        [An-1,0, ..., An-1,n+m-1, bm-1]
    output:
    x: [x0, ..., xn+m-1]
    optimun: (double)

      f0  ,   f1  , ...,   fn+m-1  | opt
    --------------------------------------
     A0,0 ,  A0,1 , ...,  A0,n+m-1 | b0
      ... ,  ...  ,    ,    ...    | ...
    An-1,0, An-1,1, ..., An-1,n+m-1| bm-1
    */
    Status SimplexNormalFormWithSlackVariables(
        vec & f,
        mat & Ab,
        //vec & b,
        vec & x,
        double & optimun)
    {
        int n_constraints = Ab.n_rows;
        int n_variables = Ab.n_cols - 1;
        int n_original_variables = n_variables - n_constraints;

        //һ�׶�Ŀ�꺯��
        vec h = sum(Ab, 0).t();
        h.subvec(n_original_variables, n_variables - 1).fill(0);
        uvec basic_x_subscript = linspace<uvec>(n_original_variables, n_variables - 1, n_constraints);
        //һ�׶�
        while (any(h.subvec(0, h.n_rows - 2) > 0))
        {
            int swap_in = FindSwapIn(h);
            int swap_out = FindSwapOut(Ab, swap_in);
            if (swap_out == -1)
                return Status::none;
            basic_x_subscript[swap_out] = swap_in;
            GaussianElimination(swap_in, swap_out, f, h, Ab);
        }
        double first_opt = h[n_variables];
        if (abs(first_opt) > EPSILON)
            return Status::none;
        //���׶�
        while (any(f.subvec(0, f.n_rows - 2) > 0))
        {
            int swap_in = FindSwapIn(f);
            int swap_out = FindSwapOut(Ab, swap_in);
            if (swap_out == -1)
                return Status::unbounded;
            basic_x_subscript[swap_out] = swap_in;
            GaussianElimination(swap_in, swap_out, f, Ab);
        }
        x.set_size(n_variables);
        x.fill(0);
        x.elem(basic_x_subscript) = Ab.col(n_variables);
        optimun = -f[n_variables];
        if (any(f.subvec(n_original_variables, n_variables - 1) > -EPSILON))  //fn+k == 0
            return Status::infinite;
        return Status::unique;
    }

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
}
