#include "Solver.h"
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
        /*int n_variables = A.n_cols;
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
            throw std::invalid_argument("Dimension mismatch.");*/
    }

    /*
    TODO: 
        两步转换：
           一般的LP问题变为单纯性标准形式；
           单纯性标准形式加上二阶段法使用的人工变量。
        为了减少矩阵复制的次数，希望一步做完。
        注：小于约束的两步转换可以共用同一个变量
    */
    void ConvertToNormalForm(
        bool is_min, 
        const vec & f,
        const mat & A, 
        const vec & b,
        const mat & Aeq, 
        const vec & beq,
        const vec & lb,
        const vec & ub)
    {
        // unimplemented
        throw std::exception();
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
        // unimplemented
        throw std::exception();
        return Status::unique;
    }

    int FindSwapIn(const vec & f)
    {
        for (int i = 0; i < f.n_rows; i++)
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

    template<typename T>
    T vec_remove_elem(const T& v, const uvec& index)
    {
        uvec sorted_index = sort_index(index);
        T result = T(v.n_rows - index.n_rows);
        int i_v = 0, i_result = 0, i_index = 0;
        for (; i_index < index.n_rows; i_v++)
        {
            if (i_v == index[sorted_index[i_index]])
                i_index++;
            else
                result[i_result++] = v[i_v];
        }
        for (; i_v < v.n_rows; i_v++)
            result[i_result++] = v[i_v];
        return result;
    }

    bool isDegenerated(const vec & f, const uvec & basic_index, const uvec & artificial_index)
    {
        vec tmp_f = f;
        tmp_f.elem(basic_index) += 1;
        tmp_f.elem(artificial_index) += 1;
        return any(tmp_f.subvec(0, f.n_rows - 1) > -EPSILON);  //fn+k == 0
    }

    /*
    normal form：
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
        uvec & basic_index,
        const uvec & artificial_index,
        vec & x,
        double & optimun)
    {
        int n_constraints = Ab.n_rows;
        int n_variables = Ab.n_cols - 1;
        int n_original_variables = n_variables - n_constraints;

        //一阶段目标函数
        vec h = sum(Ab, 0).t();
        h.elem(artificial_index).fill(0);
        f.t().print("f");
        h.t().print("h");
        Ab.print("Ab");
        if (artificial_index.n_rows > 0)
        {
            //一阶段
            while (any(h.subvec(0, n_variables - 1) > 0))
            {
                int swap_in = FindSwapIn(h);
                int swap_out = FindSwapOut(Ab, swap_in);
                if (swap_out == -1)
                    return Status::none;
                basic_index[swap_out] = swap_in;
                GaussianElimination(swap_in, swap_out, f, h, Ab);
                f.t().print("f");
                h.t().print("h");
                Ab.print("Ab");
            }
            double first_opt = h[n_variables];
            if (abs(first_opt) > EPSILON)
                return Status::none;
        }
        //TODO:去除人工变量后做2阶段
        //     人工变量可能出现在前面
        //二阶段
        uvec origin_index = linspace<uvec>(0, n_variables - 1, n_variables);
        origin_index = vec_remove_elem(origin_index, artificial_index);
        while (any(f.elem(origin_index) > 0))  //TODO：要改
        {
            int swap_in = FindSwapIn(f);
            int swap_out = FindSwapOut(Ab, swap_in);
            if (swap_out == -1)
                return Status::unbounded;
            basic_index[swap_out] = swap_in;
            GaussianElimination(swap_in, swap_out, f, Ab);
            f.t().print("f");
            Ab.print("Ab");
        }
        x.set_size(n_variables);
        x.fill(0);
        x.elem(basic_index) = Ab.col(n_variables);
        optimun = -f[n_variables];
        if (isDegenerated(f, basic_index, artificial_index))
            return Status::infinite;
        return Status::unique;
    }
}

