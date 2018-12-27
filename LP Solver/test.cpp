#include <iostream>
#include <armadillo>
using namespace arma;
int main()
{
    //生成一个随机矩阵A,大小为5x5，矩阵每个元素的范围为：(0,10)
    mat A = randu<mat>(5,2) * 10;
    uvec x = linspace<uvec>(2, 5, 4);
    x.print();
    (x > 0).print();
    x.elem(find(x > 0)).print();
    std::cout << (x > 0).index_max();
    int n_variables = A.n_cols;
    int n_constraints = A.n_rows;

    A.print();
    arma::mat new_A = arma::mat(n_constraints, n_variables + n_constraints);
    new_A.submat(0, 0, n_constraints - 1, n_variables - 1) = A;
    new_A.submat(0, n_variables, n_constraints - 1, n_variables + n_constraints - 1).eye();
    new_A.print();
    new_A.rows(uvec({ 2,4 })) *= -1;
    sum(new_A, 0).print();
}
