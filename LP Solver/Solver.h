#include <armadillo>
#include <exception>

namespace LPSolver
{
    enum Status{
        unique,    //唯一解
        infinite,  //无穷个解
        none,      //无解
        unbounded, //无界解
    };

    Status simplex(
        bool is_min,
        const arma::vec & f, 
        const arma::mat & A,
        const arma::vec & b,
        const arma::mat & Aeq,
        const arma::vec & beq,
        const arma::vec & lb,
        const arma::vec & ub,
        arma::vec& x);

    /*
    normal form：
    max  f'x
    s.t. Ax == b > 0
          x >= 0
    */
    Status SimplexNormalFormWithSlackVariables(
        const arma::vec & f,
        const arma::mat & A,
        const arma::vec & b,
        arma::vec& x);
}