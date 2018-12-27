#include <armadillo>
#include <exception>

namespace LPSolver
{
    const double EPSILON = 1e-7;

    enum Status{
        unique,    //唯一解
        infinite,  //无穷个解
        none,      //无可行解
        unbounded, //无界解
    };

    Status Simplex(
        bool is_min,
        const arma::vec & f, 
        const arma::mat & A,
        const arma::vec & b,
        const arma::mat & Aeq,
        const arma::vec & beq,
        const arma::vec & lb,
        const arma::vec & ub,
        arma::vec& x,
        double & optimun);

    Status SimplexNormalFormWithSlackVariables(
        arma::vec & f,
        arma::mat & Ab,
        arma::vec & x,
        double & optimun);
}