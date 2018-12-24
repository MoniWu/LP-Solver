#include <armadillo>
#include <exception>

namespace LPSolver
{
    enum Status{
        unique,    //Ψһ��
        infinite,  //�������
        none,      //�޽�
        unbounded, //�޽��
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
    normal form��
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