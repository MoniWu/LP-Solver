#include <iostream>
#include <armadillo>
#include "Solver.h"

using namespace arma;
using namespace LPSolver;
using namespace std;

int main()
{
    vec f = vec({ 0,0.1,0.2,0.3,0.8,0,0,0,0 });
    f.t().print("\n");
    mat Ab = mat({ {1,2,0,1,0,1,0,0,100},{0,0,2,2,1,0,1,0,100},{3,1,2,0,3,0,0,1,100} });
    Ab.print("\n");
    vec x = vec();
    double opt;
    Status status;
    status = SimplexNormalFormWithSlackVariables(f, Ab, x, opt);
    cout << opt << endl;
    cout << status << endl;
    x.t().print("\n");
    /*vec f = vec({ -3,1,1,0,0,0,0,0,0 });
    f.t().print("\n");
    mat Ab = mat({ {1,-2,1,1,0,1,0,0,11},{-4,1,2,0,-1,0,1,0,3},{-2,0,1,0,0,0,0,1,1} });
    Ab.print("\n");
    vec x = vec();
    double opt;
    Status status;
    status = SimplexNormalFormWithSlackVariables(f, Ab, x, opt);
    cout << opt << endl;
    cout << status << endl;
    x.t().print("\n");*/
}
