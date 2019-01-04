#include <iostream>
#include <armadillo>
#include "Solver.h"

using namespace arma;
using namespace LPSolver;
using namespace std;

int main()
{
    /*vec f = vec({ 2,3,0,0,0,0 });
    f.t().print("\n");
    mat Ab = mat({ {1,2,1,0,0,8},{4,0,0,1,0,16},{0,4,0,0,1,12} });
    Ab.print("\n");
    uvec basic_index = uvec({ 2,3,4 });
    vec x = vec();
    double opt;
    Status status;
    status = SimplexNormalFormWithSlackVariables(f, Ab, basic_index, uvec({}), x, opt);
    cout << opt << endl;
    cout << status << endl;
    x.t().print("\n");*/
    /*vec f = vec({ 3,-1,-1,0,0,0,0,0 });
    f.t().print("\n");
    mat Ab = mat({ {1,-2,1,1,0,0,0,11},{-4,1,2,0,-1,1,0,3},{-2,0,1,0,0,0,1,1} });
    Ab.print("\n");
    uvec basic_index = uvec({ 3,5,6 });
    vec x = vec();
    double opt;
    Status status;
    status = SimplexNormalFormWithSlackVariables(f, Ab, basic_index, uvec({ 5,6 }), x, opt);
    cout << opt << endl;
    cout << status << endl;
    x.t().print("\n");*/
    vec f = vec({ 0,3,0,-1,-1,0,0,0 });
    f.t().print("\n");
    mat Ab = mat({ {0,1,0,-2,1,1,0,11},{1,-4,0,1,2,0,-1,3},{0,-2,1,0,1,0,0,1} });
    Ab.print("\n");
    uvec basic_index = uvec({ 5,0,2 });
    vec x = vec();
    double opt;
    Status status;
    status = SimplexNormalFormWithSlackVariables(f, Ab, basic_index, uvec({ 0,2 }), x, opt);
    cout << opt << endl;
    cout << status << endl;
    x.t().print("\n");
    //sort_index(vec({ 10,1,5 })).print();
    //linspace<uvec>(0, 4,5).print();
}
