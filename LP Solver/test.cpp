#include <iostream>
#include <fstream>
#include <sstream>
#include <armadillo>
#include "Solver.h"
#include "Branch.h"

#define MLEN 200
using namespace arma;
using namespace LPSolver;
using namespace std;


void read(int& cnt, string& f, vector<string>& cst) 
{
	ifstream in("C:\\Users\\lenovo\\Desktop\\LP\\case0.txt");
	stringstream ss;
	if (in.fail())
	{
		cout << "fail to open" << endl;
		return;
	}
	string buf;
	getline(in, buf);
	ss << buf;
	ss >> cnt;
	getline(in, f);
	while (!in.eof()) {
		getline(in, buf);
		cst.push_back(buf);
	};
	in.close();
}

int main()
{
    /*vec f = vec({2,3,0,0,0,0,0,0,0 });
    f.t().print("\n");
    mat Ab = mat({ {1,2,1,0,0,1,0,0,8},{4,0,0,1,0,0,1,0,16},{0,4,0,0,1,0,0,1,12} });
    Ab.print("\n");
    vec x = vec();
    double opt;
    Status status;
    status = SimplexNormalFormWithSlackVariables(f, Ab, x, opt);
    cout << opt << endl;
    cout << status << endl;
    x.t().print("\n");*/
    /*vec f = vec({ 3,-1,-1,0,0,0,0,0,0 });
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
