#include <iostream>
#include <fstream>
#include <sstream>
#include <armadillo>
#include "Solver.h"
#include "Branch.h"

using namespace arma;
using namespace LPSolver;
using namespace std;


void read(int& cnt, string& f, vector<string>& cst) 
{
	ifstream in("C:\\Users\\lenovo\\Desktop\\LP\\case3.txt");
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
	int cnt = 0;
	string ff;
	vector<string> cst;
	read(cnt, ff, cst);

	vec f = vec(cnt).fill(0);
	mat Ab = mat(cst.size(),cnt+1).fill(0);
	vector<int> mode;
	Branch::parse(cnt, ff, cst, f, Ab, mode);

	uvec base;
	uvec arti;
	Branch::normal(f, Ab, mode, base, arti);
	f.t().print();
	Ab.print();
	base.t().print();

	Branch bh = Branch(f,Ab,base,arti);

	vec x = vec();
	double opt;
	Status status;
	status = SimplexNormalFormWithSlackVariables(bh.f, bh.Ab, bh.base, bh.arti, x, opt);
	cout << opt << endl;
	cout << status << endl;
	x.t().print();
	bh.Ab.print();
	cout << "new branch:" << endl;

	Branch low = bh.upBranch(0, x);
	cout << "after normal: " << endl;
	low.Ab.print();
	cout << "base: ";
	low.base.t().print();

	vec x2 = vec();
	Status s2 = SimplexNormalFormWithSlackVariables(low.f, low.Ab, low.base, low.arti, x2, opt);
	cout << "Branch result:" << endl;
	cout << "opt: " << opt << endl;
	cout << "x: ";
	x2.t().print();
	low.Ab.print();
	getchar();
	return 0;
}
