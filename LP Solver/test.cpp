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

int RecurBranch(const Branch& bch, double& opt, vec& opt_x, string mark);

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

	Branch root = Branch(f,Ab,base,arti);

	vec x_opt = vec();
	vec x = vec();
	double opt = 0;
	double result = 0;
	
	Status s0 = SimplexNormalFormWithSlackVariables(root.f, root.Ab, root.base, root.arti, x, result);
	int index = 0;
	switch (s0) {
	case Status::infinite: 
		index = Branch::FindFirstNotInt(x);
		if (index == -1) {
			cout << "\t[opt]\t" << result << endl;
			cout << "\t[solution]" << endl;
			x.t().print();
			return 0;
		}
		else {
			RecurBranch(root.lowBranch(index, x), opt, x_opt, "root->L");
			RecurBranch(root.upBranch(index, x), opt, x_opt, "root->R");
			cout << "\t[opt]\t" << opt << endl;
			cout << "\t[solution]" << endl;
			x_opt.t().print();
			getchar();
			return 0;
		}
		break;
	case Status::none: cout << "No feasible solutoin" << endl; return -1;
	case Status::unbounded: cout << "Unbounded solution" << endl; return -1;
	case Status::unique: 
		if (Branch::FindFirstNotInt(x) == -1) {
			cout << "The unique solution is the only Integer solution: ";
			x.t().print();
			return 0;
		}
		else {
			cout << "No integer solution, the unique solution contains non-integer:";
			x.t().print();
			return 0;
		}
		break;
	default:break;
	}
	getchar();
	return 0;
}

int RecurBranch(const Branch& bch, double& opt, vec& opt_x, string mark)
{
	Branch cpy = bch;
	double result = 0;
	vec x = vec();
	cout << "\t[Branch]\t" << mark << endl;
	cout << "\t[Matrix]" << endl;
	bch.Ab.print();
	Status s1 = SimplexNormalFormWithSlackVariables(cpy.f, cpy.Ab, cpy.base, cpy.arti, x, result);
	cout << "\t[status] " << s1 << endl;
	if (s1 == 0 || s1 == 1) {
		cout << "\t[result] " << result << endl;
		cout << "\t[solution]" << endl;
		x.t().print();
	}
	cout << endl;
		
	int index = 0;
	switch (s1) {
		
	case Status::infinite: 
		index = Branch::FindFirstNotInt(x);
		if (index == -1) {
			if (result > opt) {
				opt = result;
				opt_x = x;
			}
			return 0;
		}
		else {
			RecurBranch(cpy.lowBranch(index, x), opt, opt_x, mark+"->L");
			RecurBranch(cpy.upBranch(index, x), opt, opt_x, mark+"->R");
			return 0;
		}
		break;
	case Status::none: return -1;
	case Status::unbounded:  return -2;
	case Status::unique:
		if (Branch::FindFirstNotInt(x) == -1) {
			if (result > opt) {
				opt = result;
				opt_x = x;
			}
			return 0;
		}
		else return -3;
	default: break;
	}
	return 0;
}
