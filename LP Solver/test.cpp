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
	int cnt;
	string f;
	vector<string> cst;
	read(cnt,f,cst);
	Branch bh = Branch(cnt, f, cst);

	getchar();
	return 0;
}
