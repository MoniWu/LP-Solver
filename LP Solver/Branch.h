#pragma once
#include <armadillo>
using namespace std;
using namespace arma;
class Branch {
public:
	vec f;
	mat Ab;
	vector<int> mode;
	uvec base;
	uvec arti;
public:
	Branch(int cnt, string ff, vector<string> cst) {
		int m = (int)cst.size();
		this->f = vec(cnt);
		this->f.fill(0);
		this->Ab = mat(m,cnt+1);
		this->Ab.fill(0);
		parse(cnt, ff, cst);
		normal();
	}
	int FindFirstNotInt(vec x);
	void parse(const int& cnt, const string& ff, const vector<string>& cst);
	void normal();
};