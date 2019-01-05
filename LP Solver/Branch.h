#pragma once
#include <armadillo>
using namespace std;
using namespace arma;
class Branch {
public:
	const double EPSILON = 1e-7;
	vec f;
	mat Ab;
	//vector<int> mode;
	uvec base;
	uvec arti;
public:
	Branch() {

	}
	Branch(vec f, mat Ab, uvec base, uvec arti) {
		this->f = f;
		this->Ab = Ab;
		this->base = base;
		this->arti = arti;
	}
	static void parse(const int& cnt, const string& ff, const vector<string>& cst, vec& f0, mat& Ab0, vector<int>& mode0);
	static void normal(vec& f0, mat& Ab0, vector<int>& mode0, uvec& base0, uvec& arti0);
	static int FindFirstNotInt(vec x);
	Branch upBranch(int index, vec x);
	Branch lowBranch(int index, vec x);
};