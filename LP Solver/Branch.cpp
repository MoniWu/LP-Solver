#include "Branch.h"
#include <math.h>
#include <sstream>
using namespace std;
int Branch::FindFirstNotInt(vec x) {
	for (int i = 0; i < x.n_rows; i++) {
		int temp = (int)x[i];
		if ((x[i] - temp) > this->EPSILON)
			return i;
	}
	return -1;
}

void Branch::parse(const int& cnt, const string& ff, const vector<string>& cst, vec& f0, mat& Ab0, vector<int>& mode0)
{
	
	stringstream ss;
	stringstream part;
	string seg;
	int index = 0;
	int coef = 0;
	bool positive = false;
	
	ss << ff;
	while (!ss.eof()) {
		ss >> seg;
		positive = seg[0] == '+' ? true : false;
		if (seg.find("C") != string::npos) {
			part << seg.substr(seg.find("C")+1);
			part >> index;
			f0[index-1] = positive ? 1 : -1;
			part.clear();
		}
		else {
			part << seg.substr(1);
			part >> coef;
			part.clear();
			ss >> seg;
			part << seg.substr(seg.find("C") + 1);
			part >> index;
			f0[index-1] = positive ? coef : -coef;
			part.clear();
		}
	}
	
	int mode = 0;
	int b = 0;
	for (int i = 0; i < cst.size(); i++) {
		stringstream cs;
		cs << cst[i];
		while (!cs.eof()) {
			cs >> seg;
			if (seg[0] == 'R') continue;
			if (seg[0] == '<' || seg[0] == '>' || seg[0] == '=')
				break;
			positive = seg[0] == '-' ? false : true;
			if (seg.find("C") != string::npos) {
				part << seg.substr(seg.find("C") + 1);
				part >> index;
				Ab0.at(i,index-1) = positive ? 1 : -1;
				part.clear();
				seg.clear();
			}
			else {
				part << seg;
				part >> coef;
				part.clear();
				cs >> seg;
				part << seg.substr(seg.find("C") + 1);
				part >> index;
				Ab0.at(i,index-1) = coef;
				part.clear();
				seg.clear();
			}
		}

		switch (seg[0]) {
		case '<': mode = 0; break;
		case '>': mode = 1; break;
		case '=': mode = 2; break;
		default: break;
		}
		mode0.push_back(mode);
		
		cs >> b;
		Ab0.at(i, cnt) = b;
		cs.clear();
	}
	//this->Ab.print();
}

void Branch::normal(vec& f0, mat& Ab0, vector<int>& mode0, uvec& base0, uvec& arti0) {
	vector<arma::uword> tbase;
	vector<arma::uword> tart;
	int cnum = (int)Ab0.n_rows;
	
	for (int i = 0; i < mode0.size(); i++) {
		colvec cv = colvec(cnum).fill(0);
		colvec cv2 = colvec(cnum).fill(0);
		switch (mode0[i]) {
		case 0: cv[i] = 1; Ab0.insert_cols(Ab0.n_cols - 1, cv); 
				tbase.push_back(Ab0.n_cols-2); break;
		case 1: cv[i] = -1; Ab0.insert_cols(Ab0.n_cols - 1, cv); 
				cv2[i] = 1; Ab0.insert_cols(Ab0.n_cols - 1, cv2); 
				tbase.push_back(Ab0.n_cols - 2); 
				tart.push_back(Ab0.n_cols - 2); break;
		case 2: cv[i] = 1; Ab0.insert_cols(Ab0.n_cols - 1, cv); 
				tbase.push_back(Ab0.n_cols - 2);
				tart.push_back(Ab0.n_cols - 2); break;
		default:break;
		}
	}
	base0 = uvec(tbase);
	arti0 = uvec(tart);
	f0.insert_rows(f0.n_rows, vec(Ab0.n_cols - f0.n_rows).fill(0));
}

Branch Branch::upBranch(int index, vec x)
{
	vec new_f = this->f;
	mat new_Ab = this->Ab;
	uvec new_base = this->base;
	uvec new_arti = this->arti;

	rowvec rv = rowvec(new_Ab.n_cols).fill(0);
	rv[index] = 1;
	rv[new_Ab.n_cols - 1] = ((int)x[index]) + 1;

	for (int i = 0; i < new_Ab.n_rows; i++) {
		if (fabs(new_Ab.row(i)[index] - 0) > this->EPSILON) {
			rv -= new_Ab.row(i);
			break;
		}
	}
	cout << "new: ";
	rv.print();
	if (rv[rv.n_cols - 1] >= 0) {

		new_Ab.insert_rows(new_Ab.n_rows, rv);

		colvec cv = colvec(new_Ab.n_rows).fill(0);
		cv[new_Ab.n_rows - 1] = -1;
		colvec cv2 = colvec(new_Ab.n_rows).fill(0);
		cv2[new_Ab.n_rows - 1] = 1;
		new_Ab.insert_cols(new_Ab.n_cols - 1, cv);
		new_Ab.insert_cols(new_Ab.n_cols - 1, cv2);

		new_base.insert_rows(new_base.n_rows, uvec({ new_Ab.n_cols - 2 }));
		new_arti.insert_rows(new_arti.n_rows, uvec({ new_Ab.n_cols - 2 }));
		new_f.insert_rows(new_f.n_rows - 1, vec(2).fill(0));
	}
	else {
		rv = -rv;
		new_Ab.insert_rows(new_Ab.n_rows, rv);

		colvec cv = colvec(new_Ab.n_rows).fill(0);
		cv[new_Ab.n_rows - 1] = 1;
		new_Ab.insert_cols(new_Ab.n_cols - 1, cv);
		new_base.insert_rows(new_base.n_rows, uvec({ new_Ab.n_cols - 2 }));
		new_f.insert_rows(new_f.n_rows - 1, vec({ 0 }));
	}
	return Branch(new_f, new_Ab, new_base, new_arti);
}

Branch Branch::lowBranch(int index, vec x)
{
	vec new_f = this->f;
	mat new_Ab = this->Ab;
	uvec new_base = this->base;
	uvec new_arti = this->arti;

	rowvec rv = rowvec(new_Ab.n_cols).fill(0);
	rv[index] = 1;
	rv[new_Ab.n_cols - 1] = (int)x[index];

	for (int i = 0; i < new_Ab.n_rows; i++) {
		if (fabs(new_Ab.row(i)[index] - 0) > this->EPSILON) {
			rv -= new_Ab.row(i);
			break;
		}
	}
	cout << "new: ";
	rv.print();
	if (rv[rv.n_cols - 1] >= 0) {
		new_Ab.insert_rows(new_Ab.n_rows, rv);

		colvec cv = colvec(new_Ab.n_rows).fill(0);
		cv[new_Ab.n_rows - 1] = 1;
		new_Ab.insert_cols(new_Ab.n_cols - 1, cv);
		new_base.insert_rows(new_base.n_rows, uvec({ new_Ab.n_cols - 2 }));
		new_f.insert_rows(new_f.n_rows - 1, vec({ 0 }));
	}
	else {
		rv = -rv;
		new_Ab.insert_rows(new_Ab.n_rows, rv);

		colvec cv = colvec(new_Ab.n_rows).fill(0);
		cv[new_Ab.n_rows - 1] = -1;
		colvec cv2 = colvec(new_Ab.n_rows).fill(0);
		cv2[new_Ab.n_rows - 1] = 1;
		new_Ab.insert_cols(new_Ab.n_cols - 1, cv);
		new_Ab.insert_cols(new_Ab.n_cols - 1, cv2);

		new_base.insert_rows(new_base.n_rows, uvec({ new_Ab.n_cols - 2 }));
		new_arti.insert_rows(new_arti.n_rows, uvec({ new_Ab.n_cols - 2 }));
		new_f.insert_rows(new_f.n_rows - 1, vec(2).fill(0));
	}
	return Branch(new_f, new_Ab, new_base, new_arti);
}