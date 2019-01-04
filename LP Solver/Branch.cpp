#include "Branch.h"
#include <sstream>
using namespace std;
int Branch::FindFirstNotInt(vec x) {
	for(int i=0 ; ;)
	return -1;
}

void Branch::parse(const int& cnt, const string& ff, const vector<string>& cst)
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
			this->f[index-1] = positive ? 1 : -1;
			part.clear();
		}
		else {
			part << seg.substr(1);
			part >> coef;
			part.clear();
			ss >> seg;
			part << seg.substr(seg.find("C") + 1);
			part >> index;
			this->f[index-1] = positive ? coef : -coef;
			part.clear();
		}
	}
	this->f.t().print();
	
	int mode = 0;
	int b = 0;
	for (int i = 0; i < cst.size(); i++) {
		stringstream cs;
		cs << cst[i];
		while (!cs.eof()) {
			cs >> seg;
			if (seg[0] == '<' || seg[0] == '>' || seg[0] == '=')
				break;
			positive = seg[0] == '+' ? true : false;
			if (seg.find("C") != string::npos) {
				part << seg.substr(seg.find("C") + 1);
				part >> index;
				this->Ab.at(i,index-1) = positive ? 1 : -1;
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
				this->Ab.at(i,index-1) = coef;
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
		this->mode.push_back(mode);
		
		cs >> b;
		this->Ab.at(i, cnt) = b;
		cs.clear();
	}
	//this->Ab.print();
}

void Branch::normal() {
	vector<arma::uword> tbase;
	vector<arma::uword> tart;
	int cnum = (int)this->Ab.n_rows;
	
	for (int i = 0; i < this->mode.size(); i++) {
		colvec cv = colvec(cnum).fill(0);
		colvec cv2 = colvec(cnum).fill(0);
		switch (this->mode[i]) {
		case 0: cv[i] = 1; this->Ab.insert_cols(this->Ab.n_cols - 1, cv); 
				tbase.push_back(this->Ab.n_cols-2); break;
		case 1: cv[i] = -1; this->Ab.insert_cols(this->Ab.n_cols - 1, cv); 
				cv2[i] = 1; this->Ab.insert_cols(this->Ab.n_cols - 1, cv2); 
				tbase.push_back(this->Ab.n_cols - 2); 
				tart.push_back(this->Ab.n_cols - 2); break;
		case 2: cv[i] = 1; this->Ab.insert_cols(this->Ab.n_cols - 1, cv); 
				tbase.push_back(this->Ab.n_cols - 2);
				tart.push_back(this->Ab.n_cols - 2); break; break;
		default:break;
		}
	}
	this->base = uvec(tbase);
	this->arti = uvec(tart);
	this->base.t().print();
	this->arti.t().print();
	this->Ab.print();
}