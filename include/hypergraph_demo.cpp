#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <random>
#include <IC/AIC>
#include <time.h>
#include <fstream>

using namespace std;
using namespace Eigen;
using namespace IC;

int main() {
	size_t n = 6;
	size_t m = 4;

	MatrixXd M = MatrixXd::Zero(n, m);
	M << 1, 1, 0, 0,
		 1, 1, 0, 0,
		 1, 0, 0, 0,
		 0, 0, 1, 0,
		 0, 0, 1, 0,
		 0, 0, 0, 1;
	
	cout << "Incidence matrix of the (weighted) hypergraph:" << endl << "M=[" << M << "]" << endl;

	// generate the entropy function
	HypergraphEntropy hsf(M);

	// generate approximate solution via CL tree
	cout << "Info-clustering by CL tree approximation:" << endl;
	vector<size_t> first_node, second_node;
	vector<double> gamma;
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < i; j++) {
			first_node.push_back(i);
			second_node.push_back(j);
			double I = hsf(vector<size_t> {i}) + hsf(vector<size_t> {j}) - hsf(vector<size_t> {i, j});
			gamma.push_back(I);
		}
	}
	CL cl(n, first_node, second_node, gamma);

	vector<double> cl_gamma = cl.getCriticalValues();
	cout << "critical values : " << cl_gamma << endl;
	for (double gamma : cl.getCriticalValues()) {
		cout << "partition at threshold " << gamma << ":" << cl.getPartition(gamma) << endl;
	}

	// generate the exact info-clustering solution via AIC
	cout << "Agglomerative info-clustering:" << endl;
	AIC psp(hsf);
	{
		size_t i = 0;
		//cout << psp.getPartition(-1) << endl;
		while (psp.agglomerate(1E-8, 1E-10)) {
			//cout << "agglomerates to " << psp.getPartition(-1) << " at critical value " << psp.getCriticalValues().back() << endl;
		}
	}

	vector<double> psp_gamma = psp.getCriticalValues();
	cout << "critical values : " << psp_gamma << endl;
	for (double gamma : psp.getCriticalValues()) {
		cout << "partition at threshold " << gamma << ":" << psp.getPartition(gamma) << endl;
	}
}