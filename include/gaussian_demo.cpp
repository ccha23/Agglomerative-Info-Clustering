#include <IC/AIC>
#include <IC/ChowLiu>
#include <IC/gaussian>

using namespace std;
using namespace Eigen;
using namespace IC;      

int main() {
	size_t k = 5;
	double sigma = 1;
	size_t n = 50;

	cout << "n = " << n << endl;

	MatrixXd A = MatrixXd::Zero(n, n + k);
	for (size_t i = 0; i < A.rows(); i++) {
		A(i, i % k) = 1;
		A(i, k + i) = sigma;
	}
	MatrixXd S = A* A.transpose();
    cout << "Covariance matrix S = A * A' :" << endl;
	cout << "S= \n" << S << endl;
    cout << "where A= \n" << A << endl;

	// generate the entropy function
	GaussianEntropy gsf(S);
	// generate approximate solution via CL tree
	cout << "Info-clustering by CL tree approximation:" << endl;
	vector<size_t> first_node, second_node;
	vector<double> gamma;
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < i; j++) {
			first_node.push_back(i);
			second_node.push_back(j);
			double I = gsf(vector<size_t> {i}) + gsf(vector<size_t> {j}) - gsf(vector<size_t> {i, j});
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
	AIC psp(gsf);
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
    
    cout << std::log(10) << " vs " << log(10) << endl;
}