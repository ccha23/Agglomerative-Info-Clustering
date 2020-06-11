#include <IC/AIC>
#include <IC/ChowLiu>
#include <IC/FLS>

using namespace std;
using namespace Eigen;
using namespace IC;

int main(){
	int pField = 2;
	int kField = 1;

	size_t m = 4;
	size_t n0 = 1, n1 = 1, n2 = 1, n3 = 2;
	MatrixXd M0(n0,m);
	MatrixXd M1(n1,m);
	MatrixXd M2(n2,m);
	MatrixXd M3(n3,m);
	M0 << 1,0,0,0;
	M1 << 0,1,0,0;  
	M2 << 1,1,0,0,
	M3 << 1,1,1,0,
          1,1,1,1;
// 	cout << "M0 = \n" << M0 << endl;
// 	cout << "M1 = \n" << M1 << endl;
// 	cout << "M2 = \n" << M2 << endl;
// 	cout << "M3 = \n" << M3 << endl;

	vector<MatrixXd> V;
	V = {M0,M1,M2,M3};


  	pari_init(1000000,2);

	// generate the entropy function
	FiniteLinearEntropy fsf(V,pField,kField);
	size_t n = V.size();

	// generate approximate solution via CL tree
	cout << "Info-clustering by CL tree approximation:" << endl;
	vector<size_t> first_node, second_node;
	vector<double> gamma;
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < i; j++) {
			first_node.push_back(i);
			second_node.push_back(j);
			double I = fsf(vector<size_t> {i}) + fsf(vector<size_t> {j}) - fsf(vector<size_t> {i, j});
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
	AIC psp(fsf);
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


   pari_close();

}

 
