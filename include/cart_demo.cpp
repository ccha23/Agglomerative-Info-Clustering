#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <random>
#include <IC/AIC>
#include <time.h>
#include <fstream>
#include <string>

using namespace std;
using namespace Eigen;
using namespace IC;

int main(int argc, char **argv) {
    std::ifstream f;
	string model;
    if (argc > 2) {         /* if argument given */
        f.open (argv[1]);   /* open file with filename as argument */
        if (! f.is_open()) {    /* validate file open for reading */
            std::cerr << "error: file open failed '" << argv[1] << "'.\n";
            return 1;
        }
		model = argv[2];
    }
    else {  /* no argument given, error show usage */
        std::cerr << "error: insufficient input. <filename> <model(\"cart\" or \"svr\")required.\n";
        return 1;
    }

    std::string line, val;                  /* string for line & value */
    std::vector<std::vector<double>> array;    /* vector of vector<int>  */

    while (std::getline (f, line)) {        /* read each line */
        std::vector<double> v;                 /* row vector v */
        std::stringstream s (line);         /* stringstream line */
        while (getline (s, val, ','))       /* get each value (',' delimited) */
            v.push_back (std::stod (val));  /* add to row vector */
        array.push_back (v);                /* add row vector to array */
    }

    for (auto& row : array) {               /* iterate over rows */
        for (auto& val : row)               /* iterate over vals */
            std::cout << val << "  ";       /* output value      */
        std::cout << "\n";                  /* tidy up with '\n' */
    }
	fflush(stdout);
    CartEntropy csf(array, model);
    // HardCodeEntropy csf(0);
    // CartEntropy csf("../../data/GDS3893.soft");
    size_t n = csf.size();
	// generate approximate solution via CL tree
	cout << "Info-clustering by CL tree approximation:" << endl;
	vector<size_t> first_node, second_node;
	vector<double> gamma;
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < i; j++) {
			first_node.push_back(i);
			second_node.push_back(j);
			double I = csf(vector<size_t> {i}) + csf(vector<size_t> {j}) - csf(vector<size_t> {i, j});
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
	AIC psp(csf);
	{
		size_t i = 0;
		//cout << psp.getPartition(-1) << endl;
		while (psp.agglomerate(1E-8, 1E-10)) {
			// cout << "agglomerates to " << psp.getPartition(-1) << " at critical value " << psp.getCriticalValues().back() << endl;
		}
	}

	vector<double> psp_gamma = psp.getCriticalValues();
	cout << "critical values : " << psp_gamma << endl;
	for (double gamma : psp.getCriticalValues()) {
		cout << "partition at threshold " << gamma << ":" << psp.getPartition(gamma) << endl;
	}
}
