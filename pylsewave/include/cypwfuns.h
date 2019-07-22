#ifndef __CYPWFUNS_H__
#define __CYPWFUNS_H__

#include <vector>
#include <iostream>
#include <string>
//#include <stdlib.h>



namespace funs {
	//template<typename T>
	std::vector<double> linspace(double, double, int);
	std::vector<double> gradient(std::vector<double>, double);
	int grad(double*, double*, double, size_t);
	int tdma(double*, double*, double*,
		double*, double*, size_t);

    double std_dev(double *arr, size_t siz);
}

namespace ios_efile {
	int write2file(std::string, double *, int, int);
}

#endif