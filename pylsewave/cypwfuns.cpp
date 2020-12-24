#include <math.h>
#include <iostream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include "include/cypwfuns.h"
#include <fstream>
#include <utility>

//template<typename T>
std::vector<double> funs::linspace(double start_in, double end_in, int num_in)
{

	std::vector<double> linspaced;

	double start = static_cast<double>(start_in);
	double end = static_cast<double>(end_in);
	double num = static_cast<double>(num_in);

	if (num == 0) { return linspaced; }
	if (num == 1)
	{
		linspaced.push_back(start);
		return linspaced;
	}

	double delta = (end - start) / (num - 1);

	for (int i = 0; i < num - 1; ++i)
	{
		linspaced.push_back(start + delta * i);
	}
	linspaced.push_back(end); // I want to ensure that start and end
							  // are exactly the same as the input
	return linspaced;
}

std::vector<double> funs::gradient(std::vector<double> vinput, double dx) {
	std::vector<double> vout;

	std::size_t vinput_size = vinput.size();
	if (vinput_size == 1)
	{
		vout = vinput;
	}
	if (vinput_size == 2)
	{
		vout.push_back((vinput[1] - vinput[0]) / dx); //first forward 
		vout.push_back((vinput[1] - vinput[0]) / dx); //first backward
	}
	if (vinput_size > 2)
	{
		vout.push_back((vinput[1] - vinput[0]) / dx); //first forward
		for (std::vector<int>::size_type i = 1; i != vinput.size() - 1; i++)
		{
			vout.push_back((vinput[i+1] - vinput[i-1]) / (2*dx)); //central difference
		}
		vout.push_back((vinput[vinput_size - 1] - vinput[vinput_size - 2]) / dx); //firts backward
	}

	return vout;
}

int funs::grad(double* inarr, double* out, double dx, size_t N_n){

	if (N_n == 1)
	{
		out[0] = inarr[0];
	}
	if (N_n == 2)
	{
		out[0] = (inarr[1] - inarr[0]) / dx; //first forward
		out[1] = (inarr[1] - inarr[0]) / dx; //first backward
	}
	if (N_n > 2)
	{
		out[0] = (inarr[1] - inarr[0]) / dx; //first forward
		for (size_t i = 1; i < N_n - 1; i++)
		{
			out[i] = (inarr[i+1] - inarr[i-1]) / (2*dx); //central difference
		}
		out[N_n - 1] =  (inarr[N_n - 1] - inarr[N_n - 2]) / dx; //firts backward
	}

	return 0;

}

int funs::tdma(double* a, double* b, double* c,
	double* d, double* out, std::size_t N_n) {

	//size_t N_n = d.size();
//
	std::vector<double> c_star(N_n, 0.0);
	std::vector<double> d_star(N_n, 0.0);

	c_star[0] = c[0] / b[0];
	d_star[0] = d[0] / b[0];

	for (std::size_t i = 1; i < N_n - 1; i++)
	{
		double m = 1.0 / (b[i] - a[i - 1] * c_star[i - 1]);
		c_star[i] = c[i] * m;
		d_star[i] = (d[i] - a[i - 1] * d_star[i - 1])*m;
	}
//
//   std::cout << N_n << std::endl;
//
	d_star[N_n - 1] = (d[N_n - 1] - a[N_n - 2] * d_star[N_n - 2]) / (b[N_n - 1] - a[N_n - 2] * c_star[N_n - 2]);

	out[N_n - 1] = d_star[N_n - 1];
	for (std::size_t i = N_n - 1; i-- > 0;)
	{
//	    std::cout << i << std::endl;
		out[i] = d_star[i] - c_star[i] * out[i + 1];
	}

	return 0;

}

double funs::std_dev(double *arr, std::size_t siz) {
	double mean = 0.0;
	double sum_sq;
	double *pVal;
	double diff;

	pVal = arr;
	for (size_t i = 0; i < siz; ++i, ++pVal) {
		mean += *pVal;
	}
	mean /= siz;

	pVal = arr;
	sum_sq = 0.0;
	for (std::size_t i = 0; i < siz; ++i, ++pVal) {
		diff = *pVal - mean;
		sum_sq += diff * diff;
	}
	return sqrt(sum_sq / siz);
}

int ios_efile::write2file(std::string filename, double *u, int rws, int clms) {
	std::ofstream ofp;
	//std::cout << "yes" << std::endl;
	ofp.open(filename.c_str());
	for (int i = 0; i < rws; i++) {
		for (int j = 0; j < clms; j++) {
			//ofp.write((char*)&u[i*rws + j], sizeof(double));
			ofp << std::setprecision(18) << std::scientific <<u[i*clms + j]<< " ";
		}
		ofp << "\n";
	}
	ofp.close();
	return 0;
}