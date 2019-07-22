#include "include/cypwfdm.h"
#include <math.h>

using namespace shapes;
using namespace numfuncs;

Circle::Circle(double x0, double y0, double r0){
    x = x0;
    y = y0;
    r = r0;
    xx = 0;
}

Circle::~Circle(){
    
}

double Circle::getX()
{
    return x;
}

void Circle::setX(double x0){
    x = x0;
}

double Circle::getY(){
    return y;
}

double Circle::getRadius(){
    return r;
}

double Circle::getArea(){
    return pow(r, 2);
}

void Circle::setCenter(double x0, double y0){
    x = x0;
    y = y0;
}

void Circle::setRadius(double r0){
    r = r0;
}

double Circle::sum_mat(std::vector< std::vector<double> > & sv)
{

	double tot = 0;

	int svrows = sv.size();
	int svcols = sv[0].size();
	std::cout << "vector length " << svrows << " , " << svcols << std::endl;

	for (int ii = 0; ii<svrows; ii++)
	{
		for (int jj = 0; jj<svcols; jj++)
		{
			std::cout << "element before " << sv[ii][jj] << std::endl;
			sv[ii][jj] = 3.0;
			std::cout << "element " << sv[ii][jj] << std::endl;
		}
	}
	return tot;

}


void numfuncs::ccmultiply4d(double* array, double multiplier, int m, int n, int o, int p) {

	int i, j, k, l;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < o; k++)
				for (l = 0; l < p; l++)
				{
					array[i*m*n*o + j*n*o + k*o + l] = array[i*m*n*o + j*n*o + k*o + l]*multiplier;
					//printf("Array[%d*%d*%d*%d + %d*%d*%d + %d*%d + %d] * Multiplier  = Array[%d] * %d = %.lf\n", i, m, n, o, j, n, o, k, o, l, i*m*n*o + j*n*o + k*o + l, multiplier, array[i*m*n*o + j*n*o + k*o + l]);
				}

}

std::vector<std::vector<double>> numfuncs::advance_solution(std::vector<std::vector<double>> const& u_n, int n) {
	int svrows = u_n.size();
	int svcols = u_n[0].size();

	std::vector< std::vector<double> > tot;
	tot.resize(svrows, std::vector<double>(svcols, -1));


	std::cout << "vector length " << svrows << " , " << svcols << std::endl;

	for (int ii = 0; ii<svrows; ii++)
	{
		for (int jj = 0; jj<svcols; jj++)
		{
			tot.at(ii).at(jj) = (2 * u_n.at(ii).at(jj));
		}
	}
	return tot;
}