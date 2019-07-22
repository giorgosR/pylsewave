#ifndef __CYPWFDM_H__
#define __CYPWFDM_H__

#include <vector>
#include <iostream>
#include <string>
#include "cypwmesh.h"

namespace shapes {
    class Circle{
    public:
        double xx;
        Circle(double x0, double y0, double r0);
        ~Circle();
        double getX();
        void setX(double x0);
        double getY();
        double getRadius();
        double getArea();
        void setCenter(double x0, double y0);
        void setRadius(double r0);
		double sum_mat(std::vector< std::vector<double> > & sv);
    private:
        double x;
        double y;
        double r;
    };
}

namespace numfuncs{
	void ccmultiply4d(double* array, double multiplier, int m, int n, int o, int p);
	std::vector<std::vector<double>> advance_solution(std::vector<std::vector<double>> const& u_n, int n);

	class FDM {
	public:
		FDM(VesselNetwork vesselNetwork);
		~FDM();
		void solve(char* casename);
		VesselNetwork getVesselNetwork();
		void setVesselNetwork(VesselNetwork);
	protected:
		VesselNetwork vesNet;
		double dx;
		double dt;
		double *x;
		double *t;
		void advance(double *u);

	};

}

#endif
