#ifndef __CYPWMESH_H__
#define __CYPWMESH_H__

#include <vector>
#include <iostream>
#include <string>
#include <map>

class Vessel
{
public:
	Vessel(std::string const & name_, double L_, double R_proximal, double R_distal,
		double Wall_thickness, std::map<std::string, double> Windkessels, int id);
	virtual ~Vessel();
	//properties
	std::string getName();
	double getL();
	double getRadius_prox();
	double getRadius_dist();
	double getWall_th();
	double getdx();
	int getId();
	std::vector<double> get_x();
	std::map<std::string, double> getRLC();
	std::vector<double> get_k_vector();
	std::vector<double> getR0();
	std::vector<double> get_f_R0();
	std::vector<double> get_df_dR0();
	std::vector<double> get_df_dx();
	std::vector<double> get_f_R0_ph();
	std::vector<double> get_df_dR0_ph();
	std::vector<double> get_df_dx_ph();
	std::vector<double> get_f_R0_mh();
	std::vector<double> get_df_dR0_mh();
	std::vector<double> get_df_dx_mh();
	//members
	void setdx(double);
	void setRLC(std::map<std::string, double>);
	virtual void set_k_vector(std::vector<double>);
	std::vector<double> interpolate_R0(double value);
protected:
	std::string name;
	double L;
	double R_prox;
	double R_dist;
	double W_th;
	double dx;
	std::vector<double> R0;
	int Id;
	std::vector<double> x;
	std::map<std::string, double> RLC;
	std::vector<double> f_r0;
	std::vector<double> df_dr0;
	std::vector<double> df_dx;
	std::vector<double> f_r0_ph;
	std::vector<double> df_dr0_ph;
	std::vector<double> df_dx_ph;
	std::vector<double> f_r0_mh;
	std::vector<double> df_dr0_mh;
	std::vector<double> df_dx_mh;
	std::vector<double> k;
	void calculate_R0();
	static std::vector<double> f(std::vector<double>, std::vector<double>);
	static std::vector<double> dfdr(std::vector<double>, std::vector<double>);
private:

	
};

class VesselScaled : public Vessel
{
public:
	VesselScaled(std::string const & name_, double L_, double R_proximal, double R_distal,
		double Wall_thickness, std::map<std::string, double> Windkessels, int id, double rc_);
	virtual ~VesselScaled();
	virtual void set_k_vector(std::vector<double>, double, double, double);
private:
	double rc;
};


// ---- VESSEL NETWORK ---------- //
class VesselNetwork
{
public:
	VesselNetwork(std::vector<Vessel>, double, double, double, double, double);
	VesselNetwork(Vessel, double, double, double, double, double);
	~VesselNetwork();
	std::vector<Vessel> get_Vessels();
	double get_dx();
	void set_dx(double);
	double get_p0();
	void set_p0(double);
	double get_Re();
	void set_Re(double);
	double get_rho();
	void set_rho(double);
	double get_delta();
	void set_delta(double);
	//void set_boundary_layer_th(double, int);
private:

protected:
	std::vector<Vessel> vessels;
	double p0;
	double rho;
	double dx;
	double Re;
	int no_cycles;
	double delta;
	double Nx;

};


class VesselNetworkSc : public VesselNetwork
{
public:
	VesselNetworkSc(std::vector<Vessel>, double, double, double, double, double, double, double);
	VesselNetworkSc(Vessel, double, double, double, double, double, double, double);
	~VesselNetworkSc();
	void set_boundary_layer_th(double, int);
private:

protected:
	double rc;
	double qc;
};

#endif