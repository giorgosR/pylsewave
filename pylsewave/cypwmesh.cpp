#include <math.h>
#include <iostream>
#include <string>
#include <map>
#include "include/cypwfuns.h"
#include "include/cypwmesh.h"

const double pi = 3.14159265358979323846;

double pwround(double number)
{
	return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}

// ------------------ VESSEL ------------------------- //
Vessel::Vessel(std::string const & name_, double L_, double R_proximal, double R_distal,
		double Wall_thickness, std::map<std::string,
	double> Windkessels = std::map<std::string, double>(), int id=0)
{
	name = name_;
	L = L_;
	R_prox = R_proximal;
	R_dist = R_distal;
	W_th = Wall_thickness;
	Id = id;
	if (Windkessels.empty() == false)
	{
		RLC = Windkessels;
	}

}

Vessel::~Vessel()
{
}

int Vessel::getId(){
	return Id;
}

std::string Vessel::getName(){
	return name;
}

double Vessel::getL(){
	return L;
}

double Vessel::getdx() {
	return dx;
}

double Vessel::getRadius_prox(){
	return R_prox;
}

double Vessel::getRadius_dist(){
	return R_dist;
}

double Vessel::getWall_th()
{
	return W_th;
}

std::map<std::string, double> Vessel::getRLC() {
	return RLC;
}

std::vector<double> Vessel::get_k_vector() {
	return k;
}

std::vector<double> Vessel::get_x() {
	return x;
}

std::vector<double> Vessel::get_f_R0() {
	return f_r0;
}

std::vector<double> Vessel::get_df_dR0() {
	return df_dr0;
}

std::vector<double> Vessel::get_df_dx() {
	return df_dx;
}

std::vector<double> Vessel::get_f_R0_ph() {
	return f_r0_ph;
}

std::vector<double> Vessel::get_df_dR0_ph() {
	return df_dr0_ph;
}

std::vector<double> Vessel::get_df_dx_ph() {
	return df_dx_ph;
}

std::vector<double> Vessel::get_f_R0_mh() {
	return f_r0_mh;
}

std::vector<double> Vessel::get_df_dR0_mh() {
	return df_dr0_mh;
}

std::vector<double> Vessel::get_df_dx_mh() {
	return df_dx_mh;
}

std::vector<double> Vessel::getR0() {
	return R0;
}

//
void Vessel::setdx(double dx_input)
{
	//dx = dx_input;
	if (static_cast<int>(pwround(L / dx_input) + 1) == 1)
	{
		x.push_back(0.);
		x.push_back(L);
	}
	else
	{
		x = funs::linspace(0., L, static_cast<int>(pwround(L / dx_input) + 1));
	}
	dx = x[1] - x[0];
	// calculate R0(x)
	if (R0.empty() != true)
	{
		R0.clear();
	}
	this->calculate_R0();
	if (k.empty() != true)
	{
		f_r0 = f(R0, k);
		df_dr0 = dfdr(R0, k);
		df_dx = funs::gradient(R0, dx);
		f_r0_ph = f(interpolate_R0(0.5), k);
		df_dr0_ph = dfdr(interpolate_R0(0.5), k);
		df_dx_ph = funs::gradient(interpolate_R0(0.5), dx);
		f_r0_mh = f(interpolate_R0(-0.5), k);
		df_dr0_mh = dfdr(interpolate_R0(-0.5), k);
		df_dx_mh = funs::gradient(interpolate_R0(-0.5), dx);
	}
}

void Vessel::setRLC(std::map<std::string, double> dinput) {
	RLC = dinput;
}

void Vessel::set_k_vector(std::vector<double> k_input) {
	k = k_input;
}

void Vessel::calculate_R0()
{
	//int size = static_cast<int>(x.size());
	for (std::vector<int>::size_type i = 0; i != x.size(); i++)
	{
		R0.push_back(R_prox*exp(log(R_dist / R_prox)*(x[i] / L)));
	}
}

std::vector<double> Vessel::interpolate_R0(double value) {
	std::vector<double> vout(x.size(), 0.0);
	for (std::vector<int>::size_type i = 0; i != x.size(); i++)
	{
		vout[i] = R_prox*exp(log(R_dist / R_prox)*((x[i] + dx*value )/ L));
	}
	return vout;
}

std::vector<double> Vessel::f(std::vector<double> R0_input, std::vector<double> k_input) {
	double k1 = k_input[0];
	double k2 = k_input[1];
	double k3 = k_input[2];
	std::vector<double> vout;
	for (std::vector<int>::size_type i = 0; i != R0_input.size(); i++)
	{
		vout.push_back((4 / 3.) * (k2 * exp(k3 * R0_input[i]) + k1));
	}

	return vout;
}

std::vector<double> Vessel::dfdr(std::vector<double> R0_input, std::vector<double> k_input) {
	double k1 = k_input[0];
	double k2 = k_input[1];
	double k3 = k_input[2];
	std::vector<double> vout;
	for (std::vector<int>::size_type i = 0; i != R0_input.size(); i++)
	{
		vout.push_back((4 / 3.) * k2 * k3 * exp(k3 * R0_input[i]));
	}

	return vout;
}

// ------------------ VESSEL SCALED ------------------------- //
VesselScaled::VesselScaled(std::string const & name_, double L_, double R_proximal, double R_distal,
	double Wall_thickness, std::map<std::string,
	double> Windkessels = std::map<std::string, double>(), int id = 0, double rc_ = 1.0) : Vessel(name_, L_, R_proximal,
		R_distal, Wall_thickness, Windkessels, id)
{
	rc = rc_;
	L = L / rc;
	R_prox = R_prox / rc;
	R_dist = R_dist / rc;
	W_th = W_th / rc;
	//std::cout << L << std::endl;
}

void VesselScaled::set_k_vector(std::vector<double> k_input, double rho=0, double qc=0, double rc=0) {
	k = k_input;
	if (rho != 0 && qc != 0 && rc != 0)
	{
		double kc;
		kc = (rho*pow(qc, 2) / pow(rc, 4));
		k[0] = k[0] / kc;
		k[1] = k[1] / kc;
		k[2] = k[2] * rc;
	}
}

VesselScaled::~VesselScaled()
{
}

// ------------------ VESSEL NETWORK ------------------------- //
VesselNetwork::VesselNetwork(std::vector<Vessel> invessels, double inp0, double inrho=-1, double inRe=-1,
	double indx=-1, double inNx=-1) : vessels(invessels), p0(inp0), rho(inrho), Re(inRe), dx(indx), Nx(inNx)
{
	if (dx != -1) {
		for (std::vector<int>::size_type i = 0; i != vessels.size(); i++) {
			vessels[i].setdx(dx);
		}
	}

}

VesselNetwork::VesselNetwork(Vessel invessel, double inp0, double inrho = -1, double inRe = -1,
	double indx = -1, double inNx = -1) : p0(inp0), rho(inrho), Re(inRe), dx(indx), Nx(inNx)
{
	vessels.push_back(invessel);
	if (dx != -1) {
		vessels[0].setdx(dx);
	}
}

//void VesselNetwork::set_boundary_layer_th(double T, int no_of_cycles)
//{
//	double tt = T;
//}

std::vector<Vessel> VesselNetwork::get_Vessels() {
	return vessels;
}

double VesselNetwork::get_dx() {
	return dx;
}

void VesselNetwork::set_dx(double indx)
{
	dx = indx;
}

double VesselNetwork::get_p0(){
	return p0;
}

void VesselNetwork::set_p0(double inp0){
	p0 = inp0;
}

double VesselNetwork::get_Re(){
	return Re;
}

void VesselNetwork::set_Re(double inRe){
	Re = inRe;
}

double VesselNetwork::get_delta() {
	return delta;
}

void VesselNetwork::set_delta(double indelta) {
	Re = indelta;
}

double VesselNetwork::get_rho(){
	return rho;
}

void VesselNetwork::set_rho(double inrho){
	rho = inrho;
}

VesselNetwork::~VesselNetwork()
{
}

// ------------------ VESSEL NETWORK SCALED ------------------------- //
VesselNetworkSc::VesselNetworkSc(std::vector<Vessel> invessels, double inp0, double inrho = -1, double inRe = -1,
	double indx = -1, double inNx = -1, double inqc = 1., double inrc=1.) : VesselNetwork(invessels, inp0, inrho, inRe,
		indx, inNx), qc(inqc), rc(inrc)
{
	if (dx != -1)
	{
		dx = dx / rc;
		for (std::vector<int>::size_type i = 0; i != vessels.size(); i++) {
			vessels[i].setdx(dx);
		}

	}
	Re = (qc) / (4.6*rc);
	if (p0 != -1)
	{
		p0 = (p0*pow(rc, 4)) / (rho*pow(qc, 2));
	}

}

VesselNetworkSc::VesselNetworkSc(Vessel invessel, double inp0, double inrho = -1, double inRe = -1,
	double indx = -1, double inNx = -1, double inqc = 1., double inrc = 1.) : VesselNetwork(invessel, inp0, inrho, inRe,
		indx, inNx), qc(inqc), rc(inrc)
{
	vessels.push_back(invessel);
	if (dx != -1) {
		dx = dx / rc;
		vessels[0].setdx(dx);
	}
	Re = (qc) / (4.6*rc);
	if (p0 != -1)
	{
		p0 = (p0*pow(rc, 4)) / (rho*pow(qc, 2));
	}
}

void VesselNetworkSc::set_boundary_layer_th(double inT, int in_no_of_cycles)
{
	double _nu, T_cycle;
	_nu = (4.6 * rc) / (qc);
	T_cycle = inT / in_no_of_cycles;
	delta = pow((_nu * T_cycle) / (2.* pi), 0.5);
}

VesselNetworkSc::~VesselNetworkSc()
{
}