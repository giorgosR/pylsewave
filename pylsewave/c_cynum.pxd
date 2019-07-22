# __author__ = "Georgios E. Ragkousis"
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.map cimport map

cdef extern from "<utility>":
    vector[double]&& move(vector[double]&&)

cdef extern from "include/cypwfdm.h" namespace "shapes":
    cdef cppclass Circle:
        double xx
        Circle(double, double, double)
        double getX()
        double getY()
        void setX(double)
        double getRadius()
        double getArea()
        void setCenter(double, double)
        void setRadius(double)
        double sum_mat(vector[vector[double]]&)

cdef extern from "include/cypwmesh.h":
    cdef cppclass Vessel:
        Vessel(string, double, double, double, double,
               map[string, double], int) except +
        string getName()
        double getL()
        double getRadius_prox()
        double getRadius_dist()
        double getWall_th()
        int getId()
        double getdx()
        map[string, double] getRLC()
        vector[double] get_k_vector()
        vector[double] get_x()
        vector[double] getR0()
        vector[double] get_f_R0()
        vector[double] get_df_dR0()
        vector[double] get_df_dx()
        vector[double] get_f_R0_ph()
        vector[double] get_df_dR0_ph()
        vector[double] get_df_dx_ph()
        vector[double] get_f_R0_mh()
        vector[double] get_df_dR0_mh()
        vector[double] get_df_dx_mh()
        void setdx(double)
        void setRLC(map[string, double])
        void set_k_vector(vector[double])
        vector[double] interpolate_R0(double)

    cdef cppclass VesselScaled(Vessel):
        VesselScaled(string, double, double, double, double,
               map[string, double], int, double) except +
        void set_k_vector(vector[double], double, double, double)

cdef extern from "include/cypwmesh.h":
    cdef cppclass VesselNetwork:
        VesselNetwork(vector[Vessel], double, double, double, double, double) except +
        VesselNetwork(Vessel, double, double, double, double, double) except +
        vector[Vessel] get_Vessels()
        double get_dx()
        void set_dx(double)
        double get_p0()
        void set_p0(double)
        double get_Re()
        void set_Re(double)
        double get_rho()
        void set_rho(double)
        double get_delta()
        void set_delta(double)
        # void set_boundary_layer_th(double, int)

    cdef cppclass VesselNetworkSc(VesselNetwork):
        VesselNetworkSc(vector[Vessel], double, double, double, double, double, double, double) except +
        VesselNetworkSc(Vessel, double, double, double, double, double, double, double) except +
        void set_boundary_layer_th(double, int)

cdef extern from "include/cypwfuns.h" namespace "funs":
    int tdma(double *, double *, double *,
             double *, double *, size_t)
    double std_dev(double *, size_t)
    vector[double] gradient(vector[double], double)
    int grad(double*, double*, double, size_t)
#    vector[double] linspace(double, double, double)

cdef extern from "include/cypwfuns.h" namespace "ios_efile":
   int write2file(string, double *, int, int)

# cdef extern from "include/cypwfdm.h" namespace "numfuncs":
#     void ccmultiply4d(double*, double, int, int, int, int)
#     vector[vector[double]] advance_solution(vector[vector[double]] &, int)
