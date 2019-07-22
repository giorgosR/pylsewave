__author__ = 'Georgios E. Ragkousis'
import numpy as np
import os
import matplotlib.pyplot as plt
from pylsewave.mesh import Vessel, VesselNetwork
from pylsewave.pdes import PDEm, PDEsWat
from pylsewave.bcs import BCs, BCsWat
from pylsewave.viz import PlotAndStoreSolution
from pylsewave.interpolate import CubicSpline
from pylsewave.nonlinearsolvers import Newton_system_conj_points
from pylsewave.fdm import BloodWaveMacCormack, BloodWaveLaxWendroff
from pylsewave.pwconsts import *
from pylsewave.pwutils import convert_data_periodic, compute_c, linear_extrapolation
# ---- SOLUTION WITH CYTHON CLASSES ----- #
from pylsewave.cynum import cPDEsWat, BCsADAN56, cMacCormackSolver
import argparse as argprs


class partVessel(Vessel):
    """
    A class that inherits from Vessel.
    It has been created for the AAA case.
    """
    @property
    def r0(self):
        return self._R0

    @r0.setter
    def r0(self, value):
        self._R0 = value

    @property
    def id(self):
        return self._Id

    @id.setter
    def id(self, value):
        self._Id = value

    def calculate_R0(self, x):
        return self._R0

    def interpolate_R0(self, value):
        return self._R0


class Aneurysm(Vessel):
    """
    Class inheritted from base Vessel.
    With this class, a AAA geometry can be defined.
    """
    def __init__(self, *args, **kwargs):
        super(Aneurysm, self).__init__(*args, **kwargs)
        self.sigma = None
        self.beta = None
        self.max_radius = None
        self.centre = None
        self.r_anter = None
        self.r_post = None

    def set_params(self, sigma=30.0, beta=1.0, aneurysm_factor=0.5):
        """Class method for Gaussian like bell definition
        Args:
            sigma (float): sigma factor of Gaussian distribution
            beta (float): beta parameter
            aneurysm_factor (float): aneurysm scaling factor

        Returns:
            None

        """
        self.sigma = sigma
        self.beta = beta
        self.max_radius = self.r_prox * (aneurysm_factor + 1.)

    def calculate_R0(self, x):
        """
        Class method which overides the respective calculate_R0 of Vessel
        :param x: spatial data
        :return: Reference radius size of the AAA
        """
        z = np.zeros(self.x.shape[0])
        for i in range(self.x.shape[0]):
            z[i] = -0.5 * self.length + i * self.length / (self.x.shape[0] - 1)
        # z = 0.5*self.length - x*self.length/(self.x.shape[0] - 1)
        #         print z
        point_gauss = (1 / ((2. * np.pi * self.sigma * self.sigma) ** 0.5)) * (
        np.exp(1) ** (-0.5 * z ** 2 / self.sigma ** 2))
        gauss_peak = np.max(point_gauss)
        gauss_min = np.min(point_gauss)
        scale = 2. * (self.max_radius - self.r_prox) / (gauss_peak + self.beta * gauss_peak)
        # CALCULATION OF POINTS for ANTERIOR PROFILE
        self.r_anter = (point_gauss * scale + self.r_prox)
        #         self.r_anter = ((point_gauss - gauss_min)*(self.max_radius - self.r_prox))/((gauss_peak - gauss_min) + self.r_prox)
        #         print self.r_anter
        # CALCULATION OF POINTS for POSTERIOR PROFILE
        self.r_post = -(point_gauss * scale * self.beta + self.r_prox)
        #         self.r_post = -((point_gauss - gauss_min)*(self.max_radius - self.r_prox))/((gauss_peak - gauss_min) + self.r_prox)*self.beta
        # CALCULATION OF POINTS for CENTRAL LINE
        self.centre = 0.5 * (self.r_anter + self.r_post)
        # Z_TRANSLATION
        #         z = z + 0.5*self.length

        return 0.5 * (self.r_anter - self.r_post)

    def splitSegments(self, no_of_segments, list_borders, uniform=False):
        """
        Class method in case the user wants to split the AAA geometry in different vessels.
        :param no_of_segments (int): the number of segments that the geometry will be split
        :param list_borders (list [int]): a list with indices indication where the split will be performed.
        :param uniform (bool): true if the vessel will be split uniformingly False otherwise
        :return (list [Vessel]): A list containing all the vessels comprising the AAA geometry.
        """
        segments = []
        pos = 0
        previous_pos = 0
        if uniform == False:
            for i in range(no_of_segments):
                len_segment = list_borders[i][1] - list_borders[i][0]
                previous_pos += pos
                r_prox = self.r0[previous_pos]
                pos = int(round(len_segment / self.dx))
                r_dist = self.r0[pos]
                wall_thick = self._Wall_th

                if i == 0:
                    r0s = np.zeros(pos + 1, np.float)
                    for j in range(pos + 1):
                        r0s[j] = self.r0[previous_pos + j]
                elif i != 0:
                    r0s = np.zeros(pos, np.float)
                    for j in range(pos):
                        r0s[j] = self.r0[previous_pos + j]

                segments.append(partVessel(name=self.name + "_PART_%d" % (i),
                                           L=len_segment, R_proximal=r_prox,
                                           R_distal=r_dist,
                                           Wall_thickness=wall_thick, Id=None))
                segments[i].r0 = r0s
                segments[i].set_k_vector(self._k)
                segments[i].dx = self.dx

                if i != 0:
                    segments[i].r0 = np.insert(segments[i].r0, 0, segments[i - 1].r0[-1])
                    segments[i].dx = self.dx

        return segments


class ADANBC(BCsWat):
    """
    Class which inherits from BCsWat and it is used to prescribe inlet BCs.
    """
    def U_0(self, u, t, dx, dt, vessel_index, out):
        theta = dt / dx
        dt2 = 0.5 * dt

        q_pres = self.inlet_fun(t)
        q_1 = out[1, 1]

        A = u[0, 0] - (2 * theta) * (q_1 - q_pres)

        out[0, 0] = A
        out[1, 0] = q_pres

def run_Adan_77_case(idatafile, ibcsinflowfile, oresfile, language,
                     verbose=True):
    print(type(language))
    filename = idatafile
    data = np.loadtxt(filename, delimiter="&", dtype=np.str)
    if verbose is True:
        print(" \\\\\n".join([" & ".join(map(str, line)) for line in data]))

    PRINT_STATUS = False
    WRITE_STATUS = True

    # Mynard
    k = np.array([33.7e-03, 0.3, -0.9])

    nu = CONSTANT_mu / CONSTANT_rho

    T_cycle = 1.0
    tc = 4
    T = T_cycle * tc

    dt = 1e-4

    p0 = 0.01  # Mpa

    # -------  LOAD ARTERIAL SEGMENTS  ------- #
    print k
    segments = []
    for i in range(data.shape[0]):
        segments.append(Vessel(name=data[i, 1], L=float(data[i, 2]) * 10.,
                               R_proximal=float(data[i, 3]) * 10.,
                               R_distal=float(data[i, 4]) * 10.,
                               Wall_thickness=float(data[i, 5]) * 10., Id=i))
        # set k vector
        segments[i].set_k_vector(k=k)

    # -------  INFLOW (IN VIVO) WAVE  ------- #
    invivo_data = np.loadtxt(ibcsinflowfile, delimiter=" ")
    time_measured = invivo_data[:, 0]
    flow_measured = invivo_data[:, 1] * 1000.

    time_periodic, flow_periodic = convert_data_periodic(time_measured,
                                                         flow_measured,
                                                         cycles=tc, plot=True)
    q_inlet_bc = CubicSpline(time_periodic, flow_periodic)

    # -------  TERMINAL VESSELS  ------- #
    terminal_vessels = {5: [18104., 72417., 3.129e-06], 9: [11539., 46155., 4.909e-06],
                        12: [47813., 191252., 1.185e-06], 13: [11749., 46995., 4.821e-06],
                        14: [9391., 37563., 6.032e-06], 15: [5760., 23041., 9.833e-06],
                        18: [9424., 37696., 6.011e-06], 19: [5779., 23118., 9.801e-06],
                        23: [19243., 76972., 2.944e-06], 27: [11332., 45329., 4.998e-06],
                        30: [47986., 191945., 1.180e-06], 31: [11976., 47905., 4.730e-06],
                        32: [249127., 996508., 2.274e-07], 34: [255583., 1022333., 2.216e-07],
                        36: [232434., 929735., 2.437e-07], 38: [234425., 937702., 2.416e-07],
                        43: [3349., 13394., 1.692e-05], 45: [343394., 1373574., 1.650e-07],
                        46: [4733., 18933., 1.197e-05], 47: [2182., 8728., 2.596e-05],
                        49: [2263., 9051., 2.503e-05], 51: [2270., 9082., 2.495e-05],
                        53: [23913., 95652., 2.369e-06], 59: [4146., 16582., 1.366e-05],
                        60: [3427., 13707., 1.653e-05], 63: [24525., 98100., 2.310e-06],
                        66: [21156., 84625., 2.677e-06], 69: [4158., 16632., 1.362e-05],
                        70: [3429., 13715., 1.652e-05], 73: [24533., 98131., 2.309e-06],
                        76: [21166., 84662., 2.676e-06]}

    for i in terminal_vessels.keys():
        terminal_vessels[i][0] = terminal_vessels[i][0] * 1e-010
        terminal_vessels[i][1] = terminal_vessels[i][1] * 1e-010
        terminal_vessels[i][2] = terminal_vessels[i][2] * 1e+010

    # -------  BIFURCATIONS  ------- #
    bif_vessels = [[0, 1, 2], [1, 3, 4], [3, 5, 6], [4, 14, 15],
                   [8, 9, 10], [10, 11, 13], [2, 16, 17], [16, 18, 19],
                   [17, 20, 21], [20, 23, 24], [26, 27, 28], [28, 29, 31],
                   [22, 32, 33], [33, 34, 35], [35, 36, 37], [37, 38, 39],
                   [40, 41, 42], [41, 43, 44], [44, 45, 46], [42, 47, 48],
                   [48, 49, 50], [50, 51, 52], [52, 53, 54], [54, 55, 56],
                   [55, 57, 59], [58, 60, 61], [62, 63, 64], [56, 67, 69],
                   [68, 70, 71], [72, 73, 74]]

    # -------  CONJUCTIONS  ------- #
    conj_points = [[6, 7], [7, 8], [11, 12], [24, 25], [25, 26], [29, 30],
                   [21, 22], [39, 40], [57, 58], [61, 62], [64, 65], [65, 66],
                   [67, 68], [71, 72], [74, 75], [75, 76]]

    for i in terminal_vessels.keys():
        c0_distal = compute_c(segments[i].r_dist, k)
        #     print c0_distal
        A0_distal = np.pi * ((segments[i].r_dist) ** 2)
        # R1 should be the same with the input characteristic impedance
        Z1_distal = (CONSTANT_rho * c0_distal) / A0_distal

        R1 = terminal_vessels[i][0]
        R2 = terminal_vessels[i][1]
        C_t = terminal_vessels[i][2]
        # add RLC data in each terminal vessel
        segments[i].RLC = {"R_1": Z1_distal, "R_t": R2, "C_t": C_t}

    # create the Arterial Network domain/mesh
    # Reflecting BCs
    Nx = None
    vesssel_network = VesselNetwork(vessels=segments,
                                    rho=CONSTANT_rho, Re=0.,
                                    p0=p0, dx=4.5, Nx=Nx)

    # give a name for the output database file
    casename = "/results/" + oresfile

    # check CFL and set dx accordingly
    # check CFL and set dx accordingly
    siz_ves = len(vesssel_network.vessels)
    compare_l_c0 = []
    for i in range(siz_ves):
        c_max = np.max(compute_c(vesssel_network.vessels[i].r0, k))
        A = np.pi*(vesssel_network.vessels[i].r_prox*vesssel_network.vessels[i].r_prox)
        compare_l_c0.append(vesssel_network.vessels[i].length / c_max)
    
    min_value = min(compare_l_c0)
    index_min_value = np.argmin(compare_l_c0)
    print("The min length to wave speed radio has been computed to Vessel: '%s' " % vesssel_network.vessels[index_min_value].name)
        
    # Nx_i = 1
    min_time = []
    for i in range(siz_ves):
        Nx_i = 4*np.floor((vesssel_network.vessels[i].length / compute_c(vesssel_network.vessels[i].r_prox, k))/(min_value))
        dx_i = vesssel_network.vessels[i].length / Nx_i
        vesssel_network.vessels[i].dx = dx_i
        min_time.append(dx_i/np.max(compute_c(vesssel_network.vessels[i].r0, k)))
    
    CFL = 0.5
    dt = CFL*(min(min_time))
    print(dt)

    # callback function to store solution
    number_of_frames = 200
    skip = int(round(T / dt)) / number_of_frames
    umin = 0.1
    umax = 1.5
    myCallback = PlotAndStoreSolution(casename=casename, umin=umin,
                                      umax=umax, skip_frame=skip,
                                      screen_movie=True, backend=None,
                                      filename='/results/tmpdata')
    if language == 'cy':
        myPDEs = cPDEsWat(vesssel_network)
        myBCs = BCsADAN56(myPDEs, q_inlet_bc.eval_spline)
    
        U0_vessel = np.array([0],dtype=np.int)
        UL_vessel = np.array(terminal_vessels.keys())
        UBif_vessel = np.array(bif_vessels)
        UConj_vessel = np.array(conj_points)
    
        mySolver = cMacCormackSolver(myBCs)
        mySolver.set_T(dt=dt, T=T, no_cycles=tc)
        mySolver.set_BC(U0_vessel, UL_vessel, UBif_vessel, UConj_vessel)
        mySolver.solve(casename, myCallback)
        myCallback.close_file(casename)

    elif language == 'py':
        # PDEs #
        myPDEs = PDEsWat(vesssel_network)
        # BCS #
        myBCs = ADANBC(myPDEs, q_inlet_bc.eval_spline)
        U0_vessel = np.array([0],dtype=np.int)
        UL_vessel = np.array(terminal_vessels.keys())
        UBif_vessel = np.array(bif_vessels)
        UConj_vessel = np.array(conj_points)
        
        ### ----- PYTHON ------ ###
        mySolver = BloodWaveMacCormack(myBCs)
        mySolver.set_T(dt=dt, T=T, no_cycles=tc)
        mySolver.set_BC(U0_vessel, UL_vessel, UBif_vessel, UConj_vessel)
        ### ------- PYTHON ONLY ------------------ ###
        ### ------- SOLVE AND TIME --------------- ###
        mySolver.solve(casename, myCallback)
        myCallback.close_file(casename)

def main(*args):
    run_Adan_77_case(*args)

    return STATUS_OK

if __name__ == "__main__":
    import sys
    parser = argprs.ArgumentParser(description="Case study for a detailed arterial network" +
                                               "consisted of 77 arterial segments")
    parser.add_argument('-ivesseldatafile', help='input data file (.dat/.txt) with the vessel information')
    parser.add_argument('-ibcinflowfile', help='input data file (.dat/.txt) with the inflow information')
    parser.add_argument('-oresfile', help='output result zipped file name (.npz)')
    parser.add_argument('-language', type=str, help='Solve with Python or C (py or cy)')
    parser.add_argument('-ovtkfile', help='Output visualisation file (vtk/Paraview)')
    args = parser.parse_args()
    print(args)
    if (args.ivesseldatafile is None) or (args.ibcinflowfile is None) or (args.oresfile is None):
        print('File should be exected as:\n' + sys.argv[0] +
              " -ivesseldatafile <ivesseldatafile> -ibcinflowfile <ibcinflowfile>" +
              " -oresfile <oresfile>")
        print(STATUS_ERROR)
    else:
        status = main(args.ivesseldatafile, args.ibcinflowfile, args.oresfile,
                      args.language)
        print(status)
