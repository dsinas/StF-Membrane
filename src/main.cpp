/********************************************************************************
 *
 * Usage: A program that performs Monte Carlo computations of the model for the
 * 2D membrane composition and its out-of-plane deformation. The model considers
 * a coupling between local composition and local curvature of the membrane.
 *
 * to compile: g++ main.cpp fftsrc/fftw++.cc -lfftw3 -lm
 *
 * This program is made freely available with the understanding that every copy
 * of this file must include src files and headers and that it comes without any
 * WITHOUT ANY WARRANTY.
 *
 ********************************************************************************/
 
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <getopt.h>

#include "CLASS_MEM.h"

using namespace std;

/** function to print a help message for the program to run ****************** */
void print_usage() {
	cerr << "Purpose:\t Compute the structure factor of the membrane compositions coupled to membrane curvature" << endl;
	cerr << "Usage:\t\t MEMBRANE [options]" << endl;
	cerr << "Available options (default in []):" << endl;
	cerr << "\t\t-h\t print this help message" << endl;
	cerr << "\t\t-L arg\t set size of the lattice L - default: [100]" << endl;
	cerr << "\t\t-J arg\t set Ising coupling J - default: [0.5)]" << endl;
	cerr << "\t\t-g arg\t set curvature coupling g - default: [1.0)]" << endl;
}


/** MAIN FUNCTION *********************************************************** **/

int main( int argc, char* argv[] ) {

    int Lx = 100, Ly = Lx; // the default lattice size
    double temperature = 1.; // set constant: Ising coupling changes
    double mesh_size = 1.; // lattice size a scales to one
    double J_const = 0.5; // the default value for Ising coupling
    double kappa_const = 20.0; // bending rigidity
    double sigma_const = 1.0; // surface tension
    double gamma_const = 1.0; // the default value for curvature coupling

    double sD = 0.5; // initial composiotion of the membrane: 50-50
    double sU = 1.-sD;

    srand((unsigned) time(NULL));
    seed((unsigned long) (rand()*102011));
	
    char opt;
	while ((opt = getopt(argc, argv, "hL:J:g:")) != -1) {
		switch (opt) {
		case 'h':
			print_usage();
			return 0;
		case 'L':
			Lx = atoi(optarg);
			Ly = Lx;
			break;
		case 'J':
			J_const = atof(optarg);
			break;
		case 'g':
			gamma_const = atof(optarg);
			break;
		default:
			print_usage();
			return 1;
		}
	}
	
	cout << "\nSystem runs at" << endl;
	cout << "L = " << Lx << endl;
    cout << "J = " << J_const << endl;
    cout << "g = " << gamma_const << endl;
	
    MembraneClass m_system(Lx, Ly);
    m_system.set_temperature(temperature);
    m_system.set_mesh_size(mesh_size);
    m_system.set_J_const(J_const);
    m_system.set_kappa(kappa_const);
    m_system.set_sigma(sigma_const);
    m_system.set_gamma(gamma_const);
    m_system.set_spin_composition(sD,sU);

    m_system.initialization();
    m_system.set_pre_coeffs();
    m_system.set_hq_coeffs();

    m_system.get_data_files();

}
