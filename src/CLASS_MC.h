/*******************************************************************************
* 
* CLASS_MC.h: A class for implementing Monte Carlo move on 2D membrane and
* update height profile in the Fourier space stochastically.
*
* This program is made freely available with the understanding that every 
* copy of this file must include this header and that it comes without any 
* WITHOUT ANY WARRANTY.
*******************************************************************************/

#ifndef CLASS_MC_H
#define CLASS_MC_H

#include "RndGen.h"
#include "CLASS_GRID.h"
#include "fftsrc/Array.h"
#include "fftsrc/fftw++.h"

using namespace std;
using namespace Array;
using namespace fftwpp;

const double Boltzmann_const = 1.;
const int SpinD = -1;
const int SpinU = +1;

vector<vector<double> > sr_matrix;
vector<vector<double> > hr_matrix;
vector<vector<double> > sq_norm;
vector<vector<double> > hq_norm;
vector<vector<double> > sq_hq;
vector<vector<complex<double> > > hq_complex;
vector<vector<complex<double> > > sq_complex;

vector<vector<double> > preK;
vector<vector<double> > preS;
vector<vector<double> > pre_mu;
vector<vector<double> > sig_hq;

class MC_SystemClass {

protected:
    LatticeClass m_grid;
    int m_lx, m_ly;
    int nyp, nxp;
    double m_beta, m_temperature;
    double m_energy;
    int m_orderParameter;
    double mesh_size, J_const, KAPPA, SIGMA, GAMMA;

    double coupling_E, delta_coupling_E;
    double mem_E, delta_mem_E;

    double spinD_frac;
    double spinU_frac;
    int initial_OP;

    double mc_criterion;
    unsigned long int mc_counter;

public:
    MC_SystemClass(const int _Lx, const int _Ly):
            m_grid(_Lx,_Ly), m_lx(_Lx), m_ly(_Ly)
            , mesh_size(0.), J_const(0.), KAPPA(0.), SIGMA(0.), GAMMA(0.)
            , mc_counter(0)

    {
        nyp = m_ly/2+1;
        nxp = (m_lx+1)/2;
    }

    ~MC_SystemClass() {}

    void initialization() {

        initialize_lattices();

        m_orderParameter = get_order_parameter();
        m_energy = get_energy();

        int nSweep = 10; // on average every "10" sweeps height profile updated
        mc_criterion = 1./(double)(nSweep*m_grid.get_volume()); // this is given experimentally

    }

    /** initializing lattices ********************************************* **/
    void initialize_lattices() {

        /** initializing the lattice states **/
        for (int i=0;i<m_grid.get_volume();i++) {
            int spin;
            if (drandom()<spinD_frac)
                spin = SpinD;
            else
                spin = SpinU;

            m_grid.set_cellSpin(i,spin);
        }

        /** initializing the lattice heights **/
        for (int i=0;i<m_grid.get_volume();i++) {
            m_grid.set_cellHeight(i,0.*drandom());
        }

        set_order_parameter(); // it is not required for the grand canonical simulation.

        hq_complex.resize(m_lx,vector<complex<double> > (nyp));
        sq_norm.resize(m_lx,vector<double > (m_ly));
        hq_norm.resize(m_lx,vector<double > (m_ly));
        sq_hq.resize(m_lx,vector<double > (m_ly));

    }

    /** set order parameter corresponding to initial compositions ********* **/
    void set_order_parameter() {

        while (get_order_parameter()<initial_OP) {
            int idx = (int) (drandom()*(double)m_grid.get_volume());
            int cell_state = m_grid.get_cellSpin(idx);
            if (cell_state<0)
                m_grid.set_cellSpin(idx,SpinU);
        }

        while (get_order_parameter()>initial_OP) {
            int idx = (int) (drandom()*(double)m_grid.get_volume());
            int cell_state = m_grid.get_cellSpin(idx);
            if (cell_state>0)
                m_grid.set_cellSpin(idx,SpinD);
        }

    }

    /** ******************************************************************* **/
    int get_order_parameter() {
        int result(0);
        for (int i=0;i<m_grid.get_volume();i++)
            result += m_grid.get_cellSpin(i);
        return result;
    }

    /** ******************************************************************* **/
    /** calculating the total energy of the system in Fourier space ******* **/
    double get_energy_F() {

        double E0(0.), E1(0.), E2(0.), E3(0.);
        int spin, s_sum;

        for (int idx=0;idx<m_grid.get_volume();idx++) {

            /** energy contribution of spin-spin interaction **/
            spin = m_grid.get_cellSpin(idx);
            s_sum = m_grid.get_sum_neighbors_spin(idx);
            E0 += (double)(spin*s_sum);
        }

        update_sq_hq_norms();

        for (int i=0;i<m_lx;i++)
            for (int j=0;j<m_ly;j++) {
                E1 += preK[i][j] * preK[i][j] * hq_norm[i][j];
                E2 += preK[i][j] * sq_hq[i][j];
                E3 += preS[i][j] * hq_norm[i][j];
            }

        mem_E = -0.5 * E0; // initial value for Ising energy term
        coupling_E = -E2 / (double)m_grid.get_volume(); // initial value for coupling energy term

        E0 *= -0.5 * J_const;
        E1 *= 0.5 * KAPPA / (double)m_grid.get_volume();
        E2 *= -GAMMA / (double)m_grid.get_volume();
        E3 *= 0.5 * SIGMA / (double)m_grid.get_volume();

        return E0 + E1 + E2 + E3;

    }

    /** ******************************************************************* **/
    /** calculating the total energy of the system ************************ **/
    double get_energy() {

        double E0(0.), E1(0.), E2(0.), E3(0.);
        int spin, s_sum;
        double cell_height;
        double h_sum, curv;

        for (int idx=0;idx<m_grid.get_volume();idx++) {

            /** energy contribution of spin-spin interaction **/
            spin = m_grid.get_cellSpin(idx);
            s_sum = m_grid.get_sum_neighbors_spin(idx);
            E0 += (double)(spin*s_sum);

            /** energy contribution of bending term **/
            cell_height = m_grid.get_cellHeight(idx);
            h_sum = m_grid.get_sum_neighbors_height(idx);
            curv = h_sum - 4.*cell_height;
            E1 += curv*curv;

            /** energy contribution of spin-curvature coupling **/
            E2 += (double)(spin)*curv;

            /** energy contribution of tension term **/
            E3 += m_grid.get_delH2(idx);

        }

        mem_E = -0.5 * E0; // initial value for Ising energy term
        coupling_E = E2; // initial value for coupling energy term

        E0 *= -0.5 * J_const;
        E1 *= 0.5 * KAPPA;
        E2 *= GAMMA;
        E3 *= 0.5 * SIGMA;

        return E0 + E1 + E2 + E3;

    }

    /** ******************************************************************* **/
    /** Monte Carlo move                                                     */
    /** ******************************************************************* **/
    void MC_move() {

        double Pacc, delta_E(0.);

        int cell_index;
        do {
        cell_index = (int) (drandom()*(double)m_grid.get_volume());
        } while (fixed_number(cell_index));

        int spin = m_grid.get_cellSpin(cell_index);
        int s_sum = m_grid.get_sum_neighbors_spin(cell_index);

        delta_mem_E = 2.*spin*s_sum;
        delta_coupling_E = -2.*spin*curvature(cell_index);
        delta_E = J_const*delta_mem_E + GAMMA*delta_coupling_E;

        Pacc = exp(-m_beta*delta_E);
        if (Pacc>=1. || drandom()<Pacc) {
            m_grid.flip_spin(cell_index);
            m_orderParameter -= 2*spin;
            m_energy += delta_E;
            coupling_E += delta_coupling_E;
            mem_E += delta_mem_E;
        }

        if (drandom()<mc_criterion) {

            update_heights();
            m_energy = get_energy();

        }

        mc_counter++;

    }

    /** ******************************************************************* **/
    void update_heights() {
        hq_complex.resize(m_lx,vector<complex<double> > (nyp));

        sr_matrix = get_2dArraySpin();
        sq_complex = do_fft(sr_matrix);

        double mu_r, mu_i, sig;

        for (int i=0;i<m_lx;i++)
            for (int j=0;j<nyp;j++)
                if (i+j!=0) {

                    mu_r = pre_mu[i][j] * real(sq_complex[i][j]);
                    mu_i = pre_mu[i][j] * imag(sq_complex[i][j]);
                    sig = sig_hq[i][j];

                    vector<double> rg = gauss_random();

                    hq_complex[i][j] = complex<double> (mu_r+sig*rg[0], mu_i+sig*rg[1]);

                    if (j==0 && i>nxp) {
                        hq_complex[i][j].real(+hq_complex[m_lx-i][j].real());
                        hq_complex[i][j].imag(-hq_complex[m_lx-i][j].imag());
                    }

                }

        hq_complex[0][0] = complex<double> (0.,0.);

        if (m_lx%2==0 && m_ly%2==0) {
            hq_complex[nxp][0].imag(0.);
            hq_complex[0][nyp-1].imag(0.);
            hq_complex[nxp][nyp-1].imag(0.);
            for (int i=nxp+1;i<m_lx;i++) {
                hq_complex[i][nyp-1].real(+hq_complex[m_lx-i][nyp-1].real());
                hq_complex[i][nyp-1].imag(-hq_complex[m_lx-i][nyp-1].imag());
            }
        }

        if (m_lx%2==0 && m_ly%2!=0) {
            hq_complex[nxp][0].imag(0.);
        }

        if (m_lx%2!=0 && m_ly%2==0) {
            hq_complex[0][nyp-1].imag(0.);
            for (int i=nxp+1;i<m_lx;i++) {
                hq_complex[i][nyp-1].real(+hq_complex[m_lx-i][nyp-1].real());
                hq_complex[i][nyp-1].imag(-hq_complex[m_lx-i][nyp-1].imag());
            }
        }

        // do inverse Fourier transform to get the real heights

        do_ifft();

    }

    /** ******************************************************************* **/
    void set_pre_coeffs() {
        preK.resize(m_lx,vector<double > (m_ly));
        preS.resize(m_lx,vector<double > (m_ly));

        for (int i=0;i<m_lx;i++)
            for (int j=0;j<m_ly;j++) {

                preK[i][j] = 2.0 * (2.-cos(2.*M_PI*i/(double)m_lx)-cos(2.*M_PI*j/(double)m_ly));
                preS[i][j] = 0.5 * (2.-cos(4.*M_PI*i/(double)m_lx)-cos(4.*M_PI*j/(double)m_ly));

            }
    }

    /** ******************************************************************* **/
    void set_hq_coeffs() {

        vector<vector<double> > tmpF;
        tmpF.resize(m_lx,vector<double > (nyp));

        pre_mu.resize(m_lx,vector<double > (nyp));
        sig_hq.resize(m_lx,vector<double > (nyp));

        for (int i=0;i<m_lx;i++)
            for (int j=0;j<nyp;j++) {

                tmpF[i][j] = 1.0/(KAPPA*preK[i][j]*preK[i][j]+SIGMA*preS[i][j]);

                pre_mu[i][j] = GAMMA*preK[i][j]*tmpF[i][j];
                sig_hq[i][j] = sqrt(0.5*tmpF[i][j]*(double)(m_grid.get_volume())/m_beta);

            }

        if (m_lx%2==0 && m_ly%2==0) {
            sig_hq[nxp][0] *= sqrt(2.);
            sig_hq[0][nyp-1] *= sqrt(2.);
            sig_hq[nxp][nyp-1] *= sqrt(2.);
        }

        if (m_lx%2==0 && m_ly%2!=0) {
            sig_hq[nxp][0] *= sqrt(2.);
        }

        if (m_lx%2!=0 && m_ly%2==0) {
            sig_hq[0][nyp-1] *= sqrt(2.);
        }

    }

    /** ******************************************************************* **/
    bool fixed_number(int idx) {
        double frac = 0.1;

        int minOP = -(int)(frac * (double)m_grid.get_volume());
        int maxOP = +(int)(frac * (double)m_grid.get_volume());
        minOP += initial_OP;
        maxOP += initial_OP;
        int cell_state = m_grid.get_cellSpin(idx);
        int temp_orderParameter = m_orderParameter - 2*cell_state;
        if (minOP <= temp_orderParameter && temp_orderParameter <= maxOP)
            return false;
        else
            return true;
    }

    /** ******************************************************************* **/
    double curvature(int cell_index) {
        double result(0.);
        double cell_height = m_grid.get_cellHeight(cell_index);
        double h_sum = m_grid.get_sum_neighbors_height(cell_index);
        result = (h_sum - 4.*cell_height);
        return result;
    }

    /** ========================================================= **/
    double curvature9(int cell_index) {
        double result(0.);
        double cell_height = m_grid.get_cellHeight(cell_index);
        double h_sum = m_grid.get_sum_neighbors_height(cell_index);
        double h1_sum = m_grid.get_sum_next1_neighbors_height(cell_index);
        result = (4.*h_sum + h1_sum - 20.*cell_height);
        return result;
    }

    /** performing KAWASAKI Monte Carlo ************************************/
    void Kawasaki() {

        double Pacc, delta_E(0.);
        int idx[2], spin[2], s_sum[2];
        double curv[2];

        do {
            idx[0] = (int) (drandom()*(double)m_grid.get_volume());
            idx[1] = (int) (drandom()*(double)m_grid.get_volume());

            spin[0] = m_grid.get_cellSpin(idx[0]);
            spin[1] = m_grid.get_cellSpin(idx[1]);

            s_sum[0] = m_grid.get_sum_neighbors_spin(idx[0]);
            s_sum[1] = m_grid.get_sum_neighbors_spin(idx[1]);
        }
        while (spin[0]==spin[1] || s_sum[0]==s_sum[1]);

        curv[0] = curvature(idx[0]);
        curv[1] = curvature(idx[1]);

        delta_E = J_const*(spin[1]-spin[0])*(s_sum[1]-s_sum[0]);
        delta_E += GAMMA*(spin[1]-spin[0])*(curv[0]-curv[1]);

        Pacc = exp(-m_beta*delta_E);
        if (Pacc>=1 || drandom()<Pacc) {
            m_grid.swap_spins(idx[0],idx[1]);
            m_energy += delta_E;
        }

        if (drandom()<mc_criterion) {

            update_heights();
            m_energy = get_energy();

        }

        mc_counter++;

    }

    /** ******************************************************************* **/
    /** ******************************************************************* **/
    vector<vector<double> > get_2dArraySpin() {
        vector<vector<double> > result;
        result.resize(m_lx,vector<double> (m_ly));
        for (int i=0;i<m_lx;i++)
            for (int j=0;j<m_ly;j++)
                result[i][j] = (double)m_grid.get_cellSpin(i,j);
        return result;
    }

    /** ******************************************************************* **/
    vector<vector<double> > get_2dArrayHeight() {
        vector<vector<double> > result;
        result.resize(m_lx,vector<double> (m_ly));
        for (int i=0;i<m_lx;i++)
            for (int j=0;j<m_ly;j++)
                result[i][j] = (double)m_grid.get_cellHeight(i,j);
        return result;
    }

    /** Fast Fourier Transform ******************************************** **/
    vector<vector<complex<double> > > do_fft(vector<vector<double> > ff_vector) {
        size_t align=sizeof(Complex);

        array2<double> f_array(m_lx,m_ly,align);
        array2<Complex> g_array(m_lx,nyp,align);

        rcfft2d Forward(m_lx,m_ly,f_array,g_array);

        for (int i=0; i < m_lx; i++)
            for (int j=0; j < m_ly; j++)
                f_array(i,j) = ff_vector[i][j];

        Forward.fft(f_array,g_array);

        vector<vector<complex<double> > > result;
        result.resize(m_lx, vector<complex <double> > (nyp));

        for (int i=0;i<m_lx;i++)
            for (int j=0;j<nyp;j++)
                result[i][j] = g_array(i,j);
        return result;
    }

    /** inverse Fast Fourier Transform ************************************ **/
    void do_ifft() {
//         vector<vector<complex<double> > > gg_vecotr;
//         gg_vecotr = g_vector;
        size_t align=sizeof(Complex);

        array2<double> f_array(m_lx,m_ly,align);
        array2<Complex> g_array(m_lx,nyp,align);

        crfft2d Backward(m_lx,m_ly,g_array,f_array);

        for (int i=0; i < m_lx; i++)
            for (int j=0; j < nyp; j++)
                g_array(i,j) = hq_complex[i][j];

        Backward.fftNormalized(g_array,f_array);

//         vector<vector<double> > result;
//         result.resize(m_lx, vector<double> (m_ly));

        for (int i=0;i<m_lx;i++)
            for (int j=0;j<m_ly;j++)
                m_grid.set_cellHeight(i, j, f_array(i,j));

//         return result;
    }

    /** updating Fourier amplitudes of the height & spin profile ********** **/
    void update_sq_hq_norms() {
        hq_complex.resize(m_lx,vector<complex<double> > (nyp));
        hq_norm.resize(m_lx,vector<double > (m_ly));
        sq_hq.resize(m_lx,vector<double > (m_ly));

        hr_matrix = get_2dArrayHeight();
        hq_complex = do_fft(hr_matrix);

        sr_matrix = get_2dArraySpin();
        sq_complex = do_fft(sr_matrix);

        for (int i=0;i<m_lx;i++)
            for (int j=0;j<nyp;j++) {
                hq_norm[i][j] = norm(hq_complex[i][j]);
                sq_hq[i][j] = real(conj(sq_complex[i][j])*hq_complex[i][j]);
            }
        /** *************************************************************** **/
        for (int i=0;i<m_lx;i++)
            for (int j=nyp;j<m_ly;j++) {
                hq_norm[i][j] = hq_norm[i][m_ly-j];
                sq_hq[i][j] = sq_hq[i][m_ly-j];
            }
    }


    /** ******************************************************************* **/
    /** BOX MULLER algorithm                                                **/
    /** random number from standard gaussian distribution                   **/
    /** P(x) ~ exp(-0.5*x*x); mean = 0 and SD = 1                           **/
    /** ******************************************************************* **/
    vector<double> gauss_random(void) {
        vector<double> result;
        double r1,r2,tmp;

        r1 = 1.0 - drandom();
        r2 = 1.0 - drandom();
        tmp = sqrt(-2.0 * log(r1));
        result.push_back(tmp * cos(2.0*M_PI*r2));
        result.push_back(tmp * sin(2.0*M_PI*r2));
        return result;
    }

    /** ******************************************************************* **/
    double log_gauss_random(double x)
    {
        return -0.5*x*x;
    }

    /** ******************************************************************* **/
    /** setting parameters ************************************************ **/
    void set_beta(const double t_beta) {
        m_beta = t_beta;
    }

    void set_temperature(const double t_temperature) {
        m_temperature = t_temperature;
        m_beta = 1./(Boltzmann_const*t_temperature);
    }

    void set_mesh_size(const double t_mesh) {
        mesh_size = t_mesh;
    }

    void set_J_const(const double t_J) {
        J_const = t_J;
    }

    void set_kappa(const double t_kappa) {
        KAPPA = t_kappa/(mesh_size*mesh_size);
    }

    void set_sigma(const double t_sigma) {
        SIGMA = t_sigma;
    }

    void set_gamma(const double t_gamma) {
        GAMMA = t_gamma;
    }

    /** ******************************************************************* **/
    void set_spin_composition(double _sd, double _su) {
        spinD_frac = _sd;
        spinU_frac = _su;
        initial_OP = (int)((spinD_frac*(double)SpinD + spinU_frac*(double)SpinU) *
                           (double)m_grid.get_volume());
    }

    /** ******************************************************************* **/
    /** height moments **************************************************** **/
    vector<double> height_moments() {
        vector<double> result;

        double h1(0.), h2(0.);
        for (int idx=0;idx<m_grid.get_volume();idx++) {
            double tmpH = m_grid.get_cellHeight(idx);
            h1 += tmpH;
            h2 += tmpH*tmpH;
        }
        h1 /= (double)m_grid.get_volume();
        h2 /= (double)m_grid.get_volume();
        result.push_back(h1);
        result.push_back(h2);

        return result;
    }

    /** ******************************************************************* **/
    vector<vector<double > > get_sq_norm() {
        sq_norm.resize(m_lx,vector<double > (m_ly));

        sr_matrix = get_2dArraySpin();
        sq_complex = do_fft(sr_matrix);

        for (int i=0;i<m_lx;i++)
            for (int j=0;j<nyp;j++) {
                sq_norm[i][j] = norm(sq_complex[i][j]);
            }
        /** *************************************************************** **/
        for (int i=0;i<m_lx;i++)
            for (int j=nyp;j<m_ly;j++) {
                sq_norm[i][j] = sq_norm[i][m_ly-j];
            }

        return sq_norm;
    }

};

#endif
