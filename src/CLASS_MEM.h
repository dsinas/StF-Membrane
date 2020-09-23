/*******************************************************************************
* 
* CLASS_MEM.h: A class for computing structure factor of the membrane 
* composition that is coupled to the local membrane curvature.
*
* This program uses fftw3 library and fft source code for c++.
*
* This program is made freely available with the understanding that every 
* copy of this file must include this header and that it comes without any 
* WITHOUT ANY WARRANTY.

*******************************************************************************/

#ifndef CLASS_MEM_H
#define CLASS_MEM_H

#include "CLASS_MC.h"

using namespace std;

/** definition of MembraneClass ******************************************* **/
class MembraneClass : public MC_SystemClass {

protected:

    double op1, op2;
    double ave_mem_E, ave_mem_E2;
    double susceptibility, specific_heat;
    double ave_kn, nrm_kn;
    double domain_size;
    char filename[50];

public:

    MembraneClass (const int _Lx, const int _Ly):
        MC_SystemClass(_Lx, _Ly)
        {
        }

    virtual ~MembraneClass() {}

    /** ******************************************************************* **/
    void get_data_files() {

        /** these parameters set experimentally up to the resourses and the accuracy **/
        int Nmax = 100; // default: 100
        int nSweep = 10000; // default: 100000

        /** *********************************** **/

        const int Kn_max = 1 + (int)(0.5*sqrt((double)(m_lx*m_lx + m_ly*m_ly)));
        double sq[Kn_max], nsq[Kn_max];

        /** *********************************** **/
        ave_mem_E = ave_mem_E2 = 0.;
        op1 = op2 = 0.;
        /** *********************************** **/

        for (int kn=0;kn<Kn_max;kn++) sq[kn] = nsq[kn] = 0.;

        const int XH = (int)(m_lx/2), YH = (int)(m_ly/2);

        /** initial MC moves to get equilibrium **/
        // for (int ave_i=0;ave_i<Nmax;ave_i++)
		    cout << "\nearly iterations to equilibrium ...\n" << endl;
		    for (int sweep=0;sweep<nSweep;sweep++)
		        for (int i=0;i<m_grid.get_volume();i++)
		            MC_move();

        /** *********************************** **/

        for (int ave_i=0;ave_i<Nmax;ave_i++) {

        	cout << "iteration: " << ave_i+1 << " of " << Nmax << "." <<endl;

            for (int sweep=0;sweep<nSweep;sweep++) {

                for (int i=0;i<m_grid.get_volume();i++)
                    MC_move();

                /** ************************************ **/
                double tmp_m = (double)abs(m_orderParameter)/(double)m_grid.get_volume();
                op1 += tmp_m;
                op2 += tmp_m*tmp_m;
                /** ************************************ **/
                ave_mem_E += mem_E;
                ave_mem_E2 += mem_E*mem_E;

            }

            /* =================================================== */
            /* computing the structure factor: S(k)                */
            /* =================================================== */

            vector<vector<double> > sq_norm_matrix = get_sq_norm();

            int kx, ky;

            for (int ii=0;ii<m_lx;ii++)
                for (int jj=0;jj<m_ly;jj++) {
                    kx = ii;
                    if (ii>XH) kx -= m_lx;
                    ky = jj;
                    if (jj>YH) ky -= m_ly;
                    int kn = (int)sqrt((double)(kx*kx + ky*ky));
                    sq[kn] += sq_norm_matrix[ii][jj];
                    nsq[kn] ++;
                }

            /* =================================================== */

            // cout<<ave_i<<endl;

        }

        /** ************************************** **/
        op1 /= (double)(nSweep*Nmax);
        op2 /= (double)(nSweep*Nmax);

        susceptibility = (op2-op1*op1)*(double)m_grid.get_volume();

        ave_mem_E /= (double)(nSweep*Nmax);
        ave_mem_E2 /= (double)(nSweep*Nmax);

        specific_heat = (ave_mem_E2-ave_mem_E*ave_mem_E)/(double)m_grid.get_volume();

        sprintf(filename,"memE_%d_j%.2f_g%.2f.dat",m_lx,J_const,GAMMA);
        ofstream ie_file(filename);
        ie_file<<GAMMA<<"\t"<<ave_mem_E<<"\t"<<ave_mem_E2<<"\t"<<specific_heat<<"\t";
        ie_file<<op1<<"\t"<<op2<<"\t"<<susceptibility<<"\t"<<J_const*ave_mem_E<<endl;
        ie_file.close();

        /** ************************************** **/
        for (int kn=0;kn<Kn_max;kn++)
            if (nsq[kn])
                sq[kn] /= (nsq[kn] * (double)Nmax * (double)m_grid.get_volume());

        // sq[0] = 0.;

        double kk_coef = 2.*M_PI/sqrt((double)m_grid.get_volume());

        sprintf(filename,"st_factor_%d_j%.2f_g%.2f.dat",m_lx,J_const,GAMMA);
        ofstream sf_file(filename);
        for (int kn=0;kn<Kn_max;kn++)
            sf_file<<(double)kn*kk_coef<<"\t"<<sq[kn]<<"\t"<<1./sq[kn]<<"\t"<<kn<<endl;
        sf_file.close();

        ave_kn = nrm_kn = 0.;
        for (int kn=1;kn<Kn_max;kn++) {
            ave_kn += (double)kn*sq[kn];
            nrm_kn += sq[kn];
        }
        ave_kn /= nrm_kn;

        domain_size = sqrt((double)m_grid.get_volume())/ave_kn;

        /** ************************************** **/

        sprintf(filename,"ds_%d_j%.2f_g%.2f.dat",m_lx,J_const,GAMMA);
        ofstream ds_file(filename);
        ds_file<<J_const<<"\t"<<KAPPA<<"\t"<<SIGMA<<"\t"<<GAMMA<<"\t"<<domain_size<<"\t";
		    ds_file<<ave_kn<<"\t"<<ave_kn*kk_coef<<endl;
        ds_file.close();

        print_data_file();
        print_height_lattice();
        print_cmps_lattice();
        print_cmps_height();

    }

    /** printing data file *********************************************** **/
    void print_data_file() {
        sprintf(filename,"sys_data_%d_j%.2f_g%.2f.dat",m_lx,J_const,GAMMA);
        ofstream data_file(filename);
        data_file<<"Lx = "<<m_lx<<endl;
        data_file<<"Ly = "<<m_ly<<endl;
        data_file<<"\nj = "<<J_const<<endl;
        data_file<<"kappa = "<<KAPPA<<endl;
        data_file<<"sigma = "<<SIGMA<<endl;
        data_file<<"gamma = "<<GAMMA<<"\n"<<endl;
        data_file<<"n* = "<<ave_kn<<endl;
        data_file<<"k* = "<<ave_kn*2.*M_PI/sqrt((double)m_grid.get_volume())<<endl;
        data_file<<"domain size = "<<domain_size<<"\n"<<endl;
        data_file<<"susceptibility = "<<susceptibility<<endl;
        data_file<<"specific_heat = "<<specific_heat<<endl;
        data_file.close();
    }

    /** printing height lattice ******************************************* **/
    void print_height_lattice() {
        sprintf(filename,"height_%d_j%.2f_g%.2f.dat",m_lx,J_const,GAMMA);
        ofstream height_file(filename);
        for (int i=0;i<m_lx;i++) {
            for (int j=0;j<m_ly;j++) {
                int cell_index = m_grid.get_cell_index(i,j);
                double cell_height = m_grid.get_cellHeight(cell_index);
                height_file<<cell_height<<"\t";
            }
            height_file<<endl;
        }
        height_file.close();
    }

    /** printing spin lattice ********************************************* **/
    void print_cmps_lattice() {
        sprintf(filename,"cmps_%d_j%.2f_g%.2f.dat",m_lx,J_const,GAMMA);
        ofstream cmps_file(filename);
        for (int i=0;i<m_lx;i++) {
            for (int j=0;j<m_ly;j++) {
                int cell_index = m_grid.get_cell_index(i,j);
                int cell_State = m_grid.get_cellSpin(cell_index);
                cmps_file<<cell_State<<"\t";
            }
            cmps_file<<endl;
        }
        cmps_file.close();
    }

    /** printing composition_height lattices ****************************** **/
    void print_cmps_height() {
        sprintf(filename,"lh_%d_j%.2f_g%.2f.xyz",m_lx,J_const,GAMMA);
        ofstream sh_file(filename);

        sh_file<<m_grid.get_volume()<< "\n" <<"Atoms"<<endl;

        int idx, spin;
        double height;

        for (int i=0; i<m_grid.get_lx(); ++i) {
            for (int j=0; j<m_grid.get_ly(); ++j) {
                idx = m_grid.get_cell_index(i,j);
                spin = (int)(0.5*(double)(3+m_grid.get_cellSpin(idx)));
                height = m_grid.get_cellHeight(idx);
                sh_file<<spin<<"\t"<<i<<"\t"<<j<<"\t"<<height<<endl;
            }
        }
        sh_file.close();
    }

};

#endif
