#include <iostream>
#include <cmath>
#include <complex>  
#include <fstream> 
#include "/home/aneesh/Desktop/Graduation Project Antenna Company/Cpp_code/Eigen/Dense"
#include "/home/aneesh/Desktop/Graduation Project Antenna Company/Cpp_code/functions.h"
#include <chrono>

using namespace Eigen;
using namespace std;
using namespace std::chrono;

/* start timing the code */
auto start = high_resolution_clock::now();

double pi = 3.14159265358979323846;
double eps0 = 8.8541878128e-12;
double mu0 = 1.256637062e-6;
double c0 = 1/sqrt(eps0*mu0);
double imp0 = sqrt(mu0/eps0);
double tau = 153e-12;
double eps_s = 50;
double eps_inf = 2;
double d_eps = eps_s - eps_inf;
double alpha = 0.9;
double dt = 1.768e-12;
double dz = 1.1e-3;

int M = 8;
int imax = 1500;   // max discretization in space
int nmax = 10000; // max discretization in time
int isource = 750;
int i_media_start = 800;
int i_media_stop = 830;

/* defined constants */
double q = pow(tau/dt, alpha);
double A = 1/(eps0*eps_inf + eps0*d_eps/(1 + q));
double qq = q/(1 + q);
double qq_eps = eps0*d_eps/(1 + q);
double tz_mu = dt/(mu0*dz);
/* boundary condition constant */
double mur_boundary_const = (c0*dt - dz)/(c0*dt + dz);

// eps_dielectric = eps0*ones(1, i_media_stop-i_media_start+2)*10;


int main() 
{
    double a_gpof[M] = {34.711296294933700, -33.087263703948906, -0.496500152254662, -0.094861637251343, -0.024878882072776, -0.006381535486200, -0.001285233632307, -1.251502861917283e-04};
    complex<double> b_gpof[M] = {{-4.830293886144696, 3.141592653589793}, {-4.220728932009782, 0.0000}, {-1.819457273975490, 0.0000}, {-0.988060793726207, 0.0000}, {-0.520664832308133, 0.0000}, {-0.245289842399021, 0.0000}, {-0.092246963047317, 0.0000}, {-0.020870775550631, 0.0000}};
    
    /* update parameters */
    complex<double> vj[M][imax] = { 0 };
    complex<double> vj_sum[imax] = { 0 };
    complex<double> Ex[imax] = { 0 };
    complex<double> P[imax] = { 0 };
    complex<double> Hy[imax-1] = { 0 };
    complex<double> vj_prev[M][imax] = { 0 };
    complex<double> Ex_prev[imax] = { 0 };
    complex<double> P_prev[imax] = { 0 };
    complex<double> Hy_prev[imax-1] = { 0 };

    /* fields to be saved/stored */
    complex<double> Ex_pos1[nmax] = { 0 };
    complex<double> Ex_pos2[nmax] = { 0 };
    complex<double> Ex_source[nmax] = { 0 };
    complex<double> Ex_at_1000[imax] = { 0 };
    complex<double> Ex_at_2000[imax] = { 0 };
    complex<double> Ex_at_4000[imax] = { 0 };
    complex<double> Ex_at_6000[imax] = { 0 };
    complex<double> Ex_at_8000[imax] = { 0 };
    complex<double> Ex_at_10000[imax] = { 0 };



    for(int n = 0; n < nmax; n++)
    {
        
        //update vj
        for(int i = 0; i < imax; i++)
        {
            for(int j = 0; j < M; j++)
            {
                vj[j][i] = a_gpof[j]*exp(b_gpof[j])*P_prev[i] + exp(b_gpof[j])*vj_prev[j][i];
                vj_prev[j][i] = vj[j][i];
                vj_sum[i] = vj_sum[i] + vj[j][i];
            }
        }

        // update electric field
        for(int i = 1; i < imax-1; i++)
        {
            // Ex[i] = Ex_prev[i] + dt/(dz*eps0) * (Hy_prev[i] - Hy_prev[i-1]);  // air media
            Ex[i] = A * (eps0*eps_inf*Ex_prev[i] + P_prev[i] + qq*vj_sum[i] + (dt/dz)*(Hy_prev[i] - Hy_prev[i-1]));  // rekanos
        }
        // update electric boundaries
        Ex[0] = Ex_prev[1] + mur_boundary_const * (Ex[1] - Ex_prev[0]);
        Ex[imax-1] = Ex_prev[imax-2] + mur_boundary_const * (Ex[imax-2] - Ex_prev[imax-1]);
        // electric field source
        Ex[isource-1] = Ex[isource-1] + complex<double>(source_function(n), 0.00);
        // update previous time-step of electric field
        for(int i = 0; i < imax; i++)
        {
            Ex_prev[i] = Ex[i];
        }

        // update polarization vector
        for(int i = 0; i < imax; i++)
        {
            P[i] = qq_eps*Ex[i] - qq*vj_sum[i];
            P_prev[i] = P[i];
        }
        
        // update magnetic field
        for(int i = 0; i < imax-1; i++)
        {
            Hy[i] = Hy_prev[i] + tz_mu * (Ex[i+1] - Ex[i]);
            Hy_prev[i] = Hy[i];
        }

        // reset vj_sum
        for(int i = 0; i < imax; i++)
        {
            vj_sum[i] = 0;
        }

        /* save Ex fields at position 1 and 2 */
        Ex_pos1[n] = Ex[i_media_start];
        Ex_pos2[n] = Ex[i_media_stop];

        // if(n == 1000) {
        //     for(int i = 0; i < imax; i++)
        //     {
        //         Ex_at_1000[i] = Ex[i];
        //     }
        // }
        // if(n == 2000) {
        //     for(int i = 0; i < imax; i++)
        //     {
        //         Ex_at_2000[i] = Ex[i];
        //     }
        // }
        // if(n == 4000) {
        //     for(int i = 0; i < imax; i++)
        //     {
        //         Ex_at_4000[i] = Ex[i];
        //     }
        // }
        // if(n == 6000) {
        //     for(int i = 0; i < imax; i++)
        //     {
        //         Ex_at_6000[i] = Ex[i];
        //     }
        // }
        // if(n == 8000) {
        //     for(int i = 0; i < imax; i++)
        //     {
        //         Ex_at_8000[i] = Ex[i];
        //     }
        // }
        // if(n == 9999) {
        //     for(int i = 0; i < imax; i++)
        //     {
        //         Ex_at_10000[i] = Ex[i];
        //     }
        // }
    }

    /* stop timing the code */
    auto stop = high_resolution_clock::now();

    /* save fields in text folder */
    ofstream file_Ex_pos1_pos2("Ex_pos1_pos2.txt");
    if (!file_Ex_pos1_pos2) 
    {
        cout << "Could not open file";
        return -1;
    }
    for(int i = 0; i < nmax; i++) 
    {
        file_Ex_pos1_pos2 << Ex_pos1[i].real() << " ";
        file_Ex_pos1_pos2 << Ex_pos1[i].imag() << " ";
        file_Ex_pos1_pos2 << Ex_pos2[i].real() << " ";
        file_Ex_pos1_pos2 << Ex_pos2[i].imag() << endl;
    }

    /* check source function */
    ofstream source_func("source_func.txt");
    if (!source_func) 
    {
        cout << "Could not open file";
        return -1;
    }
    for(int i = 0; i < nmax; i++) 
    {
        Ex_source[i] = complex<double>(source_function(i), 0.00);
        source_func << Ex_source[i].real() << " ";
        source_func << Ex_source[i].imag() << endl;
    }

    /* check wave in media */
    ofstream wave_in_media("wave_in_media.txt");
    if (!wave_in_media) 
    {
        cout << "Could not open file";
        return -1;
    }
    for(int i = 0; i < imax; i++) 
    {
        wave_in_media << Ex_at_1000[i].real() << " ";
        wave_in_media << Ex_at_1000[i].imag() << " ";
        wave_in_media << Ex_at_2000[i].real() << " ";
        wave_in_media << Ex_at_2000[i].imag() << " ";
        wave_in_media << Ex_at_4000[i].real() << " ";
        wave_in_media << Ex_at_4000[i].imag() << " ";
        wave_in_media << Ex_at_6000[i].real() << " ";
        wave_in_media << Ex_at_6000[i].imag() << " ";
        wave_in_media << Ex_at_8000[i].real() << " ";
        wave_in_media << Ex_at_8000[i].imag() << " ";
        wave_in_media << Ex_at_10000[i].real() << " ";
        wave_in_media << Ex_at_10000[i].imag() << endl;
    }

    cout << "Done!" << endl;

    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Time = " << duration.count() << " ms" << endl;

    return 0;
}


double source_function(int t)
{
  double fc = 6e9;
  double fbw = 5.0/6.0;
  double bwr = -3.0;
  double n = (t-150)*dt;

  double fraction_of_max_peak = pow(10, bwr/20.0);
  double frequency_variance = -fbw*fbw*fc*fc/(8*log(fraction_of_max_peak));
  double time_variance = 1/(4*pi*pi*frequency_variance);

  double exponent_of_exp = -n*n/(2*time_variance);
  double time_domain_pulse_envelope = exp(exponent_of_exp);

  double quadrature_phase_modulated_envelope = time_domain_pulse_envelope * sin(2*fc*n*pi);
  return quadrature_phase_modulated_envelope;
}

