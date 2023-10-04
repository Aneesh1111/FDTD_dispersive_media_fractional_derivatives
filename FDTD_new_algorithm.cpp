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

/* constants */
const double pi = 3.14159265358979323846;
const double eps0 = 8.8541878128e-12;
const double mu0 = 1.256637062e-6;
const double c0 = 1 / sqrt(eps0 * mu0);
const double imp0 = sqrt(mu0 / eps0);
const double eps_s = 50;//50;
const double eps_inf = 2;//2;
const double d_eps = eps_s - eps_inf;

/* input parameters */
const int N = 1;
const int M = 3;
const double tau = 318e-12;
const double alpha[N] = { 0.2 };
const double beta1[M] = { 0.3, 0.5, 0.9 };
double a[N] = { 1.0 };
double b[M] = { 9.0, 2.0, 10.0 };  /* the tau^power get set in main */
const double eps_rs = 60;
const double conductivity = 0.0;

/* discretization */
const double dt = 1.768e-12;
const double dz = 1.1e-3;
const int gpof_max = 8;
const int imax = 1500;  // max discretization in space
unsigned long long int nmax = 10000;  // max discretization in time
const int isource = 700;
const int i_media_start = 800;
const int i_media_stop = 830;

/* boundary condition constant */
double mur_boundary_const = (c0 * dt - dz) / (c0 * dt + dz);


int main()
{
    /* set a_k and b_l constants */
    for (int i = 0; i < N; i++)
    {
        a[i] = a[i]*pow(tau, alpha[i]);
    }
    for (int i = 0; i < M; i++)
    {
        b[i] = b[i]*pow(tau, beta1[i]);
    }

    /* imported GPOF coefficients */
    double a_gpof_alpha1[gpof_max] = { 1.917075539048859, -0.560287204184830, -0.198639691967215, -0.085752446603181, -0.045331304427383, -0.018930049994618, -0.006964496770346, -0.001170344468412 };
    complex<double> b_gpof_alpha1[gpof_max] = {{-5.347160450874707, 3.141592653589793}, {-2.579206230697158, 0.00}, {-1.356262901949623, 0.00}, {-0.725900394409128, 0.00}, {-0.359939622406301, 0.00}, {-0.151364238040643 , 0.00}, {-0.045973859900755, 0.00}, {-0.006388607475665, 0.00}};
    double a_gpof_beta1[gpof_max] = {2.669664779268722, -1.130353379179308, -0.319124863825985, -0.127427146917053, -0.060073567274821, -0.023621097952209, -0.007853612584541, -0.001211110983464};
    complex<double> b_gpof_beta1[gpof_max] = {{-5.016172362040061 , 3.141592653589793}, {-2.708543106502772, 0.0000}, {-1.413642612584595, 0.0000}, {-0.760037597871560, 0.0000}, {-0.380999711174307 , 0.0000}, {-0.163439794019462, 0.0000}, {-0.051602128858320, 0.0000}, {-0.007889929899061, 0.00}};
    double a_gpof_beta2[gpof_max] = {5.331832122564569, -3.463142726299982, -0.575971123489243, -0.187842028805195, -0.072239715677446, -0.024849262526770, -0.006883150744650, -9.041148289020017e-04 };
    complex<double> b_gpof_beta2[gpof_max] = {{-4.682192687845149 , 3.141592653589793}, {-3.015762892824395, 0.00}, {-1.535975941329309 , 0.00}, {-0.831468861641271, 0.00}, {-0.425035368945261, 0.00}, {-0.188985508593928, 0.00}, {-0.063884248848722, 0.00}, {-0.011465455385803, 0.00}};
    double a_gpof_beta3[gpof_max] = {34.711296294933700, -33.087263703948906, -0.496500152254662, -0.094861637251343, -0.024878882072776, -0.006381535486200, -0.001285233632307, -1.251502861917283e-04};
    complex<double> b_gpof_beta3[gpof_max] = {{-4.830293886144696, 3.141592653589793}, {-4.220728932009782, 0.0000}, {-1.819457273975490, 0.0000}, {-0.988060793726207, 0.0000}, {-0.520664832308133, 0.0000}, {-0.245289842399021, 0.0000}, {-0.092246963047317, 0.0000}, {-0.020870775550631, 0.0000}};
    
    /* cole-cole parameters */
    // double a_gpof_alpha1[gpof_max] = { 0.0 };
    // complex<double> b_gpof_alpha1[gpof_max] = { 0.0 };
    // double a_gpof_beta1[gpof_max] = {34.711296294933700, -33.087263703948906, -0.496500152254662, -0.094861637251343, -0.024878882072776, -0.006381535486200, -0.001285233632307, -1.251502861917283e-04};
    // complex<double> b_gpof_beta1[gpof_max] = {{-4.830293886144696, 3.141592653589793}, {-4.220728932009782, 0.0000}, {-1.819457273975490, 0.0000}, {-0.988060793726207, 0.0000}, {-0.520664832308133, 0.0000}, {-0.245289842399021, 0.0000}, {-0.092246963047317, 0.0000}, {-0.020870775550631, 0.0000}};


    /* Precalculated GPOF coefficients */
    double a_gpof_P[gpof_max*M] = { 0 };
    complex<double> b_gpof_P[gpof_max*M] = { 0 };
    double a_gpof_E[N*gpof_max] = { 0 };
    complex<double> b_gpof_E[N*gpof_max] = { 0 };

    for (int j = 0; j < gpof_max; j++)
    {
        a_gpof_E[j] = a_gpof_alpha1[j];
        b_gpof_E[j] = b_gpof_alpha1[j];
        a_gpof_P[j] = a_gpof_beta1[j];
        b_gpof_P[j] = b_gpof_beta1[j];
        a_gpof_P[j+gpof_max] = a_gpof_beta2[j];
        b_gpof_P[j+gpof_max] = b_gpof_beta2[j];
        a_gpof_P[j+2*gpof_max] = a_gpof_beta3[j];
        b_gpof_P[j+2*gpof_max] = b_gpof_beta3[j];
    }

    /* initialise derived constants */
    double C1 = 0;
    double C2 = 0;
    double C3 = 0;
    double C4 = 0;
    double C5 = 0;

    for (int i = 0; i < N; i++)
    {
        C2 = C2 + a[i]/pow(dt, alpha[i]);
    }

    for (int i = 0; i < M; i++)
    {
        C1 = C1 + b[i]/pow(dt, beta1[i]);
    }

    C5 = conductivity/2.0 + eps_inf*eps0/dt;
    C4 = conductivity/2.0 - eps_inf*eps0/dt;
    C3 = eps0*eps_rs*(1.0 + C2)/C5;

    /* update parameters */
    complex<double> vj_E[gpof_max*N][imax] = {0};
    complex<double> vj_P[gpof_max*M][imax] = {0};
    complex<double> vj_E_sum[N][imax] = {0};
    complex<double> vj_P_sum[M][imax] = {0};
    complex<double> psi_lP[M][imax] = {0};
    complex<double> psi_kE[N][imax] = {0}; 
    complex<double> psi_lP_sum[imax] = {0};
    complex<double> psi_kE_sum[imax] = {0};
    complex<double> Ex[imax] = {0};
    complex<double> P[imax] = {0};
    complex<double> C[imax] = {0};
    complex<double> D[imax] = {0};
    complex<double> Hy[imax - 1] = {0};
    complex<double> vj_E_prev[gpof_max*N][imax] = {0};
    complex<double> vj_P_prev[gpof_max*M][imax] = {0};
    complex<double> Ex_prev[imax] = {0};
    complex<double> P_prev[imax] = {0};
    complex<double> Hy_prev[imax - 1] = {0};

    complex<double> Ex_at_1000[imax] = { 0 };
    complex<double> Ex_at_2000[imax] = { 0 };
    complex<double> Ex_at_4000[imax] = { 0 };


    /* PML vaiables */
    int PML_length = 0;
    double sigma[imax] = {0};
    double PML_const_E[imax] = { 0 };
    double PML_const_H[imax] = { 0 };

    for (int i = 0; i < imax; i++)
    {
        if (i < PML_length) 
        {
            sigma[i] = (2.0*eps0)/(3.0*dt) * pow((PML_length-i)/PML_length, 3);
        }
        else if (i > imax-PML_length)
        {
            sigma[i] = (2.0*eps0/(3.0*dt)) * pow((i-(imax-PML_length))/PML_length, 3);
        }
    }
    for (int i = 0; i < imax; i++) 
    {
        PML_const_E[i] = sigma[i]*dt/(2.0*eps0);
        PML_const_H[i] = sigma[i]/eps0;
    }

    /* fields to be saved/stored */
    complex<double> Ex_pos1[nmax] = { 0 };
    complex<double> Ex_pos2[nmax] = { 0 };
    // complex<double> Ex_source[nmax] = { 0 };

    /* FDTD scheme */
    for (int n = 0; n < nmax; n++)
    {
        /* update running sums (auxilliary vectors) */
        for (int i = i_media_start; i < i_media_stop+1; i++)
        {
            for (int ii = 0; ii < N; ii++)
            {
                for (int j = gpof_max*ii; j < gpof_max*(ii + 1); j++)
                {
                    vj_E[j][i] = a_gpof_E[j] * exp(b_gpof_E[j]) * Ex_prev[i] + exp(b_gpof_E[j]) * vj_E_prev[j][i];
                    vj_E_prev[j][i] = vj_E[j][i];
                    vj_E_sum[ii][i] = vj_E_sum[ii][i] + vj_E[j][i];
                }
            }

            for (int ii = 0; ii < M; ii++)
            {
                for (int j = gpof_max*ii; j < gpof_max*(ii + 1); j++)
                {
                    vj_P[j][i] = a_gpof_P[j] * exp(b_gpof_P[j]) * P_prev[i] + exp(b_gpof_P[j]) * vj_P_prev[j][i];
                    vj_P_prev[j][i] = vj_P[j][i];
                    vj_P_sum[ii][i] = vj_P_sum[ii][i] + vj_P[j][i];
                }
            }

            for (int l = 0; l < M; l++)
            {
                psi_lP_sum[i] = psi_lP_sum[i] + b[l]/pow(dt, beta1[l]) * vj_P_sum[l][i];
            }
            for (int k = 0; k < N; k++)
            {
                psi_kE_sum[i] = psi_kE_sum[i] + a[k]/pow(dt, alpha[k]) * vj_E_sum[k][i];
            }            
        }

        /* update polarisation */
        for (int i = 1; i < imax-1; i++)
        {
            P[i] = 1.0/(1.0 + C1 + C3/dt) * ( C3 * ( (Hy_prev[i] - Hy_prev[i-1])/dz - C4*Ex_prev[i] + P_prev[i]/dt ) - psi_lP_sum[i] + eps0*eps_rs*psi_kE_sum[i]);
        }

        /* update electric field */
        /* air */
        for (int i = 1; i < i_media_start; i++)
        {
            Ex[i] = Ex_prev[i]*(1 - PML_const_E[i])/(1 + PML_const_E[i]) + dt/(dz*eps0)*(Hy_prev[i] - Hy_prev[i-1])/(1 + PML_const_E[i]);
        }
        /* media */
        for (int i = i_media_start; i < i_media_stop+1; i++)
        {
            Ex[i] = 1.0/C5 * ( (Hy_prev[i] - Hy_prev[i-1])/dz - C4*Ex_prev[i] - (P[i] - P_prev[i])/dt );
            // Ex[i] = Ex_prev[i]*(1 - PML_const_E[i])/(1 + PML_const_E[i]) + dt/(dz*eps0)*(Hy_prev[i] - Hy_prev[i-1])/(1 + PML_const_E[i]);
        }
        /* air */
        for (int i = i_media_stop+1; i < imax-1; i++)
        {
            Ex[i] = Ex_prev[i]*(1 - PML_const_E[i])/(1 + PML_const_E[i]) + dt/(dz*eps0)*(Hy_prev[i] - Hy_prev[i-1])/(1 + PML_const_E[i]);
        }

        // // update electric field
        // for (int i = 1; i < imax - 1; i++)
        // {
        //     Ex[i] = A * (eps0 * eps_inf * Ex_prev[i] * (1 - PML_const_E[i])/(1 + PML_const_E[i]) + P_prev[i] + qq * vj_sum[i] - (dt / dz) * (Hy_prev[i] - Hy_prev[i - 1])/(1 + PML_const_E[i])); // rekanos
        // }

        /* update electric field boundaries */
        Ex[0] = Ex_prev[1] + mur_boundary_const * (Ex[1] - Ex_prev[0]);
        Ex[imax - 1] = Ex_prev[imax - 2] + mur_boundary_const * (Ex[imax - 2] - Ex_prev[imax - 1]);
        /* electric field source */
        Ex[isource - 1] = Ex[isource - 1] + complex<double>(source_function(n), 0.00);
        /* update previous time-step of electric field */
        for (int i = 0; i < imax; i++)
        {
            Ex_prev[i] = Ex[i];
        }

        /* update magnetic field */
        for (int i = 0; i < imax-1; i++)
        {
            Hy[i] = Hy_prev[i]*(1 - PML_const_H[i]*dt) + dt/(dz*mu0) * (Ex[i + 1] - Ex[i]);
            Hy_prev[i] = Hy[i];
        }

        // // update magnetic field
        // for (int i = 0; i < imax - 1; i++)
        // {
        //     Hy[i] = Hy_prev[i] - dt *(PML_const_H[i] * Hy_prev[i] + 1/(dz*mu0) * (Ex[i + 1] - Ex[i]));
        //     Hy_prev[i] = Hy[i];
        // }

        /* reset auxilliary vector summations */
        for (int i = 0; i < imax; i++)
        {
            for (int l = 0; l < M; l++)
            {
                vj_P_sum[l][i] = 0;
            }
            for (int k = 0; k < N; k++)
            {
                vj_E_sum[k][i] = 0;
            }
            psi_lP_sum[i] = 0;
            psi_kE_sum[i] = 0;
            /* update P */
            P_prev[i] = P[i];
        }

        /* save Ex fields at position 1 and 2 */
        Ex_pos1[n] = Ex[i_media_start];
        Ex_pos2[n] = Ex[i_media_stop];

        if(n == 500) {
            for(int i = 0; i < imax; i++)
            {
                Ex_at_1000[i] = Ex[i];
            }
        }
        if(n == 1000) {
            for(int i = 0; i < imax; i++)
            {
                Ex_at_2000[i] = Ex[i];
            }
        }
        if(n == 1500) {
            for(int i = 0; i < imax; i++)
            {
                Ex_at_4000[i] = Ex[i];
            }
        }
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
    for (int i = 0; i < nmax; i++)
    {
        file_Ex_pos1_pos2 << Ex_pos1[i].real() << " ";
        file_Ex_pos1_pos2 << Ex_pos1[i].imag() << " ";
        file_Ex_pos1_pos2 << Ex_pos2[i].real() << " ";
        file_Ex_pos1_pos2 << Ex_pos2[i].imag() << endl;
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
        wave_in_media << Ex_at_4000[i].imag() << endl;
    }

    cout << "Done!" << endl;

    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Time = " << duration.count() << " ms" << endl;

    return 0;
}

double source_function(int t)
{
    double fc = 6e9;
    double fbw = 5.0 / 6.0;
    double bwr = -3.0;
    double n = (t - 150) * dt;

    double fraction_of_max_peak = pow(10, bwr / 20.0);
    double frequency_variance = -fbw * fbw * fc * fc / (8 * log(fraction_of_max_peak));
    double time_variance = 1 / (4 * pi * pi * frequency_variance);

    double exponent_of_exp = -n * n / (2 * time_variance);
    double time_domain_pulse_envelope = exp(exponent_of_exp);

    double quadrature_phase_modulated_envelope = time_domain_pulse_envelope * sin(2 * fc * n * pi);
    return quadrature_phase_modulated_envelope;
}
