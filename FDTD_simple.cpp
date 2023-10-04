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

const double pi = 3.14159265358979323846;
const double eps0 = 8.8541878128e-12;
const double mu0 = 1.256637062e-6;
const double c0 = 1 / sqrt(eps0 * mu0);
const double imp0 = sqrt(mu0 / eps0);
const double tau = 153e-12;
const double eps_r = 1;
const double dt = 1.768e-12;
const double dz = 1.1e-3;

const int imax = 1500;  // max discretization in space
const int nmax = 10000; // max discretization in time
const int isource = 750;
const int i_media_start = 800;
const int i_media_stop = 830;

double mur_boundary_const = (c0 * dt - dz) / (c0 * dt + dz);


int main()
{
    /* update parameters */
    complex<double> Ex[imax] = {0};
    complex<double> Hy[imax - 1] = {0};
    complex<double> Ex_prev[imax] = {0};
    complex<double> Hy_prev[imax - 1] = {0};

    /* PML vaiables */
    int PML_length = 50;
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

    for (int n = 0; n < nmax; n++)
    {
        // update electric field
        for (int i = 1; i < imax - 1; i++)
        {
            // Ex[i] = Ex_prev[i] + dt/(dz*eps0) * (Hy_prev[i] - Hy_prev[i-1]);  // air media
            Ex[i] = (Ex_prev[i]/dt * (1 - PML_const_E[i]) - 1.0/(dz*eps_r*eps0) * (Hy_prev[i] - Hy_prev[i - 1]))/(1.0/dt * (1 + PML_const_E[i])); 
        }
        // update electric boundaries
        Ex[0] = Ex_prev[1] + mur_boundary_const * (Ex[1] - Ex_prev[0]);
        Ex[imax - 1] = Ex_prev[imax - 2] + mur_boundary_const * (Ex[imax - 2] - Ex_prev[imax - 1]);
        // electric field source
        Ex[isource - 1] = Ex[isource - 1] + complex<double>(source_function(n), 0.00);
        // update previous time-step of electric field
        for (int i = 0; i < imax; i++)
        {
            Ex_prev[i] = Ex[i];
        }

        // update magnetic field
        for (int i = 0; i < imax - 1; i++)
        {
            Hy[i] = Hy_prev[i] - dt *(PML_const_H[i] * Hy_prev[i] + 1/(dz*mu0) * (Ex[i + 1] - Ex[i]));
            Hy_prev[i] = Hy[i];
        }


        /* save Ex fields at position 1 and 2 */
        Ex_pos1[n] = Ex[i_media_start];
        Ex_pos2[n] = Ex[i_media_stop];
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
