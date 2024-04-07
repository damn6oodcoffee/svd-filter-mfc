#pragma once 

#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <ctime>
#include <fstream>
using namespace std;

struct SignalParam
{
	double a1=0, f1=0, phi1=0;
	double a2 = 0, f2 = 0, phi2 = 0;
	double a3 = 0, f3 = 0, phi3 = 0;
};

bool GenerateSineSignal(double* s, int N, double fd, SignalParam& p);

void AddNoise(double* s, int N, double a_noise);

void CorrelogramPSD(double* r_s, double* P, int N, int L, double fd);

void Autocorr(double* s, double* r_s, int N, int L); //s - сигнал, r_s - автокорр. м-ца NxN

//int svd_hestenes(int m_m, int n_n, double * a, double * u, double * v, double * sigma);
int svd_double(int m_m, int n_n, double * a, double * u, double * v, double * sigma);
void PrintSigma(double* s, int N, int M, double* u, double* v);

void Filter(double* sigma, double* u, double* v, double* r_s, double thrhold, int L);

