#pragma once
#include "stdafx.h"
#include "signal.h"
#include <fstream>
using namespace std;




bool GenerateSineSignal(double* s, int N, double fd, SignalParam& p)

{

	double dt = 1 / fd;

	for (int i = 0; i < N; i++)
	{
		s[i] = p.a1*sin(2 * M_PI*p.f1*i*dt + p.phi1) + p.a2*sin(2 * M_PI*p.f2*i*dt + p.phi2) + p.a3*sin(2 * M_PI*p.f3*i*dt + p.phi3);

	}


	return true;

}
void AddNoise(double* s, int N, double L) // L = 10 lg(E_n/E_s)
{

	double* g = new double[N];
	srand(time(NULL));
	for (int i = 0; i < N; i++)
	{
		g[i] = 0;
		for (int k = 0; k < 12; k++)
		{
			g[i] += rand() - RAND_MAX / 2;

		}
	}






	double E_s = 0;
	double E_n = 0;

	for (int i = 0; i < N; i++)
	{
		E_s += s[i] * s[i];
		E_n += g[i] * g[i];
	}

	double Beta = sqrt(E_s / E_n * pow(10, L / 10));
	
	for (int i = 0; i < N; i++)
	{
		s[i] = s[i] + Beta*g[i];
	}

	delete[] g;
}

void CorrelogramPSD(double* r_s, double* P, int N, int L, double fd)
{
	double T = 1 / fd;

	for (int i = 0; i < N; i++)
	{
		P[i] = r_s[0];
		for (int j = 1; j < L; j++)
		{

			//P[i] += r_s[abs(j)] * cos(2 * M_PI*(i*1.0 / (N - 1))*fd*j*T);
			P[i] += 2*r_s[j] * cos(2 * M_PI*(i*1.0 / (N - 1))*fd*j*T);
			
		}
		P[i] *=T;
		//P[i] = abs(P[i])*T;
	}




}


void Autocorr(double* s, double* r_s, int N, int L) //s - сигнал, r_s - автокорр. м-ца NxN
{

	fstream fileObj;
	fileObj.open("r_s.txt", ios::out| ios::trunc);


	for (int m = 0; m < L; m++)//(int m = 0; m < Odr; m++)
	{
		r_s[m] = 0;
		for (int n = 0; n < N - m; n++)
		{
			r_s[m] += (1. / (N - m)) * s[n + m] * s[n];
		}
	}


	for (int i = 1; i < L; i++)
	{
		for (int j = 0; j < L; j++)
		{

			r_s[i*L + j] = r_s[abs(i - j)];

		}


	}


	for (int i = 0; i < L; i++)
	{
		for (int j = 0; j < L; j++)
		{
			//fileObj.precision(4);
			fileObj << r_s[i*L + j] << "  ";

		}
	//	fileObj << endl;

	}



}


void PrintSigma(double* s, int N, int M, double *u, double *v)
{


	for (int k = 0; k < N; k++)
	{
		s[k*N+k] = s[k];

	}

	for (int m = 0; m < N; m++)
	{
		for (int l = 0; l < M; l++)
		{
			if (m != l) s[m*M +l] = 0;

		}


	}




	fstream sigmafile;
	sigmafile.open("sigma.txt", ios::out | ios::trunc);;


	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			sigmafile.precision(4);
			sigmafile << u[i*N + j] << "\t";

		}
		sigmafile << endl;

	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			sigmafile.precision(4);
			sigmafile << s[i*N + j] << "\t";

		}
		sigmafile << endl;

	}


	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			sigmafile.precision(4);
			sigmafile << v[i*N + j] << "\t";

		}
		sigmafile << endl;

	}


	sigmafile << endl << endl;

	double * r_new = new double[N*N];

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			r_new[i*N + j] = 0;
			for (int k = 0; k < N; k++)
			{
				r_new[i*N + j] += u[i*N + k] * s[j + k*N];
			}

		}

	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			s[i*N + j] = 0;
			for (int k = 0; k < N; k++)
			{
				s[i*N + j] += r_new[i*N + k] * v[k + j*N];
			//	s[i*N + j] += r_new[i*N + k] * v[j+k*N];
			}

		}

	}



	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			sigmafile.precision(4);
			sigmafile << s[i*N + j] << "\t";

		}
		sigmafile << endl;

	}




	delete[] r_new;

}


/*
int svd_hestenes(int m_m, int n_n, double * a, double * u, double * v, double * sigma)
{
	double thr = 1.E-4f, nul = 1.E-16f;
	int n, m, i, j, l, k, lort, iter, in, ll, kk;
	double alfa, betta, hamma, eta, t, cos0, sin0, buf, s;
	n = n_n;
	m = m_m;
	for (i = 0; i<n; i++)
	{
		in = i*n;
		for (j = 0; j<n; j++)
		if (i == j) v[in + j] = 1.;
		else v[in + j] = 0.;
	}
	for (i = 0; i<m; i++)
	{
		in = i*n;
		for (j = 0; j<n; j++)
		{
			u[in + j] = a[in + j];
		}
	}

	iter = 0;
	while (1)
	{
		lort = 0;
		iter++;
		for (l = 0; l<n - 1; l++)
		for (k = l + 1; k<n; k++)
		{
			alfa = 0.; betta = 0.; hamma = 0.;
			for (i = 0; i<m; i++)
			{
				in = i*n;
				ll = in + l;
				kk = in + k;
				alfa += u[ll] * u[ll];
				betta += u[kk] * u[kk];
				hamma += u[ll] * u[kk];
			}

			if (sqrt(alfa*betta) < nul)	continue;
			if (fabs(hamma) / sqrt(alfa*betta)<thr) continue;

			lort = 1;
			eta = (betta - alfa) / (2.f*hamma);
			t = (double)((eta / fabs(eta)) / (fabs(eta) + sqrt(1. + eta*eta)));
			cos0 = (double)(1. / sqrt(1. + t*t));
			sin0 = t*cos0;

			for (i = 0; i<m; i++)
			{
				in = i*n;
				buf = u[in + l] * cos0 - u[in + k] * sin0;
				u[in + k] = u[in + l] * sin0 + u[in + k] * cos0;
				u[in + l] = buf;

				if (i >= n) continue;
				buf = v[in + l] * cos0 - v[in + k] * sin0;
				v[in + k] = v[in + l] * sin0 + v[in + k] * cos0;
				v[in + l] = buf;
			}
		}

		if (!lort) break;
	}

	for (i = 0; i<n; i++)
	{
		s = 0.;
		for (j = 0; j<m; j++)	s += u[j*n + i] * u[j*n + i];
		s = (double)sqrt(s);
		sigma[i] = s;
		if (s < nul)	continue;
		for (j = 0; j<m; j++)	u[j*n + i] /= s;
	}
	//======= Sortirovka ==============
	for (i = 0; i<n - 1; i++)
	for (j = i; j<n; j++)
	if (sigma[i]<sigma[j])
	{
		s = sigma[i]; sigma[i] = sigma[j]; sigma[j] = s;
		for (k = 0; k<m; k++)
		{
			s = u[i + k*n]; u[i + k*n] = u[j + k*n]; u[j + k*n] = s;
		}
		for (k = 0; k<n; k++)
		{
			s = v[i + k*n]; v[i + k*n] = v[j + k*n]; v[j + k*n] = s;
		}
	}


	fstream f;
	f.open("svd.txt", ios::trunc | ios::out);
	f.precision(5);
	

	for (int i = 0; i < n; i++)
	{

		for (int j = 0; j < n; j++)
		{
			f  << u[i*n + j] << "\t";
		}
		f << endl;
	}

	f << endl;	f << endl;
	for (int i = 0; i < n; i++)
	{

		for (int j = 0; j < n; j++)
		{
			f  << v[i*n + j] << "\t";
		}
		f << endl;
	}





	return iter;
}
*/
int svd_double(int m_m, int n_n, double * a, double * u, double * v, double * sigma)
{
	fstream fileObj;
	fileObj.open("r_s_in_svd.txt", ios::out | ios::trunc);
	for (int i = 0; i < m_m; i++)
	{
		for (int j = 0; j < m_m; j++)
		{
		
			fileObj << a[i*m_m + j] << "\t";

		}
		fileObj << endl;

	}



	double thr = 0.000001f;
	int n, m, i, j, l, k, lort, iter, in, ll, kk;
	double alfa, betta, hamma, eta, t, cos0, sin0, buf, s, r;
	n = n_n;
	m = m_m;
	for (i = 0; i<n; i++)
	{
		in = i*n;
		for (j = 0; j<n; j++)
		if (i == j) v[in + j] = 1.;
		else v[in + j] = 0.;
	}
	for (i = 0; i<m; i++)
	{
		in = i*n;
		for (j = 0; j<n; j++)
		{
			u[in + j] = a[in + j];
		}
	}

	iter = 0;
	while (1)
	{
		lort = 0;
		iter++;
		for (l = 0; l<n - 1; l++)
		for (k = l + 1; k<n; k++)
		{
			alfa = 0.; betta = 0.; hamma = 0.;
			for (i = 0; i<m; i++)
			{
				in = i*n;
				ll = in + l;
				kk = in + k;
				alfa += u[ll] * u[ll];
				betta += u[kk] * u[kk];
				hamma += u[ll] * u[kk];
			}

			if (sqrt(alfa*betta) < 1.e-10)	continue;
			if (fabs(hamma) / sqrt(alfa*betta)<thr)
				continue;

			lort = 1;
			eta = (betta - alfa) / (2.f*hamma);
			t = (eta / fabs(eta)) /
				(fabs(eta) + (double)sqrt(1.f + eta*eta));
			cos0 = 1.f / (double)sqrt(1.f + t*t);
			sin0 = t*cos0;

			for (i = 0; i<m; i++)
			{
				in = i*n;
				buf = u[in + l] * cos0 - u[in + k] * sin0;
				u[in + k] = u[in + l] * sin0 + u[in + k] * cos0;
				u[in + l] = buf;

				if (i >= n) continue;
				buf = v[in + l] * cos0 - v[in + k] * sin0;
				v[in + k] = v[in + l] * sin0 + v[in + k] * cos0;
				v[in + l] = buf;
			}
		}

		if (!lort) break;
	}

	for (i = 0; i<n; i++)
	{
		s = 0.;
		for (j = 0; j<m; j++)	s += u[j*n + i] * u[j*n + i];
		s = (double)sqrt(s);
		sigma[i] = s;
		if (s < 1.e-10)	continue;
		for (j = 0; j<m; j++)	u[j*n + i] = u[j*n + i] / s;
	}
	for (i = 0; i<n - 1; i++)
	for (j = i; j<n; j++)
	if (sigma[i]<sigma[j])
	{
		r = sigma[i]; sigma[i] = sigma[j]; sigma[j] = r;
		for (k = 0; k<m; k++)
		{
			r = u[i + k*n]; u[i + k*n] = u[j + k*n]; u[j + k*n] = r;
		}
		for (k = 0; k<n; k++)
		{
			r = v[i + k*n]; v[i + k*n] = v[j + k*n]; v[j + k*n] = r;
		}
	}

	return iter;
}

void Filter(double* sigma, double* u, double* v, double* r_s, double thrhold, int L)
{


	for (int i = 1; i <= L; i++)
	{

		if (sigma[i%L] < sigma[0] * thrhold / 100.) sigma[i%L] = 0;

	}

	double* r_sbuf = new double[L*L];

	for (int i = 0; i < L; i++)
	{
		for (int j = 0; j < L; j++)
		{
			r_sbuf[i*L + j] = 0;
			for (int k = 0; k < L; k++)
			{
			//	r_sbuf[i*L + j] += u[i*L + k] * sigma[j + k*L];

				if(k==j)r_sbuf[i*L + j] += u[i*L + k] * sigma[k];
			
			}

		}

	}

	for (int i = 0; i < L; i++)
	{
		for (int j = 0; j < L; j++)
		{
			r_s[i*L + j] = 0;
			for (int k = 0; k < L; k++)
			{
				r_s[i*L + j] += r_sbuf[i*L + k] * v[k + j*L];
				//	s[i*N + j] += r_new[i*N + k] * v[j+k*N];
			}

		}

	}

	fstream sigmafile;
	sigmafile.open("filtered.txt", ios::out | ios::trunc);

	for (int i = 0; i < L; i++)
	{
		for (int j = 0; j < L; j++)
		{
			sigmafile.precision(4);
		//	sigmafile << v[i*L + j] << "\t";

		}
	//	sigmafile << endl;

	}

	delete[] r_sbuf;
	


 }