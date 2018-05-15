//
// Author: Changjie Guan (changjieguan@gmail.com)
//
//
//BSD 2-Clause License
//
//Copyright (c) 2018, Changjie
//All rights reserved.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//
//* Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
//* Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
//OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include <iostream>
#include <chrono>
#include <vector>
#include <functional>
#include <math.h>

#include "phr.hpp"
#include "bfgs.hpp"

using namespace std;

struct TestFun : public ConstraintObj
{
	double d;
	TestFun()
	{
		d = 1e-8;
		m_n = 2;
		m_l = 1;
		m_m = 1;
		m_x.resize(m_n);
		m_x.fill(0.);
	}

	double cost(VectorXd const &x)
	{
		return pow(x[0] - 2., 2) + pow(x[1] - 1., 2);
	}

	void grad(VectorXd &diff, VectorXd &x)
	{
		diff.resize(x.size());
		diff.fill(0.);
		double fit = cost(x);
		for (size_t i = 0; i < x.size(); i++)
		{
			x[i] = x[i] + d;
			diff[i] = (cost(x) - fit) / d;
			x[i] = x[i] - d;
		}
	}

	void hf(VectorXd &value, VectorXd const &x)
	{
		value.resize(m_l);
		value[0] = x[0] - 2. * x[1] + 1.;
	}

	void gf(VectorXd &value, VectorXd const &x)
	{
		value.resize(m_m);
		value[0] = -0.25 * x[0] * x[0] - x[1] * x[1] + 1.;
	}

	void dhf(MatrixXd &grad, VectorXd &x)
	{
		grad.resize(m_n, m_l);
		grad(0, 0) = 1.;
		grad(1, 0) = -2.;
	}

	void dgf(MatrixXd &grad, VectorXd &x)
	{
		grad.resize(m_n, m_m);
		grad(0, 0) = -0.5 * x[0];
		grad(1, 0) = -2. * x[1];
	}

	double mpsi(VectorXd const &x, VectorXd &mu, VectorXd &lambda, double sigma)
	{
		double f = cost(x);
		VectorXd he, gi;
		hf(he, x);
		gf(gi, x);

		size_t l = he.size();
		size_t m = gi.size();
		double s1 = 0.;
		double psi = f;

		for (size_t i = 0; i < l; i++)
		{
			psi = psi - he[i] * mu[i];
			s1 = s1 + he[i] * he[i];
		}
		psi += 0.5 * sigma * s1;
		double s2 = 0., s3 = 0.;
		for (size_t i = 0; i < m; i++)
		{
			s3 = std::max(0., lambda[i] - sigma * gi[i]);
			s2 += s3 * s3 - lambda[i] * lambda[i];
		}
		psi += s2 / (2. * sigma);
		return psi;
	}

	void dmpsi(VectorXd &dpsi, VectorXd &x, VectorXd &mu, VectorXd &lambda, double sigma)
	{
		grad(dpsi, x);
		VectorXd he, gi;
		MatrixXd dhe, dgi;
		hf(he, x);
		gf(gi, x);
		dhf(dhe, x);
		dgf(dgi, x);

		dpsi = dpsi + (sigma * he[0] - mu[0]) * dhe.col(0);

		dpsi = dpsi + (sigma * gi[0] - lambda[0]) * dgi.col(0);
	}
};

int main()
{

	TestFun f1;
	int iters = 50;
	double cost;
	auto begin = std::chrono::steady_clock::now();

	for (size_t i = 0; i < iters; i++)
	{
		f1.m_x[0] = 3.0;
		f1.m_x[1] = 3.0;
		cost = multphr<TestFun>(std::make_shared<TestFun>(f1), f1.m_x);
	}
	auto end = std::chrono::steady_clock::now();

	std::cout << "Average Time difference = "
			  << std::chrono::duration_cast<std::chrono::microseconds>(end -
																	   begin)
						 .count() /
					 1000. / (float)iters
			  << "ms \n";

	cout << "cost = " << cost << endl;
	cout << "x = " << f1.m_x.transpose() << endl;

	BFGS solver;

	begin = std::chrono::steady_clock::now();
	for (size_t i = 0; i < iters; i++)
	{
		f1.m_x[0] = 3.0;
		f1.m_x[1] = 3.0;
		cost = BFGS_V2::bfgs<TestFun>(std::make_shared<TestFun>(f1), f1.m_x);
	}
	end = std::chrono::steady_clock::now();

	std::cout << "Average Time difference = "
			  << std::chrono::duration_cast<std::chrono::microseconds>(end -
																	   begin)
						 .count() /
					 1000. / (float)iters
			  << "ms \n";

	cout << "cost = " << cost << endl;
	cout << "x = " << f1.m_x.transpose() << endl;

	return 0;
}
