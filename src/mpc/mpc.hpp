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

#ifndef __MPC_HPP__
#define __MPC_HPP__

#include <iostream>

#include<Eigen/Core>
#include<Eigen/QR>
#include<Eigen/Dense>

#include "bfgs.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

#define PI 3.1415926

struct MPC_Setting
{
	double Lf;	// COG to front wheel
	double Lr;  	// COG to rear wheel

	double dt;	// sample time
	size_t N;	// horizon

	size_t n; 	// dimension of state
	size_t m;	// dimension of input

	MPC_Setting()
	{
		Lf = 1.17;
		Lr = 1.77;

		dt = 1/100.;
		N = 25;

		n = 6;
		m = 2;
	}
};

struct MPC_Argument
{
	double x, y, psi, v, cte, epsi;	// initial state
	double ref_v;
	double latency;
	double coeffs[3];

	MPC_Argument()
	{
		x = y = psi = v = cte = epsi = 0.;
		ref_v = 40/3.6;
		latency = 1/1000. * 10;
	}
};

double getFitness(const MPC_Setting &setting, const MPC_Argument &arg, const VectorXd &input);
void getGradient(VectorXd &grad,const MPC_Setting &setting, const MPC_Argument &arg, VectorXd &input, const double fit);
void getHess(MatrixXd &hess,const MPC_Setting &setting, const MPC_Argument &arg, VectorXd &input, const double fit);

inline double polyval(const double *coeffs, double x)
{
	return coeffs[0] + coeffs[1] * x + coeffs[2] * x * x;
}

inline double pow2(double x)
{
	return x*x;
}

struct Motion_Model:public NonConstraintObj
{
	MPC_Argument m_arg;
	MPC_Setting m_setting;

	Motion_Model(){}

	double cost(VectorXd const &x)
	{
		return getFitness(m_setting, m_arg, x);
	}
	void grad(VectorXd &diff, VectorXd &x)
	{
		double fit = cost(x);
		getGradient(diff, m_setting, m_arg, x, fit);
	}
};

#endif
