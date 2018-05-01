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

#ifndef __PHR_HPP__
#define __PHR_HPP__

#include "bfgs.hpp"

// param1: returned vector
// param2: input position
typedef std::function<void(VectorXd &, VectorXd const &)> Constraint_Fun;
typedef std::function<void(MatrixXd &, VectorXd &)> Diff_Constraint_Fun;

struct ConstraintObj
{
	VectorXd m_x;
	size_t m_n;	// dimension of goal pos
	size_t m_l;	// dimension of constraint equalities
	size_t m_m; // dimension of constraint inequalities

	ConstraintObj()
	{
		m_n = 1;
		m_l = 1;
		m_m = 1;
		m_x.resize(m_n);
		m_x.fill(0.);
	}
	virtual ~ConstraintObj() = default;
	// cost function of the goal
	virtual double cost(VectorXd const&x) =0;
	// gradient of the cost function
	virtual void grad(VectorXd &grad, VectorXd &x)=0;
	// constraint equalities
	virtual void hf(VectorXd &value, VectorXd const &x) = 0;
	// constraint inequalities
	virtual void gf(VectorXd &value, VectorXd const &x) = 0;
	// gradient of equalities
	virtual void dhf(MatrixXd &grad, VectorXd &x) = 0;
	// gradient of inequalities
	virtual void dgf(MatrixXd &grad, VectorXd &x) = 0;
	// mpsi
	virtual double mpsi();
	// dmpsi
	virtual void dmpsi();
};

double mpsi(VectorXd const &x,
		Cost_Fun fun, Constraint_Fun hf, Constraint_Fun gf,
		Diff_Fun dfun, Diff_Constraint_Fun dhf, Diff_Constraint_Fun dgf,
		VectorXd &mu, VectorXd lambda, double sigma)
{
	double f = fun(x);
	VectorXd he, gi;
	hf(he, x);
	gf(gi, x);

	size_t l = he.size();
	size_t m = gi.size();
	double s1 = 0.;
	double psi = f;

	for(size_t i=0;i<l;i++)
	{
		psi = psi - he[i] * mu[i];
		s1 = s1 + he[i]*he[i];
	}
	psi += 0.5*sigma*s1;
	double s2 = 0., s3=0.;
	for(size_t i=0;i<m;i++)
	{
		s3 = std::max(0., lambda[i] - sigma*gi[i]);
		s2 += s3*s3 - lambda[i]*lambda[i];
	}
	psi += s2/(2. * sigma);
	return psi;
}

void dmpsi(VectorXd &dpsi, VectorXd const &x,
		Cost_Fun fun, Constraint_Fun hf, Constraint_Fun gf,
		Diff_Fun dfun, Diff_Constraint_Fun dhf, Diff_Constraint_Fun dgf,
		VectorXd &mu, VectorXd lambda, double sigma)
{
	dfun(dpsi, x);
	VectorXd he, gi;
	MatrixXd dhe, dgi;

	hf(he, x);
	gf(gi, x);
	dhf(dhe, x);
	dgf(dgi, x);

	size_t l,m;
	l = he.size();
	m = gi.size();

}

struct PHR
{
	size_t maxk;
	double sigma, eta, theta, epsilon;
	bool initialized;

	Cost_Fun fun;
	Constraint_Fun hf;
	Constraint_Fun gf;
	Diff_Fun dfun;
	Diff_Constraint_Fun dhf;
	Diff_Constraint_Fun dgf;

	PHR()
	{
		maxk = 500;
		sigma = 2.0;
		eta = 2.0;
		theta = 0.8;
		epsilon = 1e-5;
		initialized = false;
	}
	template<class T> void init(T &obj)
	{
		// store a call to a member function and object
		using std::placeholders::_1;
		using std::placeholders::_2;

		fun = std::bind( &T::cost, obj, _1 );
		dfun = std::bind( &T::grad, obj, _1,_2 );

		hf = std::bind( &T::hf, obj, _1,_2 );
		dhf = std::bind( &T::dhf, obj, _1,_2 );

		gf = std::bind( &T::gf, obj, _1,_2 );
		dgf = std::bind( &T::dgf, obj, _1,_2 );

		initialized = true;
	}

	double solve(Cost_Fun fun, Cost_Fun hf, Cost_Fun gf, Diff_Fun dfun, Diff_Constraint_Fun dhf, Diff_Constraint_Fun dgf, VectorXd &x0);
};

#endif
