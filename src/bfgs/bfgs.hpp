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

#ifndef __BFGS_HPP__
#define __BFGS_HPP__

#include <functional>
#include <memory>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Dense>

using namespace Eigen;

using Eigen::MatrixXd;
using Eigen::VectorXd;

// ref: http://en.cppreference.com/w/cpp/utility/functional/function

// input is the position(parameter)
// Add const to avoid error of invalid initialization of non-const reference
typedef std::function<double(VectorXd const &)> Cost_Fun;
// param1: returned gradient
// param2: input position
typedef std::function<void(VectorXd &, VectorXd &)> Diff_Fun;

struct NonConstraintObj
{
	VectorXd m_x;
	size_t m_n;

	NonConstraintObj()
	{
		m_n = 1;
		m_x.resize(m_n);
		m_x.fill(0.);
	}
	virtual ~NonConstraintObj() = default;
	virtual double cost(VectorXd const &x) = 0;
	virtual void grad(VectorXd &grad, VectorXd &x) = 0;
};

struct BFGS
{
	size_t MAX_STEP;
	size_t stop_step;
	size_t step;
	double rho;
	double sigma;
	double epsilon;

	bool improved;
	bool initialized;

	Cost_Fun cost_fun;
	Diff_Fun grad_fun;

	std::shared_ptr<NonConstraintObj> funPtr;
	BFGS()
	{
		MAX_STEP = 500;
		stop_step = 5;
		step = 0;
		rho = 0.55;
		sigma = 0.4;
		epsilon = 1e-5;
		improved = false;
		initialized = false;
		funPtr = nullptr;
	}
	template <class T>
	void init(T &obj)
	{
		// store a call to a member function and object
		using std::placeholders::_1;
		using std::placeholders::_2;
		cost_fun = std::bind(&T::cost, obj, _1);
		grad_fun = std::bind(&T::grad, obj, _1, _2);
		initialized = true;
	}
	double solve_(Cost_Fun fun, Diff_Fun dfun, VectorXd &input);
	double solve(VectorXd &input)
	{
		if (initialized)
		{
			return solve_(cost_fun, grad_fun, input);
		}
		else
		{
			printf("BFGS is not initialized! Return zeros.\n");
			return 1e9;
		}
	}
};

struct BFGS_V2
{
	//double *Bk;

	BFGS_V2()
	{
		//Bk = nullptr;
		//cudaMalloc((void **)&Bk, sizeof(double)*3);
	}
	template <class T>
	static double bfgs(std::shared_ptr<T> ptr, VectorXd &x0)
	{
		size_t maxk = 500, stop_step = 5;
		double rho = 0.55, sigma = 0.4, epsilon = 1e-5;

		size_t k = 0, n = x0.size();
		MatrixXd Bk = MatrixXd::Identity(n, n);

		VectorXd gk, dk, sk, yk, x;
		size_t m = 0, mk = 0, bad = 0;
		double fit = 0., tmp = 0.;
		double best = 1e10;
		while (k < maxk && bad < stop_step)
		{
			fit = ptr->cost(x0);
			if (fit < best - epsilon)
			{
				best = fit;
				bad = 0;
			}
			else
			{
				bad++;
			}
			ptr->grad(gk, x0);
			//std::cout<<"gk = "<<gk<<std::endl;
			if (gk.norm() < epsilon)
			{
				break;
			}
			dk = Bk.householderQr().solve(-gk); // llt ?
			m = 0, mk = 0;

			tmp = gk.transpose() * dk;
			while (m < 100)
			{
				if (ptr->cost(x0 + pow(rho, m) * dk) < (fit + sigma * pow(rho, m) * tmp))
				{
					mk = m;
					break;
				}
				m++;
			}
			x = x0 + pow(rho, mk) * dk;
			sk = pow(rho, mk) * dk;
			ptr->grad(yk, x);
			yk = yk - gk;
			if (yk.transpose() * sk > 0)
			{
				Bk = Bk - (Bk * sk * sk.transpose() * Bk) / (sk.transpose() * Bk * sk) + (yk * yk.transpose()) / (yk.transpose() * sk);
			}
			x0 = x;
			k++;
		}
		x0 = x;
		return fit;
	}
};

#endif
