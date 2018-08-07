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
	virtual double cost(VectorXd const&x) = 0;
	// gradient of the cost function
	virtual void grad(VectorXd &grad, VectorXd &x) = 0;
	// constraint equalities
	virtual void hf(VectorXd &value, VectorXd const &x) = 0;
	// constraint inequalities
	virtual void gf(VectorXd &value, VectorXd const &x) = 0;
	// gradient of equalities
	virtual void dhf(MatrixXd &grad, VectorXd &x) = 0;
	// gradient of inequalities
	virtual void dgf(MatrixXd &grad, VectorXd &x) = 0;
	// mpsi
	virtual double mpsi(VectorXd const &x, VectorXd &mu, VectorXd &lambda, double sigma) = 0;
	// dmpsi
	virtual void dmpsi(VectorXd &grad, VectorXd &x, VectorXd &mu, VectorXd &lambda, double sigma) = 0;
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

	for (size_t i = 0;i < l;i++)
	{
		psi = psi - he[i] * mu[i];
		s1 = s1 + he[i] * he[i];
	}
	psi += 0.5*sigma*s1;
	double s2 = 0., s3 = 0.;
	for (size_t i = 0;i < m;i++)
	{
		s3 = std::max(0., lambda[i] - sigma * gi[i]);
		s2 += s3 * s3 - lambda[i] * lambda[i];
	}
	psi += s2 / (2. * sigma);
	return psi;
}

void dmpsi(VectorXd &dpsi, VectorXd &x,
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

	size_t l, m;
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

		fun = std::bind(&T::cost, obj, _1);
		dfun = std::bind(&T::grad, obj, _1, _2);

		hf = std::bind(&T::hf, obj, _1, _2);
		dhf = std::bind(&T::dhf, obj, _1, _2);

		gf = std::bind(&T::gf, obj, _1, _2);
		dgf = std::bind(&T::dgf, obj, _1, _2);

		initialized = true;
	}

	double solve(Cost_Fun fun, Constraint_Fun hf, Constraint_Fun gf, Diff_Fun dfun, Diff_Constraint_Fun dhf, Diff_Constraint_Fun dgf, VectorXd &x0);
};

/**
 * Solve constraint nonlinear optimization using multiplier method of PHR.
 * @param ptr A shared pointer to a function object. The object must define eight functions:
 * 			  fun, hf, gf, dfun, dhf, dgf, mpsi, dmpsi
 * @param x0 The input and output parameter
 * @return Cost of the returned parameter
 */
template<class T>
double multphr(std::shared_ptr<T> ptr, VectorXd &x0)
{
	// define parameters
	size_t maxk = 500, k = 0, ink = 0, n, l, m;
	double sigma = 2., eta = 2., theta = 0.8, epsilon = 1e-5;

	VectorXd x = x0, he, gi;

	// get dimension
	n = ptr->m_n;
	l = ptr->m_l;
	m = ptr->m_m;

	// initialize multipliers
	VectorXd mu(l);
	mu.fill(0.1);
	VectorXd lambda(m);
	lambda.fill(0.1);

	double btak = 10.;
	double btaold = 10.;

	Cost_Fun cost_fun;
	Diff_Fun diff_fun;

	using std::placeholders::_1;
	using std::placeholders::_2;

	BFGS bfgs_solver;
	double tmp;
	while (k<maxk && btak>epsilon)
	{
		cost_fun = [ptr, &mu, &lambda, &sigma](VectorXd const x) { return ptr->mpsi(x, mu, lambda, sigma);};
		diff_fun = [ptr, &mu, &lambda, &sigma](VectorXd &grad, VectorXd &x) { ptr->dmpsi(grad, x, mu, lambda, sigma);};
		bfgs_solver.solve_(cost_fun, diff_fun, x0);

		x = x0;
		ptr->hf(he, x);
		ptr->gf(gi, x);
		btak = 0.;
		for (size_t i = 0;i < l;i++)
		{
			btak += he[i] * he[i];
		}
		for (size_t i = 0;i < m;i++)
		{
			tmp = std::min(gi[i], lambda[i] / sigma);
			btak += tmp * tmp;
		}
		btak = std::sqrt(btak);

		if (btak > epsilon)
		{
			if ((k >= 2) && (btak > (theta*btaold)))
			{
				sigma *= eta;
			}
			// update multiplier vector
			for (size_t i = 0;i < l;i++) { mu[i] -= sigma * he[i]; }
			for (size_t i = 0;i < m;i++) { lambda[i] = std::max(0., lambda[i] - sigma * gi[i]); }
		}
		btaold = btak;
		x0 = x;
		k++;
	}

	return ptr->cost(x0);
}

namespace SISS {
	// input is the position(parameter)
	// Add const to avoid error of invalid initialization of non-const reference
	typedef std::function<double(VectorXd const &, VectorXd &, VectorXd &, double)> MPsi;
	// param1: returned gradient
	// param2: input position
	typedef std::function<void(VectorXd &, VectorXd &, VectorXd &, VectorXd &, double)> DMPsi;
	// param1: returned hessian matrix
	// param2: input position
	typedef std::function<void(MatrixXd &, VectorXd &, VectorXd &, VectorXd &, double)> HMPsi;

	struct BFGS_Hess
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

		BFGS_Hess()
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

		template<class T>
		double solve_(std::shared_ptr<T> ptr, VectorXd &input, VectorXd &mu, VectorXd &lam, double sigma)
		{
			size_t n = input.size();

			VectorXd gk(n);
			VectorXd x0(n);
			x0 = input;
			double fit = ptr->mpsi(x0, mu, lam, sigma);

			MatrixXd Bk = MatrixXd::Identity(n, n);
			ptr->hmpsi(Bk, input, mu, lam, sigma);

			VectorXd x, dk, sk, yk, valid;

			int k = 0, bad = 0;
			double best = 1e10;
			while (k < MAX_STEP && bad < stop_step)
			{
				fit = ptr->mpsi(x0, mu, lam, sigma);
				if (fit < best - epsilon)
				{
					best = fit;
					bad = 0;
				}
				else
				{
					bad++;
				}
				//std::cout<<"k = "<<k<<" fit = "<< fit<<std::endl;
				ptr->dmpsi(gk, x0, mu, lam, sigma);
				if (gk.norm() < epsilon) { break; }
				// Solve linear equations. Ref: https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
				dk = Bk.householderQr().solve(-gk);	// llt ?
													// Armijo search
				int m = 0;
				int mk = 0;
				double tmp;
				tmp = gk.transpose()*dk;
				while (m < 100)
				{
					if (ptr->mpsi(x0 + pow(rho, m)*dk, mu, lam, sigma) < (fit + sigma * pow(rho, m)*tmp))
					{
						mk = m;
						break;
					}
					m++;
				}
				//std::cout << "m = " << m << std::endl;
				x = x0 + pow(rho, mk)*dk;
				sk = pow(rho, mk)*dk;
				//fit = fun(x);
				ptr->dmpsi(yk, x, mu, lam, sigma);
				yk = yk - gk;
				if (yk.transpose() * sk > 0)
				{
					Bk = Bk - (Bk*sk*sk.transpose()*Bk) / (sk.transpose()*Bk*sk) + (yk*yk.transpose()) / (yk.transpose()*sk);
				}

				x0 = x;
				k++;
			}
			input = x0;
			step = k;
			return fit;
		}
	};

	/**
	* Solve constraint nonlinear optimization using multiplier method of PHR.
	* @param ptr A shared pointer to a function object. The object must define eight functions:
	* 			  fun, hf, gf, dfun, dhf, dgf, mpsi, dmpsi
	* @param x0 The input and output parameter
	* @return Cost of the returned parameter
	*/
	template<class T>
	double multphr_hess(std::shared_ptr<T> ptr, VectorXd &x0)
	{
		// define parameters
		size_t maxk = 500, k = 0, ink = 0, n, l, m;
		double sigma = 2., eta = 2., theta = 0.8, epsilon = 1e-5;

		VectorXd x = x0, he, gi;

		// get dimension
		n = ptr->m_n;
		l = ptr->m_l;
		m = ptr->m_m;

		// initialize multipliers
		VectorXd mu(l);
		mu.fill(0.1);
		VectorXd lambda(m);
		lambda.fill(0.1);

		double btak = 10.;
		double btaold = 10.;

		Cost_Fun cost_fun;
		Diff_Fun diff_fun;

		using std::placeholders::_1;
		using std::placeholders::_2;

		BFGS_Hess bfgs_solver;
		double tmp;
		while (k<maxk && btak>epsilon)
		{
			//cost_fun = [ptr, &mu, &lambda, &sigma](VectorXd const x) { return ptr->mpsi(x, mu, lambda, sigma);};
			//diff_fun = [ptr, &mu, &lambda, &sigma](VectorXd &grad, VectorXd &x) { ptr->dmpsi(grad, x, mu, lambda, sigma);};
			bfgs_solver.solve_<T>(ptr, x0, mu, lambda, sigma);

			x = x0;
			ptr->hf(he, x);
			ptr->gf(gi, x);
			btak = 0.;
			for (size_t i = 0;i < l;i++)
			{
				btak += he[i] * he[i];
			}
			for (size_t i = 0;i < m;i++)
			{
				tmp = std::min(gi[i], lambda[i] / sigma);
				btak += tmp * tmp;
			}
			btak = std::sqrt(btak);

			if (btak > epsilon)
			{
				if ((k >= 2) && (btak > (theta*btaold)))
				{
					sigma *= eta;
				}
				// update multiplier vector
				for (size_t i = 0;i < l;i++) { mu[i] -= sigma * he[i]; }
				for (size_t i = 0;i < m;i++) { lambda[i] = std::max(0., lambda[i] - sigma * gi[i]); }
			}
			btaold = btak;
			x0 = x;
			k++;
		}

		return ptr->cost(x0);
	}
}
#endif
