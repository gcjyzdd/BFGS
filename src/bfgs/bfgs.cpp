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

#include "bfgs.hpp"

double BFGS::solve_(Cost_Fun fun, Diff_Fun dfun, Hess_Fun hfun, VectorXd &input)
{
	size_t n = input.size();

	VectorXd gk(n);
	VectorXd x0(n);
	x0 = input;
	double fit = fun(x0);

	MatrixXd Bk = MatrixXd::Identity(n, n);
	hfun(Bk, x0);

	VectorXd x, dk, sk, yk, valid;

	int k = 0, bad = 0;
	double best = 1e10;
	while (k<MAX_STEP && bad<stop_step)
	{
		fit = fun(x0);
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
		dfun(gk, x0);
		if (gk.norm() < epsilon) { break; }
		// Solve linear equations. Ref: https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
		dk = Bk.householderQr().solve(-gk);	// llt ?
											// Armijo search
		int m = 0;
		int mk = 0;
		double tmp;
		tmp = gk.transpose()*dk;
		while (m<100)
		{
			if (fun(x0 + pow(rho, m)*dk) < (fit + sigma * pow(rho, m)*tmp))
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
		dfun(yk, x);
		yk = yk - gk;
		if (yk.transpose() * sk >0)
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

double BFGS::solve_( Cost_Fun fun, Diff_Fun dfun, VectorXd &input)
{
	size_t n = input.size();

	VectorXd gk(n);
	VectorXd x0(n);
	x0 = input;
	double fit = fun(x0);

	MatrixXd Bk = MatrixXd::Identity(n, n);
	VectorXd x, dk, sk, yk, valid;

	int k=0, bad = 0;
	double best = 1e10;
	while(k<MAX_STEP && bad<stop_step)
	{
		fit = fun(x0);
		if(fit < best - epsilon)
		{
			best = fit;
			bad = 0;
		}
		else
		{
			bad++;
		}
		//std::cout<<"k = "<<k<<" fit = "<< fit<<std::endl;
		dfun(gk, x0);
		if(gk.norm() < epsilon){break;}
		// Solve linear equations. Ref: https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
		dk = Bk.householderQr().solve(-gk);	// llt ?
		// Armijo search
		int m = 0;
		int mk = 0;
		double tmp;
		tmp = gk.transpose()*dk;
		while(m<100)
		{
			if( fun(x0+pow(rho,m)*dk) < (fit+sigma*pow(rho,m)*tmp) )
			{
				mk = m;
				break;
			}
			m++;
		}
		//std::cout << "m = " << m << std::endl;
		x = x0 + pow(rho,mk)*dk;
		sk = pow(rho, mk)*dk;
		//fit = fun(x);
		dfun(yk, x);
		yk = yk - gk;
		if(yk.transpose() * sk >0)
		{
			Bk = Bk - (Bk*sk*sk.transpose()*Bk)/(sk.transpose()*Bk*sk) + (yk*yk.transpose())/(yk.transpose()*sk);
		}

		x0 = x;
		k++;
	}
	input = x0;
	step = k;
	return fit;
}
