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
#include<iostream>
#include <chrono>
#include<vector>
#include<functional>
#include <math.h>

#include "mpc.hpp"
#include "bfgs.hpp"

using namespace std;


struct Fun1: public NonConstraintObj
{
	double d;
	double t_;
	Fun1()
	{
		m_x.resize(2);
		m_x[0] = 0.;
		m_x[1] = 0.;
		d = 1e-8;
		t_ = 2.;
	}

	double cost(VectorXd const &x)
	{
		return 100*pow(x[0]*x[0] - x[1],2) + pow(x[0] - 1.,2) + t_;
	}

	void grad(VectorXd &diff, VectorXd &x)
	{
		diff.resize(x.size());
		diff.fill(0.);
		double fit = cost(x);
		for(size_t i=0;i<x.size();i++)
		{
			x[i] = x[i]+d;
			diff[i] = (cost(x) - fit)/d;
			x[i] = x[i]-d;
		}
	}
};


int main()
{
	BFGS solver;
	//NonConstraintObj *ptr = new Fun1(2);
	Fun1 f1;
	cout<<"ptr->cost(ptr->m_x) = "<< f1.cost(f1.m_x)<<endl;

	solver.init<Fun1>(f1);
	cout<<"b "<< solver.cost_fun(f1.m_x)<<endl;
	f1.t_ = 3.;
	solver.init<Fun1>(f1);
	cout<<"a "<< solver.cost_fun(f1.m_x)<<endl;

	int iters = 50;
	double cost;

	// test the template function bfgs
	f1.m_x[0] = 0;
	f1.m_x[1] = 0;
	cost = bfgs<Fun1>(std::make_shared<Fun1>(f1), f1.m_x);
	cout<<"cost = "<<cost<<endl;
	cout<<"x = "<<f1.m_x.transpose()<<endl;


	std::chrono::steady_clock::time_point begin =
			std::chrono::steady_clock::now();

	for (size_t i = 0; i < iters; i++) {
		f1.m_x[0] = 0;
		f1.m_x[1] = 0;
		cost = solver.solve(f1.m_x);
		//cost = bfgs<Fun1>(std::make_shared<Fun1>(f1), f1.m_x);
	}

	std::chrono::steady_clock::time_point end =
			std::chrono::steady_clock::now();
	std::cout << "Average Time difference = "
			<< std::chrono::duration_cast<std::chrono::microseconds>(end -
					begin).count() / 1000./(float)iters
					<< "ms \n";

	cout<<"cost = "<<cost<<", step = "<<solver.step<<endl;
	cout<<"x = "<<f1.m_x.transpose()<<endl;
	return 0;
}
