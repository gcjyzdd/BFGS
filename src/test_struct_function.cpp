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


struct Fun1
{
	VectorXd m_x;
	double d;
	Fun1()
	{
		m_x.resize(2);
		m_x[0] = 0.;
		m_x[1] = 0.;
		d = 1e-8;
	}
	double cost(VectorXd const &x)
	{
		return 100*pow(x[0]*x[0] - x[1],2) + pow(x[0] - 1.,2);
	}
	void diff(VectorXd &diff, VectorXd &x)
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
	Fun1 f1;

	// store a call to a member function and object
	using std::placeholders::_1;
	using std::placeholders::_2;
	Cost_Fun cost_fun1 = std::bind( &Fun1::cost, f1, _1 );
	//std::function<double(VectorXd &)> cost_fun1 = std::bind( &Fun1::cost, f1, _1 );
	Diff_Fun diff_fun1 = std::bind( &Fun1::diff, f1, _1,_2 );

	int iters = 50;
	double cost;
	std::chrono::steady_clock::time_point begin =
			std::chrono::steady_clock::now();

	for (size_t i = 0; i < iters; i++) {
		f1.m_x[0] = 0;
		f1.m_x[1] = 0;
		cost = solver.solve(cost_fun1,diff_fun1,f1.m_x);
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
