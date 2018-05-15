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

#include "mpc.hpp"
#include "bfgs.hpp"

using namespace std;

int main()
{
	BFGS solver;
	Motion_Model mpc;

	mpc.m_arg.coeffs[0] = -1;
	mpc.m_arg.coeffs[1] = 0;
	mpc.m_arg.coeffs[2] = 0;

	mpc.m_arg.x = -1;
	mpc.m_arg.y = 10;
	mpc.m_arg.psi = 0;
	mpc.m_arg.v = 10;
	mpc.m_arg.cte = polyval(mpc.m_arg.coeffs, mpc.m_arg.x) - mpc.m_arg.y;
	mpc.m_arg.epsi = mpc.m_arg.psi - atan(mpc.m_arg.coeffs[1] + 2 * mpc.m_arg.coeffs[2]);

	mpc.m_setting.N = 25;
	mpc.m_setting.m = 2;
	mpc.m_setting.Lf = 1.17;
	mpc.m_setting.dt = 0.05;
	mpc.m_arg.ref_v = 40;

	solver.init<Motion_Model>(mpc);

	int iters = 50;
	double cost;

	std::vector<float> result(2 + mpc.m_setting.m * mpc.m_setting.N);
	VectorXd rst(mpc.m_setting.m * (mpc.m_setting.N - 1));

	// Test speed of two versions. They have similar performance.
	auto begin = std::chrono::steady_clock::now();
	for (size_t i = 0; i < iters; i++)
	{
		rst.fill(0.);
		cost = solver.solve(rst);
	}
	auto end = std::chrono::steady_clock::now();
	std::cout << "Bind average Time difference = "
			  << std::chrono::duration_cast<std::chrono::microseconds>(end -
																	   begin)
						 .count() /
					 1000. / (float)iters
			  << "ms \n";

	cout << "cost = " << cost << ", step = " << solver.step << endl;
	BFGS_V2 bfgs2;
	begin = std::chrono::steady_clock::now();

	for (size_t i = 0; i < iters; i++)
	{
		rst.fill(0.);
		cost = bfgs2.bfgs<Motion_Model>(std::make_shared<Motion_Model>(mpc), rst);
	}
	end = std::chrono::steady_clock::now();
	std::cout << "Template average Time difference = "
			  << std::chrono::duration_cast<std::chrono::microseconds>(end -
																	   begin)
						 .count() /
					 1000. / (float)iters
			  << "ms \n";

	cout << "cost = " << cost << endl;
	return 0;
}
