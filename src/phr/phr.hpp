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

struct PHR
{
	size_t maxk;
	double sigma, eta, theta, epsilon;

	PHR()
	{
		maxk = 500;
		sigma = 2.0;
		eta = 2.0;
		theta = 0.8;
		epsilon = 1e-5;
	}
	double solve(Cost_Fun fun, Cost_Fun hf, Cost_Fun gf, Diff_Fun dfun, Diff_Constraint_Fun dhf, Diff_Constraint_Fun dgf, VectorXd &x0);
};

#endif
