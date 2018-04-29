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

#include<functional>
#include<iostream>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/Dense"

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


struct BFGS
{
	size_t MAX_STEP;
	size_t stop_step;
	double rho;
	double sigma;
	double epsilon;

	bool improved;

	BFGS()
	{
		MAX_STEP = 500;
		stop_step = 5;
		rho = 0.55;
		sigma = 0.4;
		epsilon = 1e-5;
		improved = false;
	}

	void solve( Cost_Fun fun, Diff_Fun dfun, VectorXd &input);
};


#endif
