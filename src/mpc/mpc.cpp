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
#include "mpc.hpp"
#include <math.h>

double getFitness(const MPC_Setting &setting, const MPC_Argument &arg, const VectorXd &input)
{
	VectorXd state0;
	VectorXd state1;

	double fit = 0;

	state0.resize(6);
	state0.fill(0.);
	state1.resize(6);
	state1.fill(0.);

	state0[0] = arg.x;
	state0[1] = arg.y;
	state0[2] = arg.psi;
	state0[3] = arg.v;
	state0[4] = arg.cte;
	state0[5] = arg.epsi;

	double Lf = setting.Lf;
	double dt = setting.dt;
	double ref_v = arg.ref_v;
	size_t N = setting.N;

	fit += pow2(state0[3] - ref_v);
	fit += pow2(state0[4]);
	fit += pow2(state0[5]);

	for (size_t i = 0; i < (N - 1); i++)
	{
		state1[0] = state0[0] + state0[3] * cos(state0[2]) * dt;
		state1[1] = state0[1] + state0[3] * sin(state0[2]) * dt;
		state1[2] = state0[2] + state0[3] * input[i] / Lf * dt;
		state1[3] = state0[3] + input[N - 1 + i] * dt;
		state1[4] = polyval(arg.coeffs, state0[0]) - state0[1] + state0[3] * sin(state0[5]) * dt;
		state1[5] = state0[2] - atan(arg.coeffs[1] + 2 * arg.coeffs[2]) + state0[3] * input[i] / Lf * dt;

		fit += pow2(state1[3] - ref_v);
		fit += pow2(state1[4]);
		fit += pow2(state1[5]);

		fit += pow2(input[i]);
		fit += pow2(input[N - 1 + i]);

		state0[0] = state1[0];
		state0[1] = state1[1];
		state0[2] = state1[2];
		state0[3] = state1[3];
		state0[4] = state1[4];
		state0[5] = state1[5];
	}

	for (size_t i = 0; i < (N - 2); i++)
	{
		fit += pow2(input[i] - input[i + 1]);
		fit += pow2(input[N - 1 + i] - input[N - 1 + i + 1]);
	}
	return fit;
}

/**
 * Get predicted states.
 *
 * @param _result %the size of result is 2 + 6 + 6*N, where N is the prediction
 * horizon.
 * */
void getPrediction(const MPC_Setting &setting, const MPC_Argument &arg,
				   const VectorXd &input, std::vector<float> &result)
{

	// TODO: make each element continous such that it can copy memory easily

	// Only store x and y. Other states are not required to calculate x and y.
	if (result.size() < (2 + 2 * setting.N))
	{
		std::cerr << "The size of result is not initialized correctly.\n";
	}

	result[0] = input[1]; // acc
	result[1] = input[0]; // steer

	size_t N = setting.N;
	size_t x_start = 2;
	size_t y_start = x_start + N;

	VectorXd state0;
	VectorXd state1;

	state0.resize(6);
	state0.fill(0.);
	state1.resize(6);
	state1.fill(0.);

	state0[0] = arg.x;
	state0[1] = arg.y;
	state0[2] = arg.psi;
	state0[3] = arg.v;

	result[x_start] = state0[0];
	result[y_start] = state0[1];

	double Lf = setting.Lf;
	double dt = setting.dt;
	double ref_v = arg.ref_v;

	for (size_t i = 0; i < (N - 1); i++)
	{
		state1[0] = state0[0] + state0[3] * cos(state0[2]) * dt;
		state1[1] = state0[1] + state0[3] * sin(state0[2]) * dt;
		state1[2] = state0[2] + state0[3] * input[i] / Lf * dt;
		state1[3] = state0[3] + input[N - 1 + i] * dt;

		state0[0] = state1[0];
		state0[1] = state1[1];
		state0[2] = state1[2];
		state0[3] = state1[3];

		result[x_start + i + 1] = state0[0];
		result[y_start + i + 1] = state0[1];
	}
}

void getGradient(VectorXd &grad, const MPC_Setting &setting, const MPC_Argument &arg, VectorXd &input, const double fit)
{

	double d = 1e-6;
	size_t N = (setting.N - 1) * setting.m;
	grad.resize(N);
	grad.fill(0.);

	for (size_t i = 0; i < N; i++)
	{
		input[i] = input[i] + d;
		grad[i] = (getFitness(setting, arg, input) - fit) / d;
		input[i] = input[i] - d;
	}
}

void getHess(MatrixXd &hess, const MPC_Setting &setting, const MPC_Argument &arg, VectorXd &input, const double fit)
{
	size_t N = (setting.N - 1) * setting.m;
	hess.resize(N, N);
	double d = 1e-9;
	double h2, h1, h12;

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j <= i; j++)
		{
			input[i] += d;
			h1 = getFitness(setting, arg, input);
			input[j] += d;
			h12 = getFitness(setting, arg, input);
			input[i] -= d;
			h2 = getFitness(setting, arg, input);
			input[j] -= d;

			hess(i, j) = (h12 - h1 - h2 + fit) / d / d;
			hess(j, i) = hess(i, j);
		}
	}
}
