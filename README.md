# BFGS
Implementation of Broyden, Fletcher, Goldfarb and Shanno's (BFGS) quasi Newton's method in c++. 

## Installation

```sh
git clone git@github.com:gcjyzdd/BFGS.git
cd BFGS
chmod +x test.sh
./test.sh
```

## Examples

In the [main.cpp](./src/main.cpp), there is an example of applying BFGS to solve model predictive control.

To run the example, type

```sh
./bin/mpc_test
```

The result of MPC is shown below:

<div style="text-align:center"><img src ='./images/demo_mpc.png' /></div>

It takes `4.87416ms` to solve MPC per run on my PC.

## Dependencies

Eigen-3.3 or later.

## TODO

Show a demo of comparison with the matlab example [Swing-up Control of a Pendulum Using Nonlinear Model Predictive Control](https://nl.mathworks.com/help/mpc/ug/swing-up-control-of-a-pendulum-using-nonlinear-model-predictive-control.html).
