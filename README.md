# pbetools

[![CI](https://github.com/HugoMVale/pbetools/workflows/CI/badge.svg)](https://github.com/HugoMVale/pbetools/actions)
[![codecov](https://codecov.io/gh/HugoMVale/pbetools/branch/main/graph/badge.svg?token=1XL5LQSO9P)](https://codecov.io/gh/HugoMVale/pbetools)

## Status

`pbetools` is currently being developed and not yet ready for use.

## Description

`pbetools` is a modern-Fortran package to solve population balance equations (PBE) for one- and two-component aggregation processes using the (extended) fixed pivot method. For single component systems, the code implements the method of [Kumar & Ramkrishna (1996)](https://doi.org/10.1016/0009-2509(96)88489-2), and for bivariate aggregation the method of [Vale & McKenna (2005)](https://doi.org/10.1021/ie050179s).

## Underlying PBE

If the system is spatially homogeneous, a two-component aggregation process is described by the following PBE:

<img src="https://latex.codecogs.com/svg.latex?\large&space;\begin{align*}&space;&&space;\frac{{\partial&space;n(x,y,t)}}{{\partial&space;t}}&space;=&space;\\&space;&&space;\frac{1}{2}\int_0^x&space;{\int_0^y&space;{\beta&space;(x&space;-&space;x',y&space;-&space;y';x',y';t)n(x&space;-&space;x',y&space;-&space;y',t)n(x',y',t)&space;d&space;x'&space;d&space;y'}}&space;\\&space;&&space;-n(x,y,t)\int_0^\infty&space;{\int_0^\infty&space;{\beta&space;(x,y;x',y';t)n(x',y',t)d&space;x'&space;d&space;y'}}&space;\end{align*}" title="\large \begin{align*} & \frac{{\partial n(x,y,t)}}{{\partial t}} = \\ & \frac{1}{2}\int_0^x {\int_0^y {\beta (x - x',y - y';x',y';t)n(x - x',y - y',t)n(x',y',t) d x' d y'}} \\ & -n(x,y,t)\int_0^\infty {\int_0^\infty {\beta (x,y;x',y';t)n(x',y',t)d x' d y'}} \end{align*}" />

where <img src="https://latex.codecogs.com/svg.latex?\inline&space;n(x,y,t)&space;dx&space;dy" title="n(x,y,t) dx dy" /> is the number of particles of state <img src="https://latex.codecogs.com/svg.latex?\inline&space;(x,y)" title="(x,y)" /> per unit volume at time <img src="https://latex.codecogs.com/svg.latex?t" title="t" /> and <img src="https://latex.codecogs.com/svg.latex?\beta&space;(x,y;x',y';t)" title="\beta (x,y;x',y';t)" /> is the aggregation rate coefficient. The internal coordinates <img src="https://latex.codecogs.com/svg.latex?x" title="x" /> and <img src="https://latex.codecogs.com/svg.latex?y" title="y" /> denote the amount (mass, moles, etc.) of each component in the particle.

## Getting started

### Build

The easiest way to build/test the code and run the examples is by means of [`fpm`](https://fpm.fortran-lang.org/en/index.html). To run a given example, just do:

```
fpm run --example "example-filename"
```

and the numerical results will be stored in the [`output`](/output) subfolder. You can then use the provided Python script to read the data and plot the results.

### Usage

comming soon...

## Examples

comming soon...
