
SE3 with Gaussian Distribution
Emanuele Ruffaldi, Scuola Superiore Sant'Anna 2014-2015

Overview
---------------------------------------------

Function organization for use in Simulink

Representation: 
	group = rotation and translation as full 4x4 matrix
	adjoint = adjoint 6x6 representation of the group for transforming parameters
	algebra as parametric = 6 values of angular and translational velocity
	gaussian group = group with gaussian as 6x6 covariance

Function Organization
---------------------------------------------
Useful for Matlab Simulink code generation, all under se3fx

- representation:
	group: matrix 4x4
	alg:   vector 1x6
	dist:  flattened(16+36): matrix4x4 followed by covariance 6x6)
- operations on group (se3_): conversion, mul, invert, exp, log, adj
- operations on group distribution (se3d_): conversion/build, mul, invert, fuse, sample, build from data

Class Organization
------------------------

Class se3 and se3d wrapping the above functions


TODO/Recover
------------

Function transformation (Sigma Point) over the se3d with callback and without

