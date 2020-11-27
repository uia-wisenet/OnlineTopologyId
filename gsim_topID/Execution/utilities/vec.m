function [y] = vec(X)
%
% [Xv] = vec(X)
%
% Vectorizing (vec) operator.
%
% This operator stacks the columns of the matrix X into a 
% single column vector y.

[a,b] = size(X);

y = reshape(X,a*b,1);