function obj = assignParametersByName(obj,varargin)
%
% this function is useful for constructors
%
% Input:
% OBJ       Object of an arbitrary class
% VARARGIN  Sequence of (parameter,value) pairs. The first
%           element of each pair is the name of a parameter of
%           OBJ.
%
% Output:
% OBJ       Result of setting the parameters of input OBJ
%           according to VARARGIN
%
% Example:
%     car = Parameter.assignParametersByName( car ,'numberOfWheels',4,'year',2009,'color','blue')
%
if nargin>1
	
	if mod(length(varargin),2)~=0
		error('the number of arguments must be odd');
	end
	
	for k=1:2:length(varargin)-1
		obj.(varargin{k})=varargin{k+1};
	end
	
end
end