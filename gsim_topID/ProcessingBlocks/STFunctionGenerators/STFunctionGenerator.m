classdef STFunctionGenerator<ProcessingBlock
    % GENERATION OF SPACE-TIME SERIES
    
	properties
		nTimeSamples = []; % when empty, return all data in the data set
		% (real data generator) or raise an error otherwise.
	end
	
	methods
		
		function obj = STFunctionGenerator(varargin)
			obj@ProcessingBlock(varargin{:});
		end
		
		
	end
	
	methods(Abstract)
		
		N   = numberOfNodes(obj)
		% returns the number of time series of the data set,
		% which should equal the first dimension of the matrix
		% returned by the method realization.
		
		m_X = realization(obj,nAdditionalSamples)
			%
            % Generates a realization of a space-time function 
			%
			% m_X is an M x T matrix whose (m,t)-th entry is the t-th
			% sample of the m-th time series. 
			%
			% nAdditionalSamples (optional) is an integer (default
			% nAdditionalSamples = 0). This parameter is useful for
			% estimating the prediction error: the first obj.nTimeSamples
			% columns of m_X can be used as training samples, and the last
			% nAdditionalSamples as test samples.
			%
			%
		
	end
	
end

