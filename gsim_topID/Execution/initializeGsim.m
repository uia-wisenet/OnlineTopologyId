classdef initializeGsim
	% This class contains methods for initializing the simulator
    
	properties
	end
	
	methods(Static)
		
        % Initialize MATLAB path to include all relevant classes. 
		function initializePath 
			% DEFAULT PATH
			addpath(genpath(['./Experiments']));
            addpath(genpath(['./ProcessingBlocks/OnlineGLEstimators']))                     
            addpath(genpath(['./ProcessingBlocks/STFunctionGenerators']))
			addpath(genpath(['./ProcessingBlocks/Estimators']))
 			addpath(genpath(['./ProcessingBlocks']));
 			addpath(genpath(['./Simulators']));
			
		end
	end
	
end

