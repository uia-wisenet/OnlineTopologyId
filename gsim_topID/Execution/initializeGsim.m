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
            addpath(genpath(['./ProcessingBlocks/graphGenerators']))
            addpath(genpath(['./ProcessingBlocks/STFunctionGenerators']))
            addpath(genpath(['./ProcessingBlocks/StructuralBreakDetectors']))
			addpath(genpath(['./ProcessingBlocks/STFunctionEstimators']))
			addpath(genpath(['./ProcessingBlocks/TimeVaryingVARBlocks']))           
			
			% REPOSITORY-SPECIFIC PATH
%			addpath('./ProcessingBlocks')
 			addpath(genpath(['./ProcessingBlocks']));
		    addpath(genpath(['./utilities']));			
 			addpath(genpath(['./simulators']));
			addpath(genpath(['./Plotters']));
 			addpath(['../data']);
            addpath(['../data/unsyncd_data']);
 			addpath(['../data/Lundin_data_for_normal_operations']);
% 			addpath(['../data/MHWirth_data']);
% 			addpath(['../data/IndUrb/']);
		end
	end
	
end

