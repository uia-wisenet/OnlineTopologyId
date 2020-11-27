classdef ExperimentFunctionSet
	%
	% EXPERIMENT FUNCTION SET Superclass for collections of experiment 
	% functions to be called by GSIM.
	%
	% A class inheriting from ExperimentFunctionSet contains functions
	% whose name starts by the string in funbasename.
	%
	% These functions return an object of the class GFigure. 
	%
	
	properties(Constant)
		funbasename = 'experiment_';	
	end
		
	properties
		classname
		figbasename
		runningOnServer_property = 0;
	end
	
	methods(Access=public)
		
		% constructor
		function obj = ExperimentFunctionSet
			obj.classname=class(obj);
			obj.figbasename = [obj.outputFolder obj.classname ];
		end
		
		% main method
		function gsimExecute(obj,experimentIndices,niter,onlyplot, runningOnServer)
			% This function executes the experiments in file obj with index
			% in the vector experimentIndices.
			
			% Initial settings
			global displaySettings
			displaySettings.figbasename= obj.figbasename;
				
			% Create output folder if it does not exist
			if (~exist(obj.outputFolder,'file'))
				mkdir(obj.outputFolder);
			end
			
			% -------------
			tstart = tic;
			
			for k=1:length(experimentIndices)
				fprintf('Starting experiment %d in %s.\n',experimentIndices(k),obj.classname);
				obj.runExperiment(experimentIndices(k),niter,onlyplot, runningOnServer);
			end
			
			tm = toc(tstart);
			% -------------
			fprintf('Elapsed time: %s\n',print_time(tm));
			
		end

		% auxiliary methods
		function str = outputFileName(obj,experimentIndex)
			% name without extension and location of the mat/fig/pdf/...
			% files corresponding to experiment <experimentIndex>
			filebasename=[obj.outputFolder obj.classname '_'];  
			str=sprintf('%s%d',filebasename,experimentIndex);		
		end
		
		function str = outputFolder(obj)
			str = ['Experiments/' obj.classname '_data/'];
		end
		
		function F = loadGFigure(obj,experimentIndex)
			if (~exist([obj.outputFileName(experimentIndex) '.mat'],'file'))
				error('You must run the simulation first. Use gsim(0,%d)',experimentIndex);
			end
			assert(numel(experimentIndex)==1)
			stru = load(obj.outputFileName(experimentIndex));
			F = stru.F;
		end
				
		
	end
	
	methods(Access=private)
		
		function runExperiment(obj,experimentIndex,niter,onlyplot, runningOnServer)
			% runs an experiment if onlyplot = 0; otherwise, it only plots
			% the results.
				
			funname =sprintf('%s%d',obj.funbasename,experimentIndex);
			
			if ~ismethod(obj,funname)
				error('Experiment %d does not exist in %s',experimentIndex,obj.classname);
			end
			
			if onlyplot
				% load results from previous executions
				F = obj.loadGFigure(experimentIndex);
			else
				% run the experiment and save the results
				fprintf('niter = %d\n',niter)					
				F=feval(funname,obj,niter);
				save(obj.outputFileName(experimentIndex),'F');				
			end
			
			if runningOnServer % No plotting if running on a server.
				if onlyplot
					warning ('Your configuration: onlyplot=1, runningOnServer=1');
				end
				return
			end
			
			if isa(F,'GFigure')
				F.plot(experimentIndex);
			else
				fprintf('Experiment %d does not create any figure.\n',experimentIndex)
			end
			
		end
		
	end
	
	% Methods for backwards compatibility	
	methods
		
		function gsimCompute(varargin) 
			disp('Your gsim.m file is out of date. Run gsimStartup to update.')
		end
		
	end
	
end

