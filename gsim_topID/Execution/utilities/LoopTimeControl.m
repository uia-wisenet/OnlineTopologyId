classdef LoopTimeControl < handle;
	%Loop Time Control monitor. This class allows to monitor the 
	%execution time for very time-consuming loops.
	% HOW TO USE:
	% ltm = LoopTimeControl(nIter) constructs a new object to monitor the
	% execution time of a loop whose number of iterations is nIter.
	% The object should be constructed right before starting the loop.
	%
	% ltm.go(itIndex) updates the information about the runtime and shows a
	% message in the command line informing the user about the estimated
	% time to complete the task.
	% This method should be called right before the loop's end clause.
	%
	% By default and to reduce overhead, the monitor will update its 
	% estimate in time intervals between 0.1 s and 1 s. If you want to 
	% unable this feature, create the monitor using:
	% ltm = LoopTimeControl(nIter, 'all') 
	% This can be useful, for example, when warnings appear in certain 
	% iterations. Each warning message will appear next to a line 
	% informing of the iteration count.
	% 
	% EXAMPLE:
	% % I want to compute the cholesky decomposition of a random
	% % matrix 500 times.
	% ltm = LoopTimeControl(500);
	% for i = 1:500
	%     A = rand(1000);
	%     B = inv(A);
	%     ltm.go(i);
	% end
	%
	% In future developments, this will be renamed as LoopTimeMonitor.
	
	% Author: Luis Miguel Lopez, 2018
	
	properties
		tBegin
		tLast
		elapsedTime
		eta
		myString = [];
		totalIter;
		b_displayAll = 0;
	end
	
	methods
		function obj = LoopTimeControl(niter, str_in)
			obj.tBegin = tic;
			obj.tLast = tic;
			obj.totalIter = niter;
			if exist('str_in', 'var') && isequal(str_in,'all')
				obj.b_displayAll = 1;
			end
		end
		
		function go(obj, i)
			estTotal = obj.elapsedTime+obj.eta;
			if obj.b_displayAll
				timeQuantum = 0;
			else
				try
					timeQuantum = max(0.1, min(1, 0.01*estTotal* ...
						1-cos(2*pi*obj.elapsedTime/estTotal)));
					assert(not(isempty(timeQuantum)));
				catch
					timeQuantum = 0.1;
				end
			end
			if toc(obj.tLast) > timeQuantum || i==obj.totalIter
				obj.tLast = tic;
				obj.elapsedTime = toc(obj.tBegin);
				obj.eta = (obj.totalIter - i) * obj.elapsedTime ./ i;
				fprintf(repmat('\b', 1, length(obj.myString)));
				obj.myString = sprintf(...
					'Iteration %0.4g out of %d. ETA = %s', ...
					i, obj.totalIter, print_time(obj.eta));
				% TODO: the fprintf should be printed only if the change in
				% the elapsed time is significant.
				fprintf(obj.myString);
				if i==obj.totalIter
					fprintf(';; Done!\n');
				end
			elseif i>obj.totalIter
				warning('call to GO method beyond expected nIter.')
			end
		end
	end
	
end

