classdef OnlineModularVAREstimator < handle % & ProcessingBlock
	%ONLINEVARESTIMATOR Estimator of the parameters of a VAR 
	%(vector autoregressive) process using multiple OGLEs 
	% (objects inheriting from OnlineGroupLassoEstimator)
	
	properties
		%One OnlineEstimator object for each time series
		v_scalarEstimator 
		filterOrder %P
		nTimeSeries %N
		m_processBuffer %NxP matrix
		bufferLevel = 0
		
		b_leastSquareSelfLoops = 0;
		
		% TODO: this property should be removed
		% (it will be replaced with getEstimatedCoefficients method
		% that looks inside all scalar estimators one by one)
		t_estimatedCoefficients % NxNxP tensor
		%(:,:,p) = A^{(p)}
	end
	
	methods
		% function b_OK = isConsistent(obj)
		% end
		
		% create Online Group Lasso Estimators that will allow
		% to estimate the VAR process with edge sparsity.
		% Note that the OnlineVAREstimator should also work with any 
		% other kind of Online Estimator.
		function [] = createOgles(obj, ogleTemplate_in)
			assert (isa(ogleTemplate_in, 'OnlineGroupLassoEstimator'));
			v_ogle = feval(class(ogleTemplate_in)); %vector of estimators 
			% of the same class as ogle_in
			for n = obj.nTimeSeries:-1:1
				%new object in the nth entry
				ogle = feval(class(ogleTemplate_in)); 
				% inherit properties in the (super)class
				ogle.regPar = ogleTemplate_in.regPar;
				ogle.noOfGroups = obj.nTimeSeries;
				ogle.groupSize = obj.filterOrder;
				if obj.b_leastSquareSelfLoops
					ogle.n_no_regularize = n;
				end
				ogle.t = 1;
				v_ogle(n) = ogle;
			end
			obj.v_scalarEstimator = v_ogle;
		end
		
		function [] = inheritProperties(obj, ogle_template)
			% for the scalar estimators to inherit all properties
			% that are not defined in the superclass.
			% (for example, forgetting factor)
			propertyList = setdiff(...
				properties(ogle_template), ...
				properties('OnlineGroupLassoEstimator'));
			for i = 1:length(propertyList)
				propertyName = propertyList{i};
				if not(isempty(ogle_template.(propertyName)))
					for n = obj.nTimeSeries:-1:1
						obj.v_scalarEstimator(n).(propertyName) ...
							= ogle_template.(propertyName);
					end
				end
			end
		end
		
		function [] = initializeOgles(obj)
			for n = obj.nTimeSeries:-1:1
				ogle = obj.v_scalarEstimator (n);
				initialize(ogle);				
			end
		end
		
		function setTimeIndex(obj,t_in)
			for i = 1:obj.nTimeSeries
				obj.v_scalarEstimator(i).t = t_in;
			end
		end
		
		% receive a sample corresponding to a single time instant,
		% place it in the buffer, and call the function that updates
		% all of the OnlineEstimators.
		% NOTE: the time index is not updated in this method. 
		% If this is needed, call the setTimeIndex method.
		function [] = receiveSingleSample(obj, v_sample_in)
			assert(isvector(v_sample_in) && ...
				length(v_sample_in) == obj.nTimeSeries, ...
				'v_sample_in must be a vector of length nTimeSeries');
			
            if isempty(obj.m_processBuffer) || obj.bufferLevel < 0
                % create buffer matrix and set level to 0
                obj.m_processBuffer = zeros(obj.nTimeSeries, ...
					obj.filterOrder);
                obj.bufferLevel = 0;
                
            elseif obj.bufferLevel < obj.filterOrder
                % add received sample at the end; increase buf level
                obj.m_processBuffer(:, ...
					obj.filterOrder-obj.bufferLevel) = v_sample_in;
                obj.bufferLevel = obj.bufferLevel + 1 ;
                
            elseif obj.bufferLevel == obj.filterOrder
                assert(length(obj.v_scalarEstimator) == ...
					obj.nTimeSeries, ['The number of '...
					'scalarEstimators must equal nTimeSeries']);
                
                updateEstimators(obj, v_sample_in); 
                % shift buffer and add current sample at the end
                obj.m_processBuffer(:, 2:end) = ...
                    obj.m_processBuffer(:, 1:end-1);
                obj.m_processBuffer(:, 1) = v_sample_in;
                
            else
                error('Fatal: buffer level exceeds filter order');
            end

		end
		
		function m_A = estimatedAdjacencyMatrix(obj)
			if isempty(obj.t_estimatedCoefficients)
				m_A = zeros(obj.nTimeSeries);
			else
				m_A = VARSTFunctionGenerator.VARCoefficients2AdjacencyMatrix(...
				obj.t_estimatedCoefficients);
			end
		end
		
		function t_c = getEstimatedCoefficients(obj)
			N = obj.nTimeSeries;
			t_c = zeros(N,N, obj.filterOrder);
			for n = 1:obj.nTimeSeries
				t_c(n,:,:) = reshape(obj.v_scalarEstimator(n ...
					).v_coefficients, [N, obj.filterOrder]);
			end
		end
	end
	
	methods(Access=private)
		% call update function for each scalar estimator
		function [] = updateEstimators(obj, v_sample_in)

			% We assume consistency has been checked before calling
			% this function (see receiveSingleSample)
			for n = 1:obj.nTimeSeries
				y = v_sample_in(n);
				%v_x = vec(obj.m_processBuffer');
                aux = obj.m_processBuffer';
                v_x= aux(:);
				
				update(obj.v_scalarEstimator(n), v_x, y); 
				% v_x is the same for all y's
				
				% This function had a bug in the previous version:
% 				v_x = vec(obj.m_processBuffer);
% 				update(obj.v_scalarEstimator(n), v_x, y); 
% 				obj.t_estimatedCoefficients(n, :,:) = reshape(...
% 					obj.v_scalarEstimator(n).v_coefficients, ...
% 					[obj.nTimeSeries, obj.filterOrder]);
				
				obj.t_estimatedCoefficients(n, :,:) = reshape(...
					obj.v_scalarEstimator(n).v_coefficients, ...
					[obj.filterOrder, obj.nTimeSeries])';
				
				% TODO: substitute these lines with a call to the
				% getEstimatedCoefficients method.
			end
		end
	end
	
	
	methods (Static)
		function m_predictedSignalOut = predict_static( ...
			t_Coefficients, m_pastValues, numInstantsToPredict)

		[N, N_rep, P] = size(t_Coefficients);
		T_past = size(m_pastValues,2); % number of past values
		assert(T_past >= P, ['Number of past values should be ', ...
			'at least P(order of the VAR process)']);
		assert(N==N_rep, 'slices of t_estimatedCoefs are not square');
		m_predictedSignalAux = zeros(N, numInstantsToPredict);
		m_X = [m_pastValues m_predictedSignalAux];
		m_A = reshape(t_Coefficients,N,N*P);
		for colNow = T_past+(1:numInstantsToPredict)
			m_X(:,colNow) = m_A*vec(m_X(:,colNow-1:-1:colNow-P));				
		end
		m_predictedSignalOut = m_X(:,T_past+1:end);
	end	
	end
	
end