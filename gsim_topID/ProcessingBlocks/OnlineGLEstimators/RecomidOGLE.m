classdef RecomidOGLE < OnlineGroupLassoEstimator
	% RecomidOGLE computes the solution to online group lasso problem
	% in an online fashion, by providing the COMID algorithm with a loss
	% function similar to that of RLS, exponentially windowed sum of
	% the past data (in fact, each iteration considers all past data).
	% If we do not want it recursive then we simply make
	% obj.b_selectComid = 1 
	% and
	% obj.forgettingFactor = 0.
	properties
		% state properties
		m_Phi                 % autocorrelation matrix
		v_r                   % cross-correlation vector
		forgettingFactor
		sigma                 % to initialize Phi
		
		stepSize=1;           % scale factor (if adaptive stepsize,
		% will be used only in the first iterations; if diminishing
		% stepsize, will be divided by sqrt(t)
		b_adaptiveStepsize=0; % to activate adaptive stepsize
		b_diminishingStepsize=1; % to divide the stepsize by sqrt(t)
		t_startAdaptingStepsize = 5; % before this instant,
		% the stepsize stays constant (because we do not want to start very
		% early with the adaptation)
		ss_factor = 2;        % scale factor for the (adaptive) stepsize.
		                      % Should be greater than or equal to 2.
		
		mulFactor=1           % to cancel the effect of recursive
		% objective in order to compare with Batch alg.
		% (\mu in the text)
		b_selectComid=0;      % 1 when Recomid is used as Comid
		rho=0;                % used when selecting adaptive stepsize
		% in Comid (for finding max we need the previous rho)
		
		b_maxInComid=0        % use the max operator when
		% determining the COMID stepsize
	end
	
	methods
		function initialize(obj)
			obj.v_coefficients=zeros(obj.noOfGroups*obj.groupSize,1);
			obj.m_Phi=obj.sigma.^2*eye(obj.noOfGroups*obj.groupSize);
			obj.v_r=zeros(obj.noOfGroups*obj.groupSize,1);
			if not(obj.b_selectComid)
				if isempty(obj.mulFactor)
					obj.mulFactor = 1 - obj.forgettingFactor;
					warning 'setting mulFactor = forgettinFactor'
				else
					if obj.mulFactor+obj.forgettingFactor ~= 1
						warning 'mulFactor not equal to 1-forgettingFactor'
					end
				end
			else
				assert(obj.forgettingFactor==0);
			end
				
		end
		
		function [v_newCoefficients,objectiveValue] = update(obj, v_x, y)
			updatePhiandr(obj,v_x,y)
			v_newCoefficients=zeros(size(obj.v_coefficients));
			%% determine eta (stepsize-related ReCOMID parameter)
			if obj.b_adaptiveStepsize				
				if obj.b_selectComid==1
					% first version:
					% obj.rho=max(obj.rho,0.5*(norm(v_x)^2));
					if obj.b_maxInComid
						obj.rho = max(obj.rho, obj.ss_factor*norm(v_x)^2); %!
					else
						obj.rho = obj.ss_factor*norm(v_x)^2; %!
					end
				else
					% first version:
					% obj.rho=max(eig(obj.m_Phi));
					obj.rho=obj.ss_factor*max(eig(obj.m_Phi));
					% TODO: theoretical explanation for this
				end
				
				if obj.t < obj.t_startAdaptingStepsize
					% do not change the stepsize
				else
					obj.stepSize=2/obj.rho;
				end
			end
			if obj.b_diminishingStepsize
				% stepsize = constant/sqrt(t)
				alpha = obj.stepSize/sqrt(obj.t);
			else
				alpha = obj.stepSize;
			end
			eta = 1/alpha; % For simplicity, in our calculation eta=1/alpha;
			
			%% compute subgradient and update estimate
			v_subgrad = obj.computeSubgradient(); % gradient in this case
			for g=1:obj.noOfGroups                % solve the problem for each group
				indices=((g-1)*obj.groupSize+1):(g*obj.groupSize);
				b=eta*obj.v_coefficients(indices)-v_subgrad(indices);
				a = eta;
				v_newCoefficients(indices)= obj.ECOMIDRec_group(...
					b,a, obj.regPar*(g~=obj.n_no_regularize), numel(indices));
				% implemented by the given closed form
			end
			obj.v_coefficients=v_newCoefficients;
			objectiveValue=computeObjectiveValue(obj);
		end
		
		
		function updatePhiandr(obj,x,y)
			%update Phi and r after getting data sample
			obj.m_Phi = obj.forgettingFactor*obj.m_Phi + obj.mulFactor*(x*x');
			obj.v_r   = obj.forgettingFactor*obj.v_r   + obj.mulFactor*y*x;
		end
		
		
		function sg_out = computeSubgradient(obj) %gradient in this case
			%computing the gradient of {1/2 w^TPhiw-r^Tw}
			sg_out=obj.m_Phi*obj.v_coefficients - obj.v_r;
		end
		
		function obj_out=computeObjectiveValue(obj)
			glPenality=0;
			for i=1:obj.groupSize:numel(obj.v_coefficients)
				glPenality=glPenality+norm(obj.v_coefficients(i:i+obj.groupSize-1));
			end
			obj_out=0.5*(obj.v_coefficients'*obj.m_Phi*obj.v_coefficients)...
				-obj.v_r'*obj.v_coefficients + obj.regPar*glPenality;
		end
		
		
	end
	
	methods(Static)
		function v_groupCoefficients = ECOMIDRec_group(b,a, lambda, dg)
			% closed-form solution for each group
			% This method solves the general problem of the form
			% argmin_x {1/2 ax^Tx-b^Tx+lambda||x||_2}
			% in closed form
			v_groupCoefficients = (b./a).*max(1-(lambda/norm(b)),0);
			% (note for future versions) variable group size not implemented yet
			
		end
	end
end

