classdef OnlineGroupLassoEstimator < handle
		
	properties
		
		regPar % regularization parameter
        noOfGroups
        groupSize
        t         %iteration number
		n_no_regularize = 0 % TODO: not all subclasses implement this yet!
		
		v_coefficients % estimated coefficients
	end
	
	
	methods (Abstract)
		
		% UPDATE receive training sample and update coefficients.
		% myEstimator.update(v_x, y) uses the regressor-output pair to
		% update the estimated coefficients.
		% v_newCoefs = myEstimator.update(v_x, y) outputs the new
		% coefficients after computation. Also, the coefficients can
		% always be directly accessed as myEstimator.v_coefficients.
		[v_newCoefficients,objective] = update(obj, v_x, y)
		
		% INITIALIZE create necessary structures in the object to 
		% get it ready to execute the UPDATE method.
		initialize (obj)
			
	end
	
	methods
		
		% PREDICT return predicted value for an input, given the
		% current state of the estimator object.
		% y = myEstimator.predict(v_x) computes the predicted output y
		% for the input v_x as y = obj.v_coefs^T * v_x.
		function y = predict(obj, v_x)
			y = obj.v_coefficients'*v_x;
			
		end		
	end
	
end

