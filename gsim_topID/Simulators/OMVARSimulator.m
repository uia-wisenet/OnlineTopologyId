function [sq_error,pred_error,samplePower,t_A, ...
	ratioEdgeIdError, ratioDetectedEdge, ratioFalseAlarm] =  ...
	OMVARSimulator (generator, estimatorTemplate,...
	filtOrder, noOfObservations, trueCoefficients, threshold, nonStationary)

if not(exist('nonStationary', 'var'))
    trueAdjacencyMatrix = VARSTFunctionGenerator.VARCoefficients2AdjacencyMatrix(...
        trueCoefficients);
else
    % for nonstationary experiments
    trueAdjacencyMatrix = VARSTFunctionGenerator.VARCoefficients2AdjacencyMatrices(... %!(for recomidexp1007)
        trueCoefficients);
end
if not(exist('threshold', 'var'))
	threshold = 0.15;
end
m_X = generator.realization();
noOfNodes = size(m_X,1); %number of time series (nodes)
noOfPotentialEdges = noOfNodes.^2; 
bm_trueAdj = logical(trueAdjacencyMatrix);
noOfTrueEdges = nnz(bm_trueAdj);

if strcmp(estimatorTemplate, 'genie')
	b_genie = 1; %simulate genie-aided estimator
else
	b_genie = 0;
	ove = OnlineModularVAREstimator;
	ove.filterOrder = filtOrder;
	ove.nTimeSeries = noOfNodes;
	createOgles(ove, estimatorTemplate);       
	% create v_scalarEstimator (a vector of estimators) inside the OVE
	inheritProperties(ove, estimatorTemplate); 
	% copy the rest of the properties from template
	initializeOgles(ove); % call initialize method for all estimators
end

v_coeff=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
sq_error          = zeros(1, noOfObservations);
pred_error        = zeros(1, noOfObservations);
samplePower       = zeros(1, noOfObservations);
ratioDetectedEdge = zeros(1, noOfObservations);
ratioEdgeIdError  = zeros(1, noOfObservations);
ratioFalseAlarm   = zeros(1, noOfObservations);
t_A= zeros(noOfNodes,noOfNodes,noOfObservations);
for i = 1:noOfObservations
	if numel(size(trueCoefficients))==3
		trueCoefficients_now = trueCoefficients;
	else
		trueCoefficients_now = trueCoefficients(:,:,:,i);
	end
	
	if b_genie
		sq_error(i) = nan;
		try
			m_X_predicted = OnlineModularVAREstimator.predict_static(...
				trueCoefficients_now, m_X(:, i-filtOrder:i-1),1);
			pred_error(i) = norm(m_X_predicted - m_X(:,i))^2;
			samplePower(i) = norm(m_X(:,i))^2;
		catch ME
			if strcmp(ME.identifier, 'MATLAB:badsubscript') && ...
					i-filtOrder<1
				% OK
			else
				rethrow(ME)
			end
			
		end
		
	else
		try
			setTimeIndex(ove, i);
			%     for i = 1:ove.nTimeSeries
			%         ove.v_scalarEstimator(i).t = iter;
			%     end
			
			v_coeff(:,:,:,i)=ove.t_estimatedCoefficients;
			t_A(:,:,i)= VARSTFunctionGenerator.VARCoefficients2AdjacencyMatrix (v_coeff(:,:,:,i));
			sq_error(i) = norm(vec(v_coeff(:,:,:,i)-trueCoefficients_now))^2;
			% 		if numel(size(trueCoefficients))==3 % for stationary setting (constant coefficients)
			% 			sq_error(iter)=norm(vec(v_coeff(:,:,:,iter)-trueCoefficients))^2;
			% 		else
			% 			sq_error(iter)=norm(vec(v_coeff(:,:,:,iter)-trueCoefficients(:,:,:,iter)))^2;
			% 			% Bakht: (this is not compatible with stationary VAR
			% 			%         since the true estimates are changing here)
			% 		end
			m_X_predicted = ove.predict_static(v_coeff(:,:,:,i), ...
				m_X(:,i-ove.filterOrder:i-1),1);
			pred_error(i)= norm( m_X_predicted - m_X(:,i))^2;
			samplePower(i) = norm(m_X(:,i))^2;
			ratioEdgeIdError(i) = nnz(...
				bm_trueAdj-logical(t_A(:,:,i)>threshold))/noOfPotentialEdges;
			ratioDetectedEdge(i) = nnz(...
				bm_trueAdj.*logical(t_A(:,:,i)>threshold))/noOfTrueEdges;
			ratioFalseAlarm(i)  = nnz(...
				not(bm_trueAdj).*logical(t_A(:,:,i)>threshold)...
				)/(noOfPotentialEdges-noOfTrueEdges);
					
		catch ME
			if strcmp(ME.identifier, 'MATLAB:subsassigndimmismatch') && ...
					isempty(ove.t_estimatedCoefficients)
				% OK
			else
				rethrow(ME)
			end
		end
		
		receiveSingleSample(ove, m_X(:,i));

	end
	
	
	
end
end