classdef VARSTFunctionGenerator < STFunctionGenerator
    % This class generates space time functions
    % according to the vector autoregressive model
    
    properties
        t_coefficients % N x N x ORDER, where N is the number of time 
		               % series and ORDER is the model order. 
		m_DCTerm = [] % N x S matrix, where S is the number of seasons in a 
		              % periodic VAR model.  The offset term for the t-th
		              % data sample is the mod(t,S)-th column of m_DCTerm.
		              % Thus, if  S = 1, then the offset term is always
		              % m_DCTerm and the periodic VAR model boils down to a
		              % time-invariant model
        sigma         % standard deviation of the process white Gaussian noise
		        %(rename to be noisePower)
        b_discardSamples=1 % a binary variable that controls whether to discard the 
                           % initial samples
        b_smoothTransition=0 % A binary variable to indicate that the process is Smooth 
                             % Transition nonstationary VAR process
        transitionParameter  % gamma > 0 which determines the speed of the transition
        v_cCoefficients      % cCoefficients in the Exponential Smooth-Transition model
        TransitionPoint      % The point from which the transition starts
        t_coefficientsPlus   % A+ (regime)coefficients in the model
    end
    
    methods % required by superclass ProcessingBlock
                
		function ch_out = printProperty(obj,ch_propertyName,b_justName)
			% ch_out is a string displaying the value of the property whose
			% name is ch_propertyName. 
			%
			% b_justName (optional; default = 0)
			% if b_justName = 1, just the name of the property is printed.
			% Otherwise, both the name and value must be given.
			%
			% Example: obj has a property with name numberOfNodes with
			% value 4. 
			%
			%    >> obj.printProperty('numberOfNodes')
			%
			%        'N = 4'
			%
			%    >> obj.printProperty('numberOfNodes',1)
			%
			%        'N'
			%
			%
			if nargin<3
				b_justName = 0;
			end
			
			switch(ch_propertyName)
				case 'nTimeSamples'
					if b_justName
						ch_out = 'N';
					else
						ch_out = sprintf('N = %d',obj.nTimeSamples);
					end
				case 'sigma'
					if b_justName
						ch_out = '\\sigma';
					else
						ch_out = sprintf('\\sigma = %g',obj.sigma);
					end
			end
			
		end
		
		function c_out = propertiesToPrint(obj)		
			
			c_out = {'nTimeSamples','sigma'};
			
		end
		
	end
	
	methods
		
		function obj = VARSTFunctionGenerator(varargin)
			obj@STFunctionGenerator(varargin{:});
		end 
		
        function [m_signal, m_signalTrain, m_signalTest, SNR] = realization(obj,nAdditionalSamples)
			%
			% Generates a realization of a space-time function adhering to
			% the VAR model
			%
			% Input:
			%   nAdditionalSamples : (optional) is an integer (default
			%              nAdditionalSamples = 0). This parameter is
			%              useful for estimating the prediction error: the
			%              first obj.nTimeSamples columns of m_signal can
			%              be used as training samples, and the last
			%              nAdditionalSamples as test samples. 
			%
			%
			% Output: 
			%   m_signal : size(t_coefficients,1) x (obj.nTimeSamples +
			%              nAdditionalSamples)			
			%   m_signalTrain : size(t_coefficients,1) x (obj.nTimeSamples)
			%              containing only the training samples; 
			%   m_signalTest : size(t_coefficients,1) x nAdditionalSamples
			%              containing only the test samples. If
			%              nAdditionalSamples is 0, m_signalTest will be an
			%              empty matrix. 
			%
			assert(not(isempty(obj.t_coefficients)));
			assert(not(isempty(obj.nTimeSamples)),'obj.nTimeSamples needs to be set');
            if nargin<2
                nAdditionalSamples = 0;
                m_signalTest = [];
            end
            
            nTimeSeries = size(obj.t_coefficients, 1);
            if isempty(obj.m_DCTerm)
                obj.m_DCTerm = zeros(nTimeSeries);
            end
            nSeasons = size(obj.m_DCTerm,2);
            order = size(obj.t_coefficients, 3);
            nTotalTimeSamples = obj.nTimeSamples+nAdditionalSamples;
            my_sigma = obj.sigma;
            % STVAR model
            if obj.b_smoothTransition==1
                if obj.b_discardSamples == 0
                    nSamplesToDiscard=0;
                else
                    nSamplesToDiscard = 10*nTimeSeries*order; % this number should be less than transition point
                end
                m_signal1=zeros(nTimeSeries,nTotalTimeSamples+nSamplesToDiscard);
                m_signal_beforeNoise = m_signal1;
                m_signal1(:,1:order)=my_sigma*randn(nTimeSeries,order);
                m_x=zeros(nTimeSeries,(nTotalTimeSamples+nSamplesToDiscard)-order);
                for t=order+1:(nTotalTimeSamples+nSamplesToDiscard)
                    if t< obj.TransitionPoint
                        m_x(:,t)=obj.v_cCoefficients;
                    else
                        m_x(:,t)=t.*ones(nTimeSeries,1);
                    end
                    G=diag(1-exp(-obj.transitionParameter*(obj.m_x(:,t)-obj.v_cCoefficients).^2));
                    buffer=zeros(nTimeSeries,1);
                    for p=1:order
                        buffer=buffer+obj.t_coefficients(:,:,p)*m_signal1(:,t-p)+ ...
                            G*obj.t_coefficientsPlus(:,:,p)*m_signal1(:,t-p);
                    end
                    m_signal_beforeNoise (:,t) = buffer;
                    m_signal1(:,t)=buffer+my_sigma*randn(nTimeSeries,1);
                end
                m_signal = m_signal1(:,(nSamplesToDiscard+1):end);
            % VAR model   
            else
                if obj.b_discardSamples == 0
                    nSamplesToDiscard=0;
                else 
                    nSamplesToDiscard = 10*nTimeSeries*order*size(obj.m_DCTerm,2); % it must be a multiple of the number of seasons
                end
                
                m_signal1=zeros(nTimeSeries,nTotalTimeSamples+nSamplesToDiscard);
                m_signal_beforeNoise = m_signal1;
                m_signal1(:,1:order)=my_sigma*randn(nTimeSeries,order);
                
                for t=order+1:(nTotalTimeSamples+nSamplesToDiscard)
                    buffer=zeros(nTimeSeries,1);
                    for p=1:order
                        buffer=buffer+obj.t_coefficients(:,:,p)*m_signal1(:,t-p);
                    end
                    if ~isempty(obj.m_DCTerm)
                        buffer = buffer + obj.m_DCTerm(:,mod(t,nSeasons)+1);
                    end
                    m_signal_beforeNoise (:,t) = buffer;
                    m_signal1(:,t)=buffer+my_sigma*randn(nTimeSeries,1);
                end
                m_signal = m_signal1(:,(nSamplesToDiscard+1):end);
                % We keep only the last T generated data
                % (to avoid estimating the matrices from transient data)
                m_signalTrain = m_signal(:,1:obj.nTimeSamples);
                if nAdditionalSamples>0
                    m_signalTest = m_signal(:,obj.nTimeSamples+1:end);
                end
                if nargout>3
                    SNR = 1 + mean(norms(m_signal_beforeNoise(:,(nSamplesToDiscard+1):end)).^2)...
                        /(nTimeSeries*my_sigma.^2);
                end
            end
		end
		
		function [out, score, badness] = isStable(obj, ngp, tol)
			switch nargin
				case 1, ngp = 15; tol = 0.2;
				case 2,	tol = 0.2;
			end
			[out, score, badness] = obj.coefficientsAreStable(...
				obj.t_coefficients, ngp, tol);
		end
		
		function N = numberOfNodes(obj)
			N = size(obj.t_coefficients, 1);
		end
	end
	
	methods (Static)
		function [out, score, badness] = coefficientsAreStable(A, tol)
			%
			% Input:
			%   A  : N x N x order tensor with the VAR coefficients, where
			%        N is  the number of time series
			%  
			if nargin == 1
				tol = 0.02;
			end
			assert(ndims(A)==3);
			[N, ~, P] = size(A);
			m_bigA = zeros(P*N);
			upperRow = reshape(A, N, N*P);
			m_bigA(1:N, :) = upperRow;
			m_bigA(N+1:end, 1:N*(P-1)) = eye(N*(P-1));
			my_eigs = eigs(m_bigA);
			score = max(abs(my_eigs));
			if score<1-tol
				out = true;
			else, out = false;
			end
			badness = score - 1;
		end
		
		function t_normalizedCoefficients = randomCoefficients(filterOrder,m_mask)
			
			%
			% Input:
			%    m_mask : N x N matrix. If m_mask(m,n)==0, then
			%             t_normalizedCoefficients(m,n,p) is 0 for all p.
			%             
			% Output:
			%    t_normalizedCoefficients : N x N x filterOrder matrix with
			%             the VAR coefficients
			%
			t_coefficients = zeros([size(m_mask) filterOrder]);
			
			myEigs = zeros(1, filterOrder);
			for p = 1:filterOrder
				t_coefficients(:,:,p)= randn(size(m_mask))...
					.*logical(m_mask);
				myEigs(p)=abs(eigs(t_coefficients(:,:,p),1));
			end
			t_normalizedCoefficients = t_coefficients;
 			[~, score] = VARSTFunctionGenerator.coefficientsAreStable(t_normalizedCoefficients);
			
			while score > 0.95
				t_normalizedCoefficients = t_normalizedCoefficients*0.9/score;
				[~, score] = VARSTFunctionGenerator.coefficientsAreStable(t_normalizedCoefficients);
			end
		end
		
		function t_normalizedCoefficients = randomCoefficientsFromGraph(...
				graph, filterOrder)
			% this function uses the adjacency matrix of GRAPH as a mask
			% for generating random VAR coefficients. See also
			% randomCoefficients(). 
			%
			
			m_adjacency = graph.m_adjacency; 			
			t_normalizedCoefficients = VARSTFunctionGenerator.randomCoefficients(filterOrder,m_adjacency);
			
        end
        function adjacencyMatrix = VARCoefficients2AdjacencyMatrix(tensor_A)
            % This method computes the adjacency matrix from VAR
            % coefficients (not completely debugged yet!) bakht
            nd=ndims(tensor_A);
            if nd<3
                error('Tensor A must be 3-way')
            elseif nd==3
                [N,N2,P]=size(tensor_A);
                assert(N==N2);
                adjacencyMatrix= reshape(norms(reshape(tensor_A,[N*N, P]), 2, 2), [N,N]);
            else
                error('Use VARCoefficients2AdjacencyMatrices method if ndims>3')

            end
        end
        function adjacencyMatrix = VARCoefficients2AdjacencyMatrices(tensor_A)
            % This method computes the adjacency matrix from VAR
            % coefficients (not completely debugged yet!) bakht
            nd=ndims(tensor_A);
            if nd<3
                error('Tensor A should be atleast 3-way')
            elseif nd==3
                warning('Use VARCoefficients2AdjacencyMatrix instead  ')
                [N,N2,P]=size(tensor_A);
                assert(N==N2);
                adjacencyMatrix= reshape(norms(reshape(tensor_A,[N*N, P]), 2, 2), [N,N]);
            else
%                 tensorB=permute(tensor_A,[1 2 4:nd 3]);
%                 vs= size(tensorB);
%                 assert(vs(1)==vs(2))
%                 permutedAdjacencyMatrix= reshape(norms(reshape(tensorB,[vs(1)*vs(2)*prod(vs(4:end)) vs(3)]), 2, 2), [vs(1),vs(2), vs(4:nd)]);
%                 adjacencyMatrix=ipermute(permutedAdjacencyMatrix, [1 2 4:nd]);
                [N,N2,P]=size(tensor_A);
                assert(N==N2);
                adjacencyMatrix= norms(tensor_A,2,3);
            end
        end
        
		
% 		function [out, score, badness] = coefficientsAreStable_old (...
% 				A, ngp, tol)
% 			switch nargin
% 				case 1, ngp = 15; tol = 0.2;
% 				case 2,	tol = 0.2;
% 			end
% 			z = sym ('z');
% 			polynomial = eye(size(A, 1));
% 			P = size(A,3);
% 			for p = 1:P
% 				polynomial = polynomial - A(:,:,p)*z^p;
% 			end
% 			[RE, IM] = meshgrid(linspace(-(1+tol), 1, ngp), linspace(-1, 1, ngp));
% 			fcomplex = @(re, im) subs(det(polynomial), (re+im*1i));
% 			grid = double(fcomplex(RE, IM));
% 			mask = logical(abs(RE+IM*1i)<=1+tol);
% 			score = min(abs(grid(mask)));
% 			
% 			draw = 0;
% 			if draw
% 				figure
% 				surf(double(RE), double(IM), abs(grid)./mask)
% 				%subplot(122)
% 				%surf(double(RE), double(IM), imag(grid)./mask)
% 			end
% 			
% 			badness = -log(score/tol);
% 			%out =  badness > 0;
% 			out = score > tol;
% 		end
	end
end