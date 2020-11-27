classdef TimeVaryingVARSTFunctionGenerator
		
properties
	
	N % number of time series 
	Q % model order
	T % total time instants
	
	t_coefficients % N x N x Q x T
	
	sigma % standard deviation of innovations
	
	noiseDistribution = 'Normal'
	noiseParams = {0; 1};
	
	b_backwardsCompatible = 1; % to keep the random edge generation
	% as it used to be, so old simulatios still work.
	% in NEW simulations, SET this to 0!
		
end
	
methods
	function [m_signal] = realization(obj)
		% Generates a realization of a space-time function adhering to
		% the time-varying VAR model
		assert(isequal(size(obj.t_coefficients), [obj.N obj.N obj.Q obj.T]));
		
		m_signal = zeros(obj.N, obj.T);
		try
			processNoise = random(obj.noiseDistribution, ...
				obj.noiseParams{:}, [obj.N, obj.T]);
		catch ME
			switch obj.noiseDistribution
				case 'Laplacian'
					auxNoise = rand(obj.N, obj.T);
					mu = obj.noiseParams{1}; b = obj.noiseParams{2};
					processNoise = mu - b.*sign(auxNoise).*log(1-2.*abs(auxNoise));
				otherwise
					rethrow (ME)
			end
		end
		for t =  1:obj.T
			if t <= obj.Q
				m_signal(:, t) = obj.sigma * processNoise(:, t); %randn(obj.N, 1);
			else
				buffer = zeros(obj.N, 1);
				for p = 1:obj.Q
					buffer = buffer + obj.t_coefficients(:,:,p, t)*m_signal(:, t-p);
				end
				m_signal(:, t) = buffer + obj.sigma * processNoise(:, t); %randn(obj.N, 1);
			end
		end
	end

end

methods (Static)
	
	function gamma_out = stParameter(level, tReach)
		% design smooth transition parameter
		% Gamma = stParameter(level, tReach) returns the parameter
		% Gamma such that the transition coefficient equals level
		% after tReach instants.
		% EXAMPLE:
		% % I want a Gamma parameter such that the transition is 99%
		% % complete after 1000 time instants:
		% gamma = TimeVaryingVARSTFunctionGenerator.stParameter(0.99, 1000)
		gamma_out = -log(1-level)/tReach^2;
	end
	
	function [t_coefficients, t_A,g] = smoothTransitionModel(t_coefs_before, t_coefs_after, paramGamma, paramTb, totalT)
		% generates a smooth transition model according to Killian and
		% L?tkepohl, 2018, sec. 18.3
		% This function takes as input the Coefficients of the two regimes,
		% transtion parameter, the time at which the transition starts, and
		% the total number of samples in STVAR model
		% Output: The output contains the timevarying coefficients of the
		% model and the timevarying adjacency matrix
        assert(not(isempty(t_coefs_before)));
        assert(not(isempty(t_coefs_after)));
		assert(not(isempty(totalT)),'Total samples needs to be set');
        assert(isequal(size(t_coefs_before),size(t_coefs_after)));
        nTimeSeries = size(t_coefs_before, 1);
        filtOrder=size(t_coefs_before,3);
        t_coefficients = zeros([size(t_coefs_before) totalT]);
        t_A= zeros([nTimeSeries, nTimeSeries, totalT]);
        g=zeros(totalT,1);
        for t=1:totalT
           g(t)=1-exp(-paramGamma*(max(t-paramTb,0).^2)); 
           t_coefficients(:,:,:,t)= t_coefs_before + g(t).*(t_coefs_after-t_coefs_before);
           t_A(:,:,t)= VARSTFunctionGenerator.VARCoefficients2AdjacencyMatrix(t_coefficients(:,:,:,t));
        end
		
		
	end
	function [t_coefficients, t_A] = randomCoefficientsFromGraph(g, q, bv_SBInstants, b_plotProcess, b_backwardsCompatible)
		assert(isa(g, 'Graph'));
		
		m_A_base = g.m_adjacency;
		p        = size(m_A_base,1);
		edges    = find(m_A_base(:));
		m_rowis  = (1:p)'*ones(1, p); m_colis = ones(p, 1)*(1:p);
		T        = length(bv_SBInstants);
		pr       = 0.6;

		t_coefficients  = zeros([size(m_A_base) q ]);
		m_randMask      = random('Binomial', 1, pr, size(m_A_base));
		m_A_initial     = m_A_base.*m_randMask;
		
% 		for l = 1:q
% 			t_coefficients(:,:,l,1) = randn(size(m_A_initial)).*logical(m_A_initial);
% 		end
		t_coefficients(:,:,:,1) = VARSTFunctionGenerator.randomCoefficients(q, m_A_initial);
		
		t_A = zeros([size(m_A_base) T]); % sequence of A matrices
		t_A(:,:,1) = m_A_initial;
		for t = 2:T
			m_newA = t_A(:,:,t-1);
			t_newCoefficients = t_coefficients(:,:,:,t-1);
			if bv_SBInstants(t) ~= 0 % change coefs associated with 
				% an edge selected at random
				if not (exist('b_backwardsCompatible', 'var'))
					b_backwardsCompatible = 1
				end
				if b_backwardsCompatible
					myEdge = randi(length(edges)); % this was a bug
				else
					myEdge = edges(randi(length(edges))); %bug fixed
				end
				% If the edge is inactive, activate it;
				% If the edge is active, deactivate with probability pr
				m_newA(myEdge) = not(m_newA(myEdge) && rand < pr);
								
				if m_newA(myEdge)
					% assign random coefficients to the edge
					t_newCoefficients(m_rowis(myEdge), m_colis(myEdge),:) = randn([q 1]);
				else % put the edge coefficients to zero
					t_newCoefficients(m_rowis(myEdge), m_colis(myEdge),:) = 0;
				end	
				if not(VARSTFunctionGenerator.coefficientsAreStable(t_newCoefficients))
					alph = linspace(-1, 1.5, 40);
					score = zeros(1, 40);
					for i = 1:length(alph)
						[~, score(i)] = VARSTFunctionGenerator.coefficientsAreStable...
							((1-alph(i))*t_newCoefficients + alph(i)*t_coefficients(:,:,:,t-1));
					end
					[min_val, arg_min] = min(score);			
					if exist('b_plotProcess', 'var') && b_plotProcess
						figure(9991); clf; 
						plot(alph, score)
						hold on; plot(alph(arg_min), min_val, 'go')
						[~, arg_one] = min(abs(alph-1));
						plot(1, score(arg_one), 'xr');
						drawnow
						%keyboard
					end
					t_newCoefficients = (1-alph(arg_min))*t_newCoefficients + alph(arg_min)*t_coefficients(:,:,:,t-1);
				end
			end
			t_A(:,:,t) = m_newA;
			t_coefficients(:,:,:,t) = t_newCoefficients;
			
		end
	end
	
	function [t_rgb, m_norms, v_coeff_out, v_mu_out] = drawDiagram(t_c, option, v_coeff_in, v_mu_in, excluded_edges)
		warning 'TODO: specify which edges should be excluded from the representation'
		if not(exist('excluded_edges', 'var'))
			excluded_edges = [];
		end
		v_size = size(t_c);
		if length(v_size)==3
			v_size(4) = 1;
		elseif length(v_size) ~= 4
			error 'nDims of t_c should be either 3 or 4'
		end
		
		t_c2 = permute(t_c, [3 1 2 4]);
		m_c = reshape(t_c2, [v_size(3), numel(t_c2)./v_size(3)]);
        v_norms = TimeVaryingVARSTFunctionGenerator.norms(m_c); 
        if not(exist('option', 'var')), option = 1; end
		switch option
			case 1
				[COEFF, SCORE, ~, ~, ~, v_mu_out] = pca(m_c');
				m_score = reshape(SCORE(:, 1), prod(v_size(1:2)),v_size(4));
				v_coeff_out = COEFF(:, 1);
			case 2
				% the norms and scores are probably correlated. 
				% idea to solve this: do pca over the data + the norms and
				% use the second PCA component to select the hue.
				[COEFF, SCORE, ~, ~, ~, v_mu_out] = pca([m_c;v_norms]');
				m_score = reshape(SCORE(:, 2), prod(v_size(1:2)),v_size(4));
				v_coeff_out = COEFF(:,2);
			case 3
				m_score = reshape((m_c'-v_mu_in)*v_coeff_in, ...
					prod(v_size(1:2)),v_size(4));
				v_coeff_out = v_coeff_in;
				v_mu_out = v_mu_in;
			case 4
				[COEFF, ~, ~, ~, ~, v_mu_out] = pca(m_c');
				v_coeff_out = COEFF(:,1); %direction of principal component
				% normalized cosine distance:
				unitary_vecs = (m_c-v_mu_out')./TimeVaryingVARSTFunctionGenerator.norms(m_c-v_mu_out')';
				distances = cos(unitary_vecs'*v_coeff_out);
				m_score = reshape(distances, prod(v_size(1:2)),v_size(4));
			case 5
				v_coeff_out = v_coeff_in;
				v_mu_out = v_mu_in;
				unitary_vecs = (m_c-v_mu_out')./TimeVaryingVARSTFunctionGenerator.norms(m_c-v_mu_out')';
				distances = cos(unitary_vecs'*v_coeff_out);
				m_score = reshape(distances, prod(v_size(1:2)),v_size(4));
			otherwise
				error 'invalid option'
		end
		
		m_norms = reshape(v_norms, prod(v_size(1:2)),v_size(4));		
		t_hsv = ones(prod(v_size(1:2)),v_size(4), 3);
		t_hsv(:,:,1) = (m_score - min(m_score(:)))./(max(m_score(:)) - min(m_score(:)));
		t_hsv(:,:,2) = m_norms./1.2./max(m_norms(:));
		t_hsv(:,:,3) = 1-m_norms./5./max(m_norms(:));
		t_rgb = hsv2rgb(t_hsv);
				
% 		m_normsP = zeros(size(m_norms)+1);
% 		m_normsP(1:end-1, 1:end-1) = m_norms;
% 		figure(9997); clf
% 		pcolor(m_normsP)
% 		axis ij
	end
	
	function h = drawCoefficientsInLocation(v_coefs, edge, t_begin,...
			pseudo_duration, nEdges, nObservations, a, maxCoef, exageration)
		x0 = a.Position(1);
		y0 = a.Position(2);
		b  = a.Position(3);
		h  = a.Position(4);
		if not(exist('exageration', 'var'))
			exageration = 5/3;
		end
		newPosition = [ x0 + (t_begin+1)/nObservations*b;  ...
			y0 + (nEdges-edge-(exageration-1)/2)/nEdges*h;
			pseudo_duration/nObservations*b; exageration*h/nEdges]';
		ax = axes('Position', newPosition);
		stem(v_coefs)
		ylim ([-maxCoef maxCoef])
		xlim([1/2 length(v_coefs)+1/2])
		set(ax,'Visible','off');
		%set(get(ax(j),'Title'),'Visible','on');
    end
    
    function v_norms = norms(m_x)
        v_norms = zeros(size(m_x,2), 1);

        for i = 1:size(m_x,2)
            v_norms(i) = norm(m_x(:,i));
        end
    end
end
end