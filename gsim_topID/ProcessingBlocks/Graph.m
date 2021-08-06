classdef Graph
	% Contains a graph and allows to perform several operations
	%
	
	properties
		m_adjacency % adjacency matrix
	end
	
	methods
				
		function obj = Graph(varargin)
		% Constructor	
			if nargin > 0
				obj =  assignParametersByName(obj,varargin{:});
			end	
		end
				
		function m_L = getLaplacian(obj)						
			v_degrees = sum(obj.m_adjacency,2);
			m_L = diag(v_degrees) - obj.m_adjacency;
		end
		
		function m_L = getNormalizedLaplacian(obj)						
			v_degrees = sum(obj.m_adjacency,2);						
			d_m12 = v_degrees.^(-.5);
			d_m12(v_degrees == 0 ) = 0;
			D_m12 = diag(d_m12);
			m_L = eye(length(v_degrees)) - D_m12*obj.m_adjacency*D_m12;
		    m_L = (m_L+m_L')/2;
			if sum(sum(isnan(m_L)))
				keyboard
			end
			
		end
		
		
		function m_V = getLaplacianEigenvectors(obj)									
			% m_V     Matrix whose columns are the eigenvectors of the
			%         Laplacian sorted in ascending order of eigenvalue
			
			m_L = obj.getLaplacian;
			[m_V,~] = eig(m_L);
		end
			
		function m_V = getNormalizedLaplacianEigenvectors(obj)									
			% m_V     Matrix whose columns are the eigenvectors of the
			%         Laplacian sorted in ascending order of eigenvalue
			
			m_L = obj.getNormalizedLaplacian;
			[m_V,~] = eig(m_L);

		end
		
		function s_n = getNumberOfVertices(obj)
			s_n = size(obj.m_adjacency,1);			
        end
        		
		function c_components = getComponents(obj)
			% (To be written)
			%
			% COMPONENTS is a cell array of C vectors, where C is the
			% number of components of the graph OBJ. COMPONENTS{c} is a
			% column vector containing the indices of the vertices in each
			% component. 
			%
			% The algorithm used...
            m_sparseAdjacency=sparse(obj.m_adjacency);
            [s_numberOfComponents,v_componentIndicator]=graphconncomp(m_sparseAdjacency,'DIRECTED',false,'WEAK',true);
            %v_componentIndicator is a vector that for each vertice has the
            %corresponding component number 
            %in s_numberOfComponents is the number of components
			c_components = {};
            for s_k=1:s_numberOfComponents
                c_components(:,s_k)={(find(v_componentIndicator(:)==s_k))};
            end
        end
		
		function C = getClusters(obj,s_numberOfClusters,s_Type)
            % Executes the spectral clustering algorithm defined by
            %   Type on the adjacency matrix W and returns the k cluster
            %   indicator vectors as columns in C.
            %   If L and U are also called, the (normalized) Laplacian and
            %   eigenvectors will also be returned.
            %
            % Input:
            % S_NUMBEROFCLUSTERS   Number of clusters
            %   'Type' - Defines the type of spectral clustering algorithm
            %            that should be used. Choices are:
            %      1 - Unnormalized
            %      2 - Normalized according to Shi and Malik (2000)
            %      3 - Normalized according to Jordan and Weiss (2002)
            %
            %
			% Output: 
			% C          
            %
            %   References:
            %   - Ulrike von Luxburg, "A Tutorial on Spectral Clustering",
            %     Statistics and Computing 17 (4), 2007
            %
            
            W = obj.m_adjacency;
            k = s_numberOfClusters;
            
            % calculate degree matrix
            degs = sum(W, 2);
            D    = sparse(1:size(W, 1), 1:size(W, 2), degs);
            
            % compute unnormalized Laplacian
            L = D - W;
            
            % compute normalized Laplacian if needed
            switch s_Type
                case 2
                    % avoid dividing by zero
                    degs(degs == 0) = eps;
                    % calculate inverse of D
                    D = spdiags(1./degs, 0, size(D, 1), size(D, 2));
                    
                    % calculate normalized Laplacian
                    L = D * L;
                case 3
                    % avoid dividing by zero
                    degs(degs == 0) = eps;
                    % calculate D^(-1/2)
                    D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
                    
                    % calculate normalized Laplacian
                    L = D * L * D;
            end
            
            % compute the eigenvectors corresponding to the k smallest
            % eigenvalues
            diff   = eps;
            [U, ~] = eigs(L, k, diff);
            
            % in case of the Jordan-Weiss algorithm, we need to normalize
            % the eigenvectors row-wise
            if s_Type == 3
                U = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2)));
            end
            
            % now use the k-means algorithm to cluster U row-wise
            % C will be a n-by-1 matrix containing the cluster number for
            % each data point
            C = kmeans(U, k, 'start', 'cluster', ...
                'EmptyAction', 'singleton');
            
            % now convert C to a n-by-k matrix containing the k indicator
            % vectors as columns
            C = sparse(1:size(D, 1), C, 1);
            
		end
        
		function v_f_tilde = getFourierTransform( obj , v_f )
			assert(size(v_f,1) == size(obj.m_adjacency,1));
			v_f_tilde = obj.getLaplacianEigenvectors'*v_f;
		end
		
		function plotFourierTransform( obj , v_f )
			assert(size(v_f,1) == size(obj.m_adjacency,1));
			v_f_tilde_unnormalized = obj.getLaplacianEigenvectors'*v_f;
			v_f_tilde_normalized = obj.getNormalizedLaplacianEigenvectors'*v_f;
			
			subplot(2,1,1)
			plot(abs(v_f_tilde_unnormalized)); title('Unnormalized Laplacian')
			subplot(2,1,2)
			plot(abs(v_f_tilde_normalized)); title('normalized Laplacian')

%             plot(abs([v_f_tilde_unnormalized v_f_tilde_normalized]));
% 			legend('unnormalized','normalized');

		end
		
						
		function graph = nearestNeighborsSubgraph(obj,s_neighborsNum)
			% delete all but the s_neighborsNum strongest links of each
			% edge
            W = triu(obj.m_adjacency);
            %rowIndices = 1 : size(W,1);
			
            for row = 1:size(W,1)
                % get the k-largest values and its position
                [sortedValues, sortedIndices] = sort(W(row,:),'descend');
                %maxValues = sortedValues(1:s_neighborsNum);
                minValueIndices = sortedIndices(s_neighborsNum+1:end);
                
                % update W                
                W(row, minValueIndices ) = 0;
            end
			
            %W = W .* logical(W');
			W = W + W';
            graph = Graph('m_adjacency', W);
		end
		
	end
	
	methods(Static)
		
		function G = constructViaFactorAnalysis( m_functionValues , alpha , beta )
% 			% (To be written)
% 			%
% 		    % Construct graph from signal values using [Dong et al. 2015]
% 			%
% 			% Input:
% 			% M_FUNCTIONVALUES     N x M Matrix where N is the number of
% 			%                      nodes and M is the number of
% 			%                      observations of a function on a graph
% 			% ALPHA, BETA          Regularization parameters of the
% 			%                      algorithm
% 			%
% 			% Output:
% 			% G                    Graph of N nodes
% 			%
% 			
% 			m_adjacency_est = [];
% 			G = Graph('m_adjacency',m_adjacency_est);
% 			

		end
		
		function m_adjacency= createAdjacencyFromLaplacian(m_laplacian)
            m_adjacency=eye(size(m_laplacian)).*diag(diag(m_laplacian))-m_laplacian;

		end
		
		function graph=constructGraphFromTable(m_T,ch_method)
			% m_T :            M x N matrix with NaN for missing entries
			%                  and values typically for ratings 
			% ch_method :      can be
			%    'cosine'      cosine similarity
			% 
			% graph is a graph with M vertices using a similarity measure
			% specified by ch_method
			
			
			switch ch_method
				case 'cosine'
					%Constructs a graph from a table of signals using modified
					%cosine similarity. Use only the common rated entries of the
					%vectors when computing the cosine similarity.
					m_adj=zeros(size(m_T,1));
					for k=1:size(m_T,1)
												
						for l=1:k-1
							v_xk=m_T(k,:);
							v_yk=m_T(l,:);
							v_res=v_xk.*v_yk;
							v_ind=(~isnan(v_res));
							if sum(v_ind)
								m_adj(k,l)=sum(v_res(v_ind))/(norm(v_xk(v_ind))*norm(v_yk(v_ind)));
							end
							
						end
					end
					%m_adjacency(isnan(m_adjacency)) = 0 ;
					m_adj=m_adj+m_adj';
					graph=Graph('m_adjacency',m_adj);										
					
				otherwise
					error('unrecognized option');
			end
		end
		
	end
	
	
end

