function [v_NMSD, v_avgEIER, m_VARcoeff_permuteReshaped, t_A_estimates] = ...
    topologySimulation(ch_algorithm, ...
    noOfNodes, filterOrder, noOfMCitr, noOfTimeInstants, ...
    generator, erdosRenyiEdgeProbability, regPar, threshold, ...
    forgettingFactor, niterPGD, stepsizeFactor_TISO)

N = noOfNodes;
P = filterOrder;
M = noOfMCitr; % indexed by m_exp
T = noOfTimeInstants;
p_edge_ER = erdosRenyiEdgeProbability;

gen_clone = VARSTFunctionGenerator;
gen_clone.sigma = generator.sigma;
gen_clone.nTimeSamples = T;

mu= 1-forgettingFactor;

m_VARcoeff_permuteReshaped = zeros(N,N*P,M);
t_A_estimates = zeros(N,N*P,T,M);
m_EIER = zeros(T,M);
m_numeratorErrors = zeros(T,M);
v_denominatorTerms = zeros(M,1);

rng(0); % random seed
ltc=LoopTimeControl(M);
for m_exp=1:M
    % generating VAR time series data:
    m_adjacency =or(rand(N)<p_edge_ER,eye(N));
    myGraph = Graph('m_adjacency',m_adjacency);
    t_VARCoefficients = VARSTFunctionGenerator.randomCoefficientsFromGraph(myGraph, P);
    gen_clone.t_coefficients = t_VARCoefficients;
    m_data = gen_clone.realization();
    trueWeightedAdjacency=VARSTFunctionGenerator.VARCoefficients2AdjacencyMatrix(t_VARCoefficients);   

    m_A_initial= zeros(N, N*P);
    m_VARcoeff_permuteReshaped(:,:,m_exp)= reshape(permute(t_VARCoefficients, [1,3,2]),N,N*P);
%     %%% regPar selection%%%%%%%
%     if exp==1
%         noOfRegPars= 20;
%         regPars= logspace(0.01,-4,noOfRegPars);
%         t_A_TISOestimatesRegPars=zeros(noOfNodes,noOfNodes*filterOrder,noOfTimeInstants,noOfRegPars);
%         num_errorTIRSORegPars= zeros(noOfTimeInstants,noOfRegPars);
%         errorSumTISORegPars= zeros(noOfRegPars,1);
%         lipschitzpar= max((norms(m_data,2)).^2);
%         v_stepSizeTISO= 1./(lipschitzpar.*ones(noOfTimeInstants,1));
%         %[~,v_LipschitzParTIRSO1]= PGDTIRSO (filterOrder, 0, forgettingFactor, mu, sigma, zeros(noOfTimeInstants,1), m_data, m_A_initial,0);
%         %v_stepSizeTIRSOregPar= 1./(max(v_LipschitzParTIRSO1).*ones(noOfTimeInstants,1));
%         %[~,v_strcvxparTIRSO]= TIRSO (filterOrder, 0, forgettingFactor, mu, sigma, zeros(noOfTimeInstants,1), m_data, m_A_initial,0,0);
%         %v_stepSizeTIRSOregPar= 0.01./(min(v_strcvxparTIRSO(noOfTimeInstants/2:end)).*(1:noOfTimeInstants));
%         for k=1:noOfRegPars
%             t_A_TISOestimatesRegPars(:,:,:,k)= TISO (filterOrder, regPars(k), v_stepSizeTISO, m_data, m_A_initial);
%             for t= filterOrder+1:noOfTimeInstants
%                 num_errorTIRSORegPars(t,k)= sum((norms((t_A_TISOestimatesRegPars(:,:,t,k)- m_VARcoeff_permuteReshaped(:,:,exp)),2,2)).^2);
%             end
%             errorSumTISORegPars(k)= sum(num_errorTIRSORegPars(:,k));
%         end
%         [~,minIndex]=min(errorSumTISORegPars);
%         regPar= regPars(minIndex);
%     end
%     
%     figure
%     semilogy(1:noOfRegPars,errorSumTIRSORegPars)
%     %%% end regpar selection %%%%%
                
               
% computing Lipschitz smoothness parameter and running algorithm
    switch ch_algorithm
        case 'TISO'
            lipschitzpar= max((norms(m_data,2)).^2);
            v_stepSizeTISO= stepsizeFactor_TISO.*ones(T,1)./lipschitzpar;
            t_A_estimates(:,:,:,m_exp)= TISO (...
                P, regPar, v_stepSizeTISO, m_data, m_A_initial);
        case 'OSGDTISO'
            lipschitzpar_loss= max((norms(m_data,2)).^2);
            lipschitzpar = lipschitzpar_loss; 
            % this is actually Lipschitz continuity param: regPars(regparIndex)*sqrt(N);
            v_stepSizeTISO= stepsizeFactor_TISO.*ones(T,1)./lipschitzpar;
            t_A_estimates(:,:,:,m_exp)= OSGDTISO (...
                P, regPar, v_stepSizeTISO, m_data, m_A_initial);
            % OBS! do we need to increase the lipschitz par in
            % OSGDTISO? The cost function includes now the
            % regularizer
        case 'TIRSO'
            [~,v_LipschitzParTIRSO]= PGDTIRSO (...
                P, regPar, forgettingFactor, mu, 0, zeros(T,1), m_data, m_A_initial,0);
            v_stepSizeTIRSO= 1./(max(v_LipschitzParTIRSO).*ones(T,1));
            % based on Lipschitz smoothness par
            t_A_estimates(:,:,:,m_exp)= TIRSO (...
                P, regPar, forgettingFactor, mu, 0, v_stepSizeTIRSO, m_data, m_A_initial,0,0);
        case 'PGDTIRSO'
            [~,v_LipschitzParTIRSO]= PGDTIRSO (...
                P, regPar, forgettingFactor, mu, 0, zeros(T,1), m_data, m_A_initial,0);
            v_stepSizeTIRSO= 1./(max(v_LipschitzParTIRSO).*ones(T,1));
            t_A_estimates(:,:,:,m_exp)= PGDTIRSO (...
                P, regPar, forgettingFactor, mu, 0, v_stepSizeTIRSO, m_data, m_A_initial, niterPGD);
    end
                
    for time= P+1:T
        % EIER
        m_A_exp=squeeze(t_A_estimates(:,:,time,m_exp));
        Aux= permute(m_A_exp,[2 1 3]);
        Aux1= reshape(Aux, [P,N,N]);
        Aux2= permute(Aux1,[3,2,1]);
        m_inferredAdjacency= norms(Aux2,2,3);
        [~,~,~,~,m_EIER(time,m_exp)]= ...
            TISOTIRSOExperiments.computeEierMissdetectionFalsealarm(...
            threshold,N,trueWeightedAdjacency, m_inferredAdjacency);
        % NMSD
        m_numeratorErrors(time,m_exp)= sum( ...
            (norms((t_A_estimates(:,:,time,m_exp)- m_VARcoeff_permuteReshaped(:,:,m_exp)),2,2)).^2 ...
            );
    end
    v_denominatorTerms(m_exp)= sum((norms(m_VARcoeff_permuteReshaped(:,:,m_exp),2,2)).^2);
    
    ltc.go(m_exp);
end

summed_den = sum(v_denominatorTerms);
v_NMSD = sum(m_numeratorErrors,2)./summed_den;
v_avgEIER = mean(m_EIER,2);

