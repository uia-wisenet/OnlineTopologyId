function  [lambda_xvalNMSD,lambda_xvalEIER] = ...
    toplogySimulationXvalNMSDEIER(ch_algorithm, ...
    noOfNodes, filterOrder, noOfTimeInstants, regParRange, noOfRegPars,...
    generator, erdosRenyiEdgeProbability, threshold, ...
    forgettingFactor, niterPGD, stepsizeFactor_TISO)
% Notation for convenience
N = noOfNodes;
P = filterOrder;
T = noOfTimeInstants;
p_edge_ER = erdosRenyiEdgeProbability;
mu= 1-forgettingFactor;
% data generator
gen_clone = VARSTFunctionGenerator;
gen_clone.sigma = generator.sigma;
gen_clone.nTimeSamples = T;
% declaring variables
t_A_estimates = zeros(N,N*P,T);
m_EIERregPars = zeros(T,noOfRegPars);
num_NMSDRegPars= zeros(T,noOfRegPars);
NMSDSumRegPars= zeros(noOfRegPars,1);
EIERSumRegPars= zeros(noOfRegPars,1);

rng(0); % random seed
% generating VAR time series data:
m_adjacency =or(rand(N)<p_edge_ER,eye(N));
myGraph = Graph('m_adjacency',m_adjacency);
t_VARCoefficients = VARSTFunctionGenerator.randomCoefficientsFromGraph(myGraph, P);
gen_clone.t_coefficients = t_VARCoefficients;
m_data = gen_clone.realization();
trueWeightedAdjacency=VARSTFunctionGenerator.VARCoefficients2AdjacencyMatrix(t_VARCoefficients);

m_A_initial= zeros(N, N*P);
m_VARcoeff_permuteReshaped= reshape(permute(t_VARCoefficients, [1,3,2]),N,N*P);
regPars= logspace(regParRange(1),regParRange(2),noOfRegPars);
% step size for TISO
lipschitzpar= max((norms(m_data,2)).^2);
v_stepSizeTISO= stepsizeFactor_TISO.*ones(T,1)./lipschitzpar;
% step size for OSGDTIS is computed inside the switch
% step size for TIRSO and PGDTIRSO
[~,v_LipschitzParTIRSO]= PGDTIRSO (...
    P, 0, forgettingFactor, mu, 0, zeros(T,1), m_data, m_A_initial,0);
v_stepSizeTIRSO= 1./(max(v_LipschitzParTIRSO).*ones(T,1));
ltc=LoopTimeControl(noOfRegPars);
for regparIndex=1:noOfRegPars
    % computing Lipschitz smoothness parameter and running algorithm
    switch ch_algorithm
        case 'TISO'
            t_A_estimates(:,:,:,regparIndex)= TISO (...
                P, regPars(regparIndex), v_stepSizeTISO, m_data, m_A_initial);
        case 'OSGDTISO'
            % step size for OSGDTISO
            lipschitzpar_loss= max((norms(m_data,2)).^2);
            lipschitzparOSGD = lipschitzpar_loss;
            % this is actually Lipschitz continuity param: regPars(regparIndex)*sqrt(N);
            v_stepSizeOSGDTISO= stepsizeFactor_TISO.*ones(T,1)./lipschitzparOSGD;
            t_A_estimates(:,:,:,regparIndex)= OSGDTISO (...
                P, regPars(regparIndex), v_stepSizeOSGDTISO, m_data, m_A_initial);
            % OBS! do we need to increase the lipschitz par in
            % OSGDTISO? The cost function includes now the
            % regularizer
        case 'TIRSO'
            t_A_estimates(:,:,:,regparIndex)= TIRSO (...
                P, regPars(regparIndex), forgettingFactor, mu, 0, v_stepSizeTIRSO, m_data, m_A_initial,0,0);
        case 'PGDTIRSO'
            t_A_estimates(:,:,:,regparIndex)= PGDTIRSO (...
                P, regPars(regparIndex), forgettingFactor, mu, 0, v_stepSizeTIRSO, m_data, m_A_initial, niterPGD);
    end
    for t= filterOrder+1:noOfTimeInstants
        num_NMSDRegPars(t,regparIndex)= sum((norms((t_A_estimates(:,:,t,regparIndex)- m_VARcoeff_permuteReshaped),2,2)).^2);
    end
    NMSDSumRegPars(regparIndex)= sum(num_NMSDRegPars(:,regparIndex));
    for time= P+1:T
        % EIER
        m_A_exp=squeeze(t_A_estimates(:,:,time,regparIndex));
        Aux= permute(m_A_exp,[2 1 3]);
        Aux1= reshape(Aux, [P,N,N]);
        Aux2= permute(Aux1,[3,2,1]);
        m_inferredAdjacency= norms(Aux2,2,3);
        [~,~,~,~,m_EIERregPars(time,regparIndex)]= ...
            TISOTIRSOExperiments.computeEierMissdetectionFalsealarm(...
            threshold,N,trueWeightedAdjacency, m_inferredAdjacency);
    end
    EIERSumRegPars(regparIndex)= sum(m_EIERregPars(:,regparIndex));
    ltc.go(regparIndex);
end
% computing the value of regPar with minimum NMSD
[~,minIndexNMSD]=min(NMSDSumRegPars);
lambda_xvalNMSD= regPars(minIndexNMSD);
% computing the value of regPar with minimum EIER
[~,minIndexEIER]=min(EIERSumRegPars);
lambda_xvalEIER= regPars(minIndexEIER);
% Plotting both the curves
xaxis=1:noOfRegPars;
semilogy(xaxis,NMSDSumRegPars)
hold on
semilogy(xaxis,EIERSumRegPars)
legend('NMSD','EIER')
end