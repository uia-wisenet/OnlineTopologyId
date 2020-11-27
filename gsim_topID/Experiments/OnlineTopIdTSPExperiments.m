classdef OnlineTopIdTSPExperiments < ExperimentFunctionSet
    % This class contains the experiments for the numerical results in
    % B. Zaman, L. M. Lopez-Ramos, D. Romero, B. Beferull-Lozano, 
    % "Online Topology Identification from Vector Autoregressive
    % Time-series," To appear in IEEE Transactions on Signal Processing.
    methods
        % Figure 1: not an experiment in the paper ( Fig. Tensor of VAR parameters)
        
        % Figure 2
        function F=experiment_2(obj, niter)
            % Stationary VAR process, experiments with multiple
            % regularization parameters to plot NMSD, NMSE, ROC, P_FA,
            % P_MD, and EIER for TISO and TIRSO.
            % For ROC, P_FA and P_MD are averaged over a given time interval
            % For P_FA, P_MD, and EIER, the threshold is selected based on
            % minimum P_MD - P_FA (Equal error rate point where P_MD = P_FA)
            % The probabilities are first computed for an interval and
            % then are averaged across MC experiments
            % Adaptive, non-diminishing stepsize
            %% Simulation parameters
            noOfNodes        = 12;
            filtOrder        = 2;
            noOfExperiments  = 50; 
            noOfObservations = 3000;
            EdgeProbability  = 0.2;
            syntheticSigma   = 0.005;
            initialSigma     = 0.01;
            f_factor         = 0.98;
            noOfThresholds   = 1000;
            v_thresholds = linspace(1e-8, 10e-1, noOfThresholds);
            v_regPars    = [1e-7 1e-6 1e-5];
            n_regPars    = length(v_regPars);
            %% Defining the variables
            startingPoint = 10; % t1, interval=T-t1
            interval      = noOfObservations- startingPoint;
            roc_MDtirsoAllExp =zeros(n_regPars, noOfThresholds,noOfExperiments,interval);
            roc_FAtirsoAllExp =zeros(n_regPars, noOfThresholds,noOfExperiments,interval);
            roc_MDtisoAllExp  =zeros(n_regPars, noOfThresholds,noOfExperiments,interval);
            roc_FAtisoAllExp  =zeros(n_regPars, noOfThresholds,noOfExperiments,interval);
            EIER_tirso    = zeros(n_regPars, noOfThresholds,noOfExperiments,interval);
            EIER_tiso     = zeros(n_regPars, noOfThresholds,noOfExperiments,interval);
            power         = zeros (noOfExperiments, noOfObservations);
            grelos        = cell  (1, n_regPars);
            glos          = cell  (1, n_regPars);
            t_Atirso      = zeros(noOfNodes,noOfNodes,noOfExperiments,noOfObservations,n_regPars);
            t_Atiso       = zeros(noOfNodes,noOfNodes,noOfExperiments,noOfObservations,n_regPars);
            m_adjacency   = zeros(noOfNodes, noOfNodes, noOfExperiments);
            myCoefficients= zeros(noOfNodes,noOfNodes,filtOrder,noOfExperiments);
            %ltc = LoopTimeControl(noOfExperiments);
            %% Loop for the Monte Carlo iterations
            parfor exp=1:noOfExperiments
                % Synthetic Data Generator
                %rng(3);
                m_adjacency(:,:,exp) = (rand(noOfNodes)<EdgeProbability).*not(eye(noOfNodes));
                myGraph = Graph('m_adjacency',m_adjacency(:,:,exp));
                myCoefficients1 = VARSTFunctionGenerator.randomCoefficientsFromGraph(myGraph, filtOrder);
                gen = VARSTFunctionGenerator;
                gen.nTimeSamples = noOfObservations;
                gen.sigma = syntheticSigma;
                myCoefficients(:,:,:,exp) = myCoefficients1;
                gen.t_coefficients = myCoefficients(:,:,:,exp);
                [sq_errorTirso(exp,:,:),sq_errorTiso(exp,:,:),predict_errorTirso(exp,:,:),...
                    predict_errorTiso(exp,:,:),power(exp,:),t_Atirso(:,:,exp,:,:),t_Atiso(:,:,exp,:,:)]= ...
                    computationsForTirsoTiso(noOfNodes, n_regPars,v_regPars, f_factor, gen, grelos, ...
                    glos, filtOrder, noOfObservations, squeeze(myCoefficients(:,:,:,exp)),exp,initialSigma)
                
                % Genie-aided
                [~, predict_errorGenie(exp,:), ~] =...
                    OMVARSimulator(gen, 'genie', filtOrder, noOfObservations, myCoefficients(:,:,:,exp));
                
                %ltc.go(exp)
            end
            % Computing the P_FA and P_MD for ROC in both tirso and tiso
            for exp=1:noOfExperiments
                for i_r=1:n_regPars
                    for t= 1:interval % for ROC, only interval is considered
                        for j_thr=1:noOfThresholds
                            [roc_MDtirsoAllExp(i_r,j_thr,exp,t),roc_FAtirsoAllExp(i_r,j_thr,exp,t),EIER_tirso(i_r,j_thr,exp,t)]= ...
                                obj.computeEierMissdetectionFalsealarm(v_thresholds(j_thr),...
                                noOfNodes, m_adjacency(:,:,exp), t_Atirso(:,:,exp,startingPoint+t, i_r));
                            [roc_MDtisoAllExp(i_r,j_thr,exp,t),roc_FAtisoAllExp(i_r,j_thr,exp,t),EIER_tiso(i_r,j_thr,exp,t)]= ...
                                obj.computeEierMissdetectionFalsealarm(v_thresholds(j_thr),...
                                noOfNodes, m_adjacency(:,:,exp), t_Atiso(:,:,exp,startingPoint+t, i_r));
                        end
                    end
                end
            end
            % averaging across experiments and time for ROC and for minimum P_MD-P_FA
            roc_MDtirso =mean(mean(roc_MDtirsoAllExp,4),3);
            roc_FAtirso =mean(mean(roc_FAtirsoAllExp,4),3);
            roc_MDtiso  =mean(mean(roc_MDtisoAllExp,4),3);
            roc_FAtiso  =mean(mean(roc_FAtisoAllExp,4),3);
            [minMDMinusFAtirso_values,minMDMinusFAtirso_thrIndices]=min(abs(roc_MDtirso-roc_FAtirso),[],2);
            [minMDMinusFAtiso_values,minMDMinusFAtiso_thrIndices]  =min(abs(roc_MDtiso-roc_FAtiso),[],2);
            
            % Initializing variables to plot the P_MD, P_FA,
            % EIER across time for both TIRSO and TISO
            MD_tirso =      zeros(n_regPars,noOfObservations,noOfExperiments);
            FA_tirso =      zeros(n_regPars,noOfObservations,noOfExperiments);
            minEIER_tirso = zeros(n_regPars,noOfObservations,noOfExperiments);
            MD_tiso =       zeros(n_regPars,noOfObservations,noOfExperiments);
            FA_tiso =       zeros(n_regPars,noOfObservations,noOfExperiments);
            minEIER_tiso =  zeros(n_regPars,noOfObservations,noOfExperiments);
            for exp=1:noOfExperiments
                for obs=5:noOfObservations
                    for m=1:n_regPars
                        [MD_tirso(m,obs,exp), FA_tirso(m,obs,exp),minEIER_tirso(m,obs,exp)]= ...
                            obj.computeEierMissdetectionFalsealarm...
                            (v_thresholds(minMDMinusFAtirso_thrIndices(m)),noOfNodes,m_adjacency(:,:,exp),t_Atirso(:,:,exp,obs,m));
                        [MD_tiso(m,obs,exp), FA_tiso(m,obs,exp),minEIER_tiso(m,obs,exp)]= ...
                            obj.computeEierMissdetectionFalsealarm...
                            (v_thresholds(minMDMinusFAtiso_thrIndices(m)),noOfNodes,m_adjacency(:,:,exp),t_Atiso(:,:,exp,obs,m));
                    end
                end
            end
            avgminMD_tirso=   mean(MD_tirso,3);
            avgminFA_tirso=   mean(FA_tirso,3);
            avgminEIER_tirso= mean(minEIER_tirso,3);
            avgminMD_tiso=    mean(MD_tiso,3);
            avgminFA_tiso=    mean(FA_tiso,3);
            avgminEIER_tiso=  mean(minEIER_tiso,3);
            
            %% Calculating NMSD and plotting
            sq_error= [transpose(squeeze(mean(sq_errorTirso,1)));transpose(squeeze(mean(sq_errorTiso,1)))];
            m_NMSD = zeros(2*n_regPars, size(sq_error, 2));
            normOfCoeff= ((norm(vec(myCoefficients)))^2)/noOfExperiments;
            for i = 1:2*n_regPars
                devNumerator=sq_error(i,:);
                m_NMSD(i,:)=devNumerator./normOfCoeff;
            end
            predict_error= [(squeeze(mean(predict_errorTirso,1)))'; (squeeze(mean(predict_errorTiso,1)))';mean(predict_errorGenie,1)];
            m_NMSE = zeros(2*n_regPars+1, size(predict_error,2));
            m_filteredNMSE = zeros(2*n_regPars+1, size(predict_error, 2));
            avgLength = 100;
            for k = 1:2*n_regPars+1
                m_NMSE(k,:)=predict_error(k,:) ./ mean( power,1);
                m_filteredNMSE(k,:) = filter(ones(1, avgLength)/avgLength, 1, m_NMSE(k,:));
            end
            
            M(1,1)=GFigure();
            M(1,1).b_logy=1;
            for i = 1:n_regPars
                M(1,1).c_legend{i} = sprintf(...
                    'TIRSO, \\lambda=%0.5g', v_regPars(i));
                M(1,1).c_styles{i} = '-';
                M(1,1).c_legend{n_regPars+i} = sprintf(...
                    'TISO, \\lambda=%0.5g',  v_regPars(i));
                M(1,1).c_styles{n_regPars+i} = '--';
            end
            M(1,1).colorPeriod = n_regPars;
            M(1,1).m_Y= m_NMSD;
            M(1,1).ch_ylabel='NMSD';
            M(1,1).ch_xlabel='time';
            
            M(1,2)=GFigure();
            M(1,2).b_logy=1;
            M(1,2).c_legend=M(1,1).c_legend;
            M(1,2).c_legend{2*n_regPars+1} = 'Genie';
            M(1,2).colorPeriod = n_regPars;
            M(1,2).c_styles = M(1,1).c_styles;
            M(1,2).c_styles{2*n_regPars+1} = '-.';
            M(1,2).m_Y= m_filteredNMSE;
            %M(1,2).m_Y= m_NMSE;
            M(1,2).ch_ylabel='NMSE';
            M(1,2).ch_xlabel='time';
            
            M(1,3)=GFigure();
            M(1,3).c_legend= M(1,1).c_legend;
            M(1,3).colorPeriod = n_regPars;
            M(1,3).c_styles = M(1,1).c_styles;
            M(1,3).m_Y=[roc_MDtirso;roc_MDtiso];
            M(1,3). m_X=[roc_FAtirso; roc_FAtiso];
            M(1,3).ch_ylabel='P_{MD}';
            M(1,3).ch_xlabel='P_{FA}';
            
            M(2,1) = GFigure();
            M(2,1).b_logy=1;
            M(2,1).c_legend=M(1,1).c_legend;
            M(2,1).colorPeriod = n_regPars;
            M(2,1).c_styles = M(1,1).c_styles;
            %M(2,1).m_Y = m_EIER;
            M(2,1).m_Y = [avgminEIER_tirso;avgminEIER_tiso];
            M(2,1).ch_ylabel = 'EIER';
            M(2,1).ch_xlabel = 'time';
            
            M(2,2) = GFigure();
            M(2,2).b_logy=1;
            M(2,2).c_legend=M(1,1).c_legend;
            M(2,2).colorPeriod = n_regPars;
            M(2,2).c_styles = M(1,1).c_styles;
            M(2,2).m_Y = [avgminMD_tirso;avgminMD_tiso];
            M(2,2).ch_ylabel = 'P_{MD}';
            M(2,2).ch_xlabel = 'time';
            
            M(2,3) = GFigure();
            M(2,3).b_logy=1;
            M(2,3).c_legend=M(1,1).c_legend;
            M(2,3).colorPeriod = n_regPars;
            M(2,3).c_styles = M(1,1).c_styles;
            M(2,3).m_Y = [avgminFA_tirso; avgminFA_tiso];
            M(2,3).ch_ylabel = 'P_{FA}';
            M(2,3).ch_xlabel = 'time';
            
            F=GFigure('m_multiplot',M);
            
        end
        
        % Figure 3
        function F= experiment_3(obj, niter)
            % Evaluating NMSD vs time and comparing it for TIRSO and TISO
            % with different choices of step sizes
            % Stationary settings
            noOfNodes        = 10;
            filterOrder      = 3;
            noOfTimeInstants = 2000;
            noOfMCitr        = 50; %Monte Carlo iterations
            erdosRenyiEdgeProbability = 0.2;
            noiseSigma       = 0.1;
            regPar           = 8e-4;
            forgettingFactor = 0.99;
            initial_sigma    = 0.01;
            mu = 1-forgettingFactor;
            
            % Initialization
            num_NMSDTIRSOLf           =zeros(noOfTimeInstants,noOfMCitr);
            num_NMSDTIRSOLfadaptive   =zeros(noOfTimeInstants,noOfMCitr);
            num_NMSDTIRSOLfsqrttdim   =zeros(noOfTimeInstants,noOfMCitr);
            num_NMSDTIRSOStrConvParTimest=zeros(noOfTimeInstants,noOfMCitr);
            num_NMSDTIRSOtrace        =zeros(noOfTimeInstants,noOfMCitr);
            num_NMSDTIRSOtracedim     =zeros(noOfTimeInstants,noOfMCitr);
            num_NMSDTISOLfdimRunningMax=zeros(noOfTimeInstants,noOfMCitr);
            num_NMSDTISOadaptdimTrace =zeros(noOfTimeInstants,noOfMCitr);
            num_NMSDTISOdimL          =zeros(noOfTimeInstants,noOfMCitr);
            num_NMSDTISOconstantL     =zeros(noOfTimeInstants,noOfMCitr);
            denNMSD                   =zeros(noOfTimeInstants,noOfMCitr);
            
            %% Loop for the Monte Carlo iterations
            ltc = LoopTimeControl(noOfMCitr);
            rng(0);
            for exp=1:noOfMCitr
                %% Synthetic Data Generation
                m_adjacency =or(rand(noOfNodes)<erdosRenyiEdgeProbability,eye(noOfNodes));
                myGraph = Graph('m_adjacency',m_adjacency);
                t_VARCoefficients = VARSTFunctionGenerator.randomCoefficientsFromGraph(myGraph, filterOrder);
                gen = VARSTFunctionGenerator;
                gen.nTimeSamples = noOfTimeInstants;
                gen.sigma = noiseSigma;
                gen.t_coefficients = t_VARCoefficients;
                m_data= gen.realization;
                %% TIRSO
                m_A_initial= zeros(noOfNodes, noOfNodes*filterOrder);
                [~,v_LipschitzParTIRSO]= PGDTIRSO (filterOrder, regPar, forgettingFactor, mu, initial_sigma, zeros(noOfTimeInstants,1), m_data, m_A_initial,0);
                % step size = 1/L_f, where L_f is the Lipschitz smoothness
                % parameter of the loss function f for all t
                v_constantStepSizeTIRSO= 1./(max(v_LipschitzParTIRSO).*ones(noOfTimeInstants,1)); % based on LIpschitz smoothness par
                t_A_TIRSOLf= TIRSO (filterOrder, regPar, forgettingFactor, mu, initial_sigma, v_constantStepSizeTIRSO, m_data, m_A_initial,0,0);
                % step size= 1/(L_f(t)), instantaneous L parameter
                v_adaptiveStepSizeTIRSO= 1./(v_LipschitzParTIRSO); %adaptive
                t_A_TIRSOLfadaptive= TIRSO (filterOrder, regPar, forgettingFactor, mu, initial_sigma, v_adaptiveStepSizeTIRSO, m_data, m_A_initial,0,0);
                % step size= 1/(L_f*sqrt(t))
                v_dimLStepSizeTIRSO= 1./(max(v_LipschitzParTIRSO).*sqrt((1:noOfTimeInstants))); % 1/(L_f.sqrt(t)
                t_A_TIRSOLfsqrttdim= TIRSO (filterOrder, regPar, forgettingFactor, mu, initial_sigma, v_dimLStepSizeTIRSO, m_data, m_A_initial,0,0);
                % step size = 1/(H_f*t) where H_f is the strong convexity
                % parameter of the loss function f for all t
                [~,v_strcvxparTIRSO]= TIRSO (filterOrder, regPar, forgettingFactor, mu, initial_sigma, zeros(noOfTimeInstants,1), m_data, m_A_initial,0,0);
                v_stepSizeTIRSO_O_1_over_t= min( 1./v_LipschitzParTIRSO, ...
                    1./(max(v_strcvxparTIRSO).*(1:noOfTimeInstants)') );
                t_A_TIRSOStrConvParTimest= TIRSO (filterOrder, regPar, forgettingFactor, mu, initial_sigma, v_stepSizeTIRSO_O_1_over_t, m_data, m_A_initial,0,0);
                % step size = 1/(trace(Phi_t))
                b_trace=1;
                b_diminishing= 0;
                v_adaptiveStepSizeTIRSO= zeros(noOfTimeInstants,1); % assigned inside TIRSO function
                t_A_TIRSOtrace= TIRSO (filterOrder, regPar, forgettingFactor, mu, initial_sigma, v_adaptiveStepSizeTIRSO, m_data, m_A_initial,b_trace,b_diminishing);
                % step size = 1/(trace(Phi_t)*sqrt(t))
                b_trace=1;
                b_diminishing= 1;
                t_A_TIRSOtracedim= TIRSO (filterOrder, regPar, forgettingFactor, mu, initial_sigma, v_adaptiveStepSizeTIRSO, m_data, m_A_initial,b_trace,b_diminishing);
                %% TISO
                % step size = 1/L_f, where L_f is the Lipschitz smoothness parameter
                v_lipschitzpar= (norms(m_data,2)).^2;
                v_runningMaxStepTISO_old= zeros(1,noOfTimeInstants);
                for par=filterOrder+1:length(v_lipschitzpar)
                    v_runningMaxStepTISO_old(par)= max(v_runningMaxStepTISO_old(par-1),v_lipschitzpar(par));
                end
                v_runningMaxStepTISO = cummax(v_runningMaxStepTISO_old);
                max_lipschitzpar= max(v_lipschitzpar);
                % step size = 1/L_f
                v_constantLstepSizeTISO= 1./(max_lipschitzpar.*ones(1,noOfTimeInstants));
                t_A_constantLTISOestimates= TISO (filterOrder, regPar, v_constantLstepSizeTISO, m_data, m_A_initial);
                % step size = 1/(L_f*sqrt(t))
                v_dimLstepSizeTISO= 1./(max_lipschitzpar.*sqrt((1:noOfTimeInstants)));
                t_A_dimLTISOestimates= TISO (filterOrder, regPar, v_dimLstepSizeTISO, m_data, m_A_initial);
                % step size = 1/(L_f(t)*sqrt(t))
                v_adaptdimTraceStepSizeTISO= 1./(v_lipschitzpar.*sqrt((1:noOfTimeInstants)));
                t_A_adaptdimTraceTISOestimates= TISO (filterOrder, regPar, v_adaptdimTraceStepSizeTISO, m_data, m_A_initial);
                % step size = running max of L_f(t)
                v_dimRunningMaxStepSizeTISO= 1./(v_runningMaxStepTISO.*sqrt(1:noOfTimeInstants));
                t_A_dimRunningMaxTISOestimates= TISO (filterOrder, regPar, v_dimRunningMaxStepSizeTISO, m_data, m_A_initial);
                m_VARcoeff_permuteReshaped= reshape(permute(t_VARCoefficients, [1,3,2]),noOfNodes,noOfNodes*filterOrder);
                
                %% computing NMSD numerator and denominator for all variants
                for time= filterOrder+1:noOfTimeInstants-1
                    num_NMSDTIRSOLf(time,exp)= ...
                        sum((norms((t_A_TIRSOLf(:,:,time)- m_VARcoeff_permuteReshaped),2,2)).^2);
                    num_NMSDTIRSOLfadaptive(time,exp)= ...
                        sum((norms((t_A_TIRSOLfadaptive(:,:,time)- m_VARcoeff_permuteReshaped),2,2)).^2);
                    num_NMSDTIRSOLfsqrttdim(time,exp)= ...
                        sum((norms((t_A_TIRSOLfsqrttdim(:,:,time)- m_VARcoeff_permuteReshaped),2,2)).^2);
                    num_NMSDTIRSOStrConvParTimest(time,exp)= ...
                        sum((norms((t_A_TIRSOStrConvParTimest(:,:,time)- m_VARcoeff_permuteReshaped),2,2)).^2);
                    num_NMSDTIRSOtrace(time,exp)= ...
                        sum((norms((t_A_TIRSOtrace(:,:,time)- m_VARcoeff_permuteReshaped),2,2)).^2);
                    num_NMSDTIRSOtracedim(time,exp)= ...
                        sum((norms((t_A_TIRSOtracedim(:,:,time)- m_VARcoeff_permuteReshaped),2,2)).^2);
                    num_NMSDTISOconstantL(time,exp)= ...
                        sum((norms((t_A_constantLTISOestimates(:,:,time)- m_VARcoeff_permuteReshaped),2,2)).^2);
                    num_NMSDTISOdimL(time,exp)= ...
                        sum((norms((t_A_dimLTISOestimates(:,:,time)- m_VARcoeff_permuteReshaped),2,2)).^2);
                    num_NMSDTISOadaptdimTrace(time,exp)= ...
                        sum((norms((t_A_adaptdimTraceTISOestimates(:,:,time)- m_VARcoeff_permuteReshaped),2,2)).^2);
                    num_NMSDTISOLfdimRunningMax(time,exp)= ...
                        sum((norms((t_A_dimRunningMaxTISOestimates(:,:,time)- m_VARcoeff_permuteReshaped),2,2)).^2);
                    denNMSD(time,exp)= ...
                        sum((norms(m_VARcoeff_permuteReshaped,2,2)).^2);
                end
                ltc.go(exp);
            end
            
            % averaging over experiments to obtain error measure
            avg_denNMSD=mean(denNMSD,2);
            errorNMSDTIRSOLf= mean(num_NMSDTIRSOLf,2)./avg_denNMSD;
            errorNMSDTIRSOLfadaptive= mean(num_NMSDTIRSOLfadaptive,2)./avg_denNMSD;
            errorNMSDTIRSOLfsqrttdim= mean(num_NMSDTIRSOLfsqrttdim,2)./avg_denNMSD;
            errorNMSDTIRSOStrConvParTimest= mean(num_NMSDTIRSOStrConvParTimest,2)./avg_denNMSD;
            errorNMSDTIRSOtrace= mean(num_NMSDTIRSOtrace,2)./avg_denNMSD;
            errorNMSDTIRSOtracedim= mean(num_NMSDTIRSOtracedim,2)./avg_denNMSD;
            errorNMSDTISOconstantL= mean(num_NMSDTISOconstantL,2)./avg_denNMSD;
            errorNMSDTISOdimL= mean(num_NMSDTISOdimL,2)./avg_denNMSD;
            errorNMSDTISOadaptdimTrace= mean(num_NMSDTISOadaptdimTrace,2)./avg_denNMSD;
            errorNMSDTISOLfdimRunningMax= mean(num_NMSDTISOLfdimRunningMax,2)./avg_denNMSD;
            
            %% Displaying figures
            figure(1); clf
            timeAxis=filterOrder+1:noOfTimeInstants-1;
            indicesToDisplay= filterOrder+1:200:noOfTimeInstants-1;
            semilogy(timeAxis,errorNMSDTIRSOtracedim(filterOrder+1:end-1),'-or',...
                timeAxis,errorNMSDTIRSOtrace(filterOrder+1:end-1),'-+b',...
                timeAxis,errorNMSDTIRSOStrConvParTimest(filterOrder+1:end-1),'-xg',...
                timeAxis,errorNMSDTIRSOLfsqrttdim(filterOrder+1:end-1),'- square b',...
                timeAxis,errorNMSDTIRSOLfadaptive(filterOrder+1:end-1),'- diamond m', ...
                timeAxis,errorNMSDTIRSOLf(filterOrder+1:end-1),'-^c',...
                timeAxis,errorNMSDTISOconstantL(filterOrder+1:end-1),'--^c',...
                timeAxis,errorNMSDTISOdimL(filterOrder+1:end-1),'-- square r',...
                timeAxis,errorNMSDTISOadaptdimTrace(filterOrder+1:end-1),'-- > b',...
                timeAxis,errorNMSDTISOLfdimRunningMax(filterOrder+1:end-1), '-- < g', ...
                'MarkerIndices',indicesToDisplay, 'LineWidth', 1.5)
            leg1=legend({'TIRSO $\alpha_t = 1/(\mathrm{trace}( \mathbf \Phi[t])\sqrt{t})$', ...
                'TIRSO $\alpha_t=1/(\mathrm{trace}(\mathbf \Phi[t]))$',...
                'TIRSO $\alpha_t=1/(\beta_{\tilde \ell}\cdot t)$', ...
                'TIRSO $\alpha_t=1/(L\sqrt{t})$',...
                'TIRSO $\alpha_t=1/L(t)$',...
                'TIRSO $\alpha_t=1/L$',...
                'TISO $\alpha_t=1/L$',...
                'TISO $\alpha_t=1/(L\sqrt{t})$',...
                'TISO $\alpha_t=1/(L(t)\sqrt{t})$',...
                'TISO $\alpha_t=1/(L_{\mathrm{max}}(t)\cdot \sqrt{t})$'});
            set(leg1,'Interpreter','latex');
            set(leg1,'FontSize',8);
            ylim([0.009 10])
            xlabel('Time, t (samples)')
            ylabel('NMSD')
            F = GFigure.captureCurrentFigure();
            
        end
        
        % Figures 4 and 5
        function F = experiment_4(obj, niter)
            % Comparing NMSD and EIER of TISO and TIRSO with other
            % competing algorithms (OSGDTISO and PGDTIRSO).
            % using topologySimulation.m, toplogySimulationXvalNMSDEIER.m,
            % and plotting function plotFigures4and5.m
            % Stationary setting
            addpath Personal/bakht/TISO_TIRSO % TODO: move the functions to
            % their right place (ProcessingBlocks)
            
            %% Configuration
            noOfMCitr   =  200;
            noOfNodes   = 10;
            filterOrder = 2;
            niterPGD    = 5;
            noOfTimeInstants = 3000;
            forgettingFactor = 0.99;
            noiseSigma  = 0.01;
            regParRange=[0.1,-7];
            noOfRegpars= 50;
            erdosRenyiEdgeProbability = 0.2;
            threshold   = 0;
            stepsizeFactor = 0.1;
            ch_algorithms = {'TISO', 'OSGDTISO', 'TIRSO', 'PGDTIRSO'};
            
            %% Computation            
            n_algs = length(ch_algorithms); % initialization
            c_NMSD          = cell(n_algs, 1);
            c_EIER          = cell(n_algs, 1);
            c_m_VARcoefs    = cell(n_algs, 1);
            c_t_A_estimates = cell(n_algs, 1);
            
            generator = VARSTFunctionGenerator;
            generator.sigma = noiseSigma;
            
            for i_algo=1:n_algs
                ch_alg = ch_algorithms{i_algo};
                [regPar_xvalNMSD(i_algo),regpar_xvalEIER(i_algo)]= toplogySimulationXvalNMSDEIER(...
                    ch_alg, noOfNodes, filterOrder, noOfTimeInstants, ...
                    regParRange, noOfRegpars,...
                    generator, erdosRenyiEdgeProbability, threshold,...
                    forgettingFactor, niterPGD, stepsizeFactor); %#ok<AGROW>
            end
            for idx_alg = 1:n_algs
                ch_alg = ch_algorithms{idx_alg};
                [c_NMSD{idx_alg}, ~, c_m_VARcoefs{idx_alg},...
                    c_t_A_estimates{idx_alg}] = ...
                    topologySimulation( ch_alg, ...
                    noOfNodes, filterOrder, noOfMCitr, noOfTimeInstants, ...
                    generator, erdosRenyiEdgeProbability, regPar_xvalNMSD(idx_alg), threshold,...
                    forgettingFactor, niterPGD, stepsizeFactor);
                [~, c_EIER{idx_alg}, c_m_VARcoefs{idx_alg}, ~] = ...
                    topologySimulation( ch_alg, ...
                    noOfNodes, filterOrder, noOfMCitr, noOfTimeInstants, ...
                    generator, erdosRenyiEdgeProbability, regpar_xvalEIER(idx_alg), threshold,...
                    forgettingFactor, niterPGD, stepsizeFactor);
            end
            
            %% Plotting
            markerDistance = 15;
            ch_styles = {'- ^ r', ': s b', '- o g', ': x k'};
            addpath Plotters
            F = plotFigures4and5_TIRSO(ch_algorithms, c_m_VARcoefs{1}(:,:,end), ...
                noOfTimeInstants, filterOrder, markerDistance, c_EIER, ...
                c_NMSD, c_t_A_estimates, ch_styles);
        end
        
        % Figure 6
        function F = experiment_6(obj, niter)
            % Plotting the averaged estimated graph via TISO and TIRSO
            % for one realization at different time instants.
            % Adaptive, non-diminishing stepsize
            % Stationary VAR Process
            % Fixed threshold for identifying nonzero edge
            %% Parameters %%
            noOfNodes        = 12;
            filtOrder        = 2;
            noOfExperiments  = 1;
            noOfObservations = 1000;
            erdosRenyiEdgeProbability = 0.2;
            syntheticSigma   = 0.005;
            initialSigma     = 0.01;
            f_factor         = 0.98;
            
            %% Synthetic Data Generator
            rng(3);
            m_adjacency = (rand(noOfNodes)<erdosRenyiEdgeProbability)...
                .*not(eye(noOfNodes));
            myGraph = Graph('m_adjacency',m_adjacency);
            myCoefficients1 = VARSTFunctionGenerator.randomCoefficientsFromGraph(myGraph, filtOrder);
            gen = VARSTFunctionGenerator;
            gen.nTimeSamples = noOfObservations;
            gen.sigma = syntheticSigma;
            myCoefficients = myCoefficients1;
            gen.t_coefficients = myCoefficients;
            trueWeightedAdjacency=VARSTFunctionGenerator.VARCoefficients2AdjacencyMatrix(myCoefficients);
            
            %% Estimator construction and call to simulator
            v_regPars  = [1e-6]; %#ok<NBRAK>
            n_regPars  = length(v_regPars);
            t_estimatedCoefficients= zeros(noOfNodes,noOfNodes,noOfExperiments,noOfObservations,2*n_regPars);
            sq_error      = zeros (noOfExperiments, noOfObservations, 2*n_regPars);
            predict_error = zeros (noOfExperiments, noOfObservations, 2*n_regPars+1);
            power         = zeros (noOfExperiments, noOfObservations, 2*n_regPars+1);
            r_EIE         = zeros (noOfExperiments, noOfObservations, 2*n_regPars+1);
            r_D           = zeros (noOfExperiments, noOfObservations, 2*n_regPars+1);
            r_FA          = zeros (noOfExperiments, noOfObservations, 2*n_regPars+1);
            tirsoes       = cell  (1, n_regPars);
            tisoes        = cell  (1, n_regPars);
            
            ltc = LoopTimeControl(noOfExperiments);
            for exp=1:noOfExperiments
                for i_r = 1:n_regPars
                    tirsoes{i_r} = RecomidOGLE;
                    tirsoes{i_r}.b_selectComid = 0;
                    tirsoes{i_r}.b_adaptiveStepsize = 1;
                    tirsoes{i_r}.sigma = initialSigma;
                    tirsoes{i_r}.b_diminishingStepsize = 0;
                    tirsoes{i_r}.ss_factor = 8;
                    tirsoes{i_r}.forgettingFactor = f_factor;
                    tirsoes{i_r}.mulFactor=1-f_factor;
                    %
                    tirsoes{i_r}.regPar = v_regPars(i_r);
                    %Calling Simulator
                    [sq_error(exp,:,i_r), predict_error(exp,:,i_r), ...
                        power(exp,:,i_r), t_estimatedCoefficients(:,:,exp,:,i_r), ...
                        r_EIE(exp,:,i_r), r_D(exp,:,i_r), r_FA(exp,:,i_r)]...
                        =  OMVARSimulator(gen, tirsoes{i_r}...
                        , filtOrder, noOfObservations, myCoefficients);
                    tisoes{i_r} = RecomidOGLE;
                    tisoes{i_r}.b_selectComid = 1;
                    tisoes{i_r}.forgettingFactor = 0;
                    tisoes{i_r}.b_maxInComid = 0;
                    tisoes{i_r}.b_adaptiveStepsize = 1;
                    tisoes{i_r}.sigma = initialSigma;
                    tisoes{i_r}.b_diminishingStepsize = 0;
                    tisoes{i_r}.ss_factor = 8;
                    %
                    tisoes{i_r}.regPar = v_regPars(i_r);
                    %Calling Simulator
                    [sq_error(exp,:,n_regPars+i_r), predict_error(...
                        exp,:,n_regPars+i_r),power(exp,:,n_regPars+i_r), ...
                        t_estimatedCoefficients(:,:,exp,:,n_regPars+i_r), ...
                        r_EIE(exp,:,n_regPars+i_r), r_D(exp,:,n_regPars+i_r), ...
                        r_FA(exp,:,n_regPars+i_r)] ...
                        = OMVARSimulator(gen, tisoes{i_r}, filtOrder, ...
                        noOfObservations, myCoefficients);
                end
                
                % Genie-aided
                [~, predict_error(exp,:,2*n_regPars+1), power(exp,:,2*n_regPars+1)] =...
                    OMVARSimulator(gen, 'genie', filtOrder, noOfObservations, myCoefficients);
                
                ltc.go(exp)
            end
            
            %% Calculating NMSD and plotting
            m_NMSD = zeros(2*n_regPars, size(sq_error, 2));
            m_EIER = zeros(2*n_regPars, size(sq_error, 2));
            normOfCoeff=zeros(1,noOfObservations);
            if ndims(myCoefficients)==3
                normOfCoeff(:) = norm(vec(myCoefficients)).^2;
            else
                for t=1:noOfObservations
                    normOfCoeff(t)= norm(vec(myCoefficients (:,:,:,t)))^2;
                end
            end
            for i = 1:2*n_regPars
                devNumerator=mean(sq_error(:,:,i),1);
                m_NMSD(i,:)=devNumerator./normOfCoeff;
                m_EIER(i,:) = mean(r_EIE(:,:,i),1);
            end
            m_NMSE = zeros(2*n_regPars+1, size(predict_error,2));
            m_filteredNMSE = zeros(2*n_regPars+1, size(predict_error, 2));
            avgLength = 100;
            for k = 1:2*n_regPars+1
                m_NMSE(k,:)=mean(predict_error(:,:,k),1) ./ mean( power(:,:,k),1);
                m_filteredNMSE(k,:) = filter(ones(1, avgLength)/avgLength, 1, m_NMSE(k,:));
            end
            G=digraph(trueWeightedAdjacency','omitselfloops'); % transpose the adjacency matrix
            Weights=G.Edges.Weight;
            positionVector = [0.0312, 0.264, 0.240, 0.49];
            subplot('Position',positionVector)
            plot(G,'layout','circle','EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
            set(gca,'XTick',[],'YTick',[])
            colorbar;
            title('True Graph')
            time=int16([ceil(noOfObservations/3) ceil(2*noOfObservations)/3 noOfObservations]);
            noOfEdges=ceil(erdosRenyiEdgeProbability*noOfNodes*(noOfNodes-1));
            %% Graphs estimated via TIRSO
            % 1
            A_tirso1= mean(t_estimatedCoefficients(:,:,1,1:time(1),1),4);
            %Computing threshold
            sortedA=sort(A_tirso1(:),'descend');
            edgeThresholdtirso=sortedA(noOfEdges);
            A_tirso1thresholded=(A_tirso1 > edgeThresholdtirso).*A_tirso1;
            Gtirso1=digraph(A_tirso1thresholded','omitselfloops');
            Weights=Gtirso1.Edges.Weight;
            subplot(2,4,2);
            plot(Gtirso1,'layout','circle','EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
            title('Average Graph via TIRSO from 0 to T/3')
            set(gca,'XTick',[],'YTick',[])
            % 2
            A_tirso2= mean(t_estimatedCoefficients(:,:,1,time(1):time(2),1),4);
            A_tirso2thresholded=(A_tirso2 > edgeThresholdtirso).*A_tirso2;
            Gtirso2=digraph(A_tirso2thresholded','omitselfloops');
            Weights=Gtirso2.Edges.Weight;
            subplot(2,4,3);
            plot(Gtirso2,'layout','circle','EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
            title('Average Graph via TIRSO from T/3 to 2T/3')
            set(gca,'XTick',[],'YTick',[])
            % 3
            A_tirso3= mean(t_estimatedCoefficients(:,:,1,time(2):time(3),1),4);
            A_tirso3thresholded=(A_tirso3 > edgeThresholdtirso).*A_tirso3;
            Gtirso3=digraph(A_tirso3thresholded','omitselfloops');
            Weights=Gtirso3.Edges.Weight;
            subplot(2,4,4);
            plot(Gtirso3,'layout','circle','EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
            title('Average Graph via TIRSO from 2T/3 to T')
            set(gca,'XTick',[],'YTick',[])
            %% graphs estimated via TISO
            % 1
            A_tiso1= mean(t_estimatedCoefficients(:,:,1,1:time(1),2),4);%Computing threshold
            sortedA=sort(vec(A_tirso1),'descend');
            edgeThresholdtiso=sortedA(noOfEdges);
            A_tiso1thresholded=(A_tiso1 > edgeThresholdtiso).*A_tiso1;
            Gtiso1=digraph(A_tiso1thresholded','omitselfloops');
            Weights=Gtiso1.Edges.Weight;
            subplot(2,4,6);
            plot(Gtiso1,'layout','circle','EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
            title('Average Graph via TISO from 0 to T/3')
            set(gca,'XTick',[],'YTick',[])
            % 2
            A_tiso2= mean(t_estimatedCoefficients(:,:,1,time(1):time(2),2),4);
            A_tiso2thresholded=(A_tiso2 > edgeThresholdtiso).*A_tiso2;
            Gtiso2=digraph(A_tiso2thresholded','omitselfloops');
            Weights=Gtiso2.Edges.Weight;
            subplot(2,4,7);
            plot(Gtiso2,'layout','circle','EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
            title('Average Graph via TISO from T/3 to 2T/3')
            set(gca,'XTick',[],'YTick',[])
            % 3
            A_tiso3= mean(t_estimatedCoefficients(:,:,1,time(2):time(3),2),4);
            A_tiso3thresholded=(A_tiso3 > edgeThresholdtiso).*A_tiso3;
            Gtiso3=digraph(A_tiso3thresholded','omitselfloops');
            Weights=Gtiso3.Edges.Weight;
            subplot(2,4,8);
            plot(Gtiso3,'layout','circle','EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
            title('Average Graph via TISO from 2T/3 to T')
            set(gca,'XTick',[],'YTick',[])
            F(1)=GFigure.captureCurrentFigure();
        end
        
        % Figure 7
        function F = experiment_7(obj, niter)
            % Plotting NMSD and NMSE for different values of forgetting
            % factor \gamma
            % Adaptive, non-diminishing stepsize
            % Smooth-Transition nonstationary VAR Process
            % NMSD and NMSE
            noOfNodes        = 12;
            filtOrder        = 2;
            noOfObservations = 3000; 
            transitionTime   = noOfObservations/3;
            transitionPar    = TimeVaryingVARSTFunctionGenerator.stParameter(...
                0.99, transitionTime);
            % determines the speed of the transition(the higher, the faster)
            noOfExperiments  = 50; % 300 Monte Carlo iterations
            erdosRenyiEdgeProbability = 0.2;
            syntheticSigma   = 0.005;
            regPar           = 1e-6;
            initialSigma     = 0.01;
            
            %% Estimator construction and call to simulator
            v_fFactors = [0.85 0.95 0.98 0.99];
            n_fFactors = length(v_fFactors);
            sq_error = zeros(noOfExperiments, noOfObservations, 4);
            predict_error= zeros(noOfExperiments, noOfObservations, 4);
            power= zeros (noOfExperiments, noOfObservations, 4);
            tirsoes =cell(1, n_fFactors);
            ltc = LoopTimeControl(noOfExperiments);
            for exp=1:noOfExperiments
                %% Synthetic Data Generator
                m_adjacency = (rand(noOfNodes)<erdosRenyiEdgeProbability)...
                    .*not(eye(noOfNodes)); %directed graph, noSelfLoops
                myGraph = Graph('m_adjacency',m_adjacency);
                myCoefficients1 = VARSTFunctionGenerator.randomCoefficientsFromGraph(myGraph, filtOrder);
                myCoefficients2 = VARSTFunctionGenerator.randomCoefficientsFromGraph(myGraph, filtOrder);
                gen = TimeVaryingVARSTFunctionGenerator;
                gen.N = noOfNodes;
                gen.Q = filtOrder;
                gen.T = noOfObservations;
                gen.sigma = syntheticSigma;
                [myCoefficients, ~, ~] = gen.smoothTransitionModel(myCoefficients1,myCoefficients2, ...
                    transitionPar, transitionTime, noOfObservations);
                gen.t_coefficients = myCoefficients;
                for i_f = 1:n_fFactors
                    tirsoes{i_f} = RecomidOGLE;
                    tirsoes{i_f}.regPar = regPar;
                    tirsoes{i_f}.b_selectComid = 0;
                    tirsoes{i_f}.b_adaptiveStepsize = 1;
                    tirsoes{i_f}.sigma = initialSigma;
                    tirsoes{i_f}.b_diminishingStepsize = 0;
                    tirsoes{i_f}.ss_factor = 8;
                    %
                    tirsoes{i_f}.forgettingFactor = v_fFactors(i_f);
                    tirsoes{i_f}.mulFactor=1-v_fFactors(i_f);
                    %Calling Simulator
                    [sq_error(exp,:,i_f), predict_error(exp,:,i_f), ...
                        power(exp,:,i_f)] =  OMVARSimulator(gen, tirsoes{i_f}...
                        , filtOrder, noOfObservations, myCoefficients,0.15,1);
                end
                tiso = RecomidOGLE;
                tiso.regPar = regPar;
                tiso.b_selectComid = 1;
                tiso.forgettingFactor = 0;
                tiso.b_maxInComid = 0;
                tiso.b_adaptiveStepsize = 1;
                tiso.sigma = initialSigma;
                tiso.b_diminishingStepsize = 0;
                tiso.ss_factor = 8;
                %Calling Simulator
                [sq_error(exp,:,n_fFactors+1), predict_error(...
                    exp,:,n_fFactors+1),power(exp,:,n_fFactors+1)] ...
                    = OMVARSimulator(gen, tiso, filtOrder, ...
                    noOfObservations, myCoefficients,0.15,1);
                
                % Genie-aided
                [~, predict_error(exp,:,n_fFactors+2), power(exp,:,n_fFactors+2)] =...
                    OMVARSimulator(gen, 'genie', filtOrder, noOfObservations, myCoefficients,0.15,1);
                
                ltc.go(exp)
            end
            %% Calculating NMSD and plotting
            m_NMSD = zeros(n_fFactors+1, size(sq_error, 2));
            normOfCoeff=zeros(1,noOfObservations);
            for t=1:noOfObservations
                normOfCoeff(t)= norm(vec(myCoefficients (:,:,:,t)))^2;
            end
            for i = 1:n_fFactors+1
                devNumerator=mean(sq_error(:,:,i),1);
                m_NMSD(i,:)=devNumerator./normOfCoeff;
            end
            m_NMSE = zeros(n_fFactors+2, size(predict_error,2));
            m_filteredNMSE = zeros(n_fFactors+2, size(predict_error, 2));
            avgLength = 100;
            for k = 1:n_fFactors+2
                m_NMSE(k,:)= mean(predict_error(:,:,k),1) ./ mean( power(:,:,k),1);
                m_filteredNMSE(k,:) = filter(ones(1, avgLength)/avgLength, 1, m_NMSE(k,:));
            end
            
            M(1,1)=GFigure();
            M(1,1).b_logy=1;
            for i = 1:n_fFactors
                M(1,1).c_legend{i} = sprintf('TIRSO, \\gamma=%0.5g',...
                    tirsoes{i}.forgettingFactor);
            end
            M(1,1).c_legend{n_fFactors+1} = 'TISO';
            M(1,1).c_styles{n_fFactors+1}= '--';
            M(1,1).m_Y= m_NMSD;
            M(1,1).ch_ylabel='NMSD';
            M(1,1).ch_xlabel='Time';
            M(2,1)=GFigure();
            M(2,1).b_logy=1;
            M(2,1).c_legend=M(1,1).c_legend;
            M(2,1).c_legend{n_fFactors+2} = 'Genie';
            M(2,1).c_styles{n_fFactors+1}= '--';
            M(2,1).m_Y= m_filteredNMSE;
            M(2,1).ch_ylabel='NMSE';
            M(2,1).ch_xlabel='Time';
            F=GFigure('m_multiplot',M);
        end
        
        % Figure 8
        function F = experiment_8(obj, niter)
            % Assessment of impact of regPar
            % Plotting NMSD and NMSE 
            % Adaptive, non-diminishing stepsize
            % Smooth-Transition nonstationary VAR Process
            noOfNodes        = 12;
            filtOrder        = 2;
            noOfExperiments  = 200; % 200
            noOfObservations = 3000;
            transitionTime   = noOfObservations/3;
            transitionPar    = TimeVaryingVARSTFunctionGenerator.stParameter(...
                0.99, transitionTime);
            erdosRenyiEdgeProbability = 0.2;
            syntheticSigma            = 0.005;
            
            initialSigma = 0.01;
            f_factor     = 0.98;
            
            %% Estimator construction and call to simulator
            v_regPars  = [1e-5 1e-6 1e-7];
            n_regPars  = length(v_regPars);
            
            sq_error      = zeros (noOfExperiments, noOfObservations, 2*n_regPars);
            predict_error = zeros (noOfExperiments, noOfObservations, 2*n_regPars+1);
            power         = zeros (noOfExperiments, noOfObservations, 2*n_regPars+1);
            tirsos        = cell  (1, n_regPars);
            tisos         = cell  (1, n_regPars);
            
            ltc = LoopTimeControl(noOfExperiments);
            for exp=1:noOfExperiments
                %% Synthetic Data Generator
                m_adjacency = (rand(noOfNodes)<erdosRenyiEdgeProbability)...
                    .*not(eye(noOfNodes)); %directed graph, noSelfLoops
                myGraph = Graph('m_adjacency',m_adjacency);
                
                myCoefficients1 = VARSTFunctionGenerator.randomCoefficientsFromGraph(myGraph, filtOrder);
                myCoefficients2 = VARSTFunctionGenerator.randomCoefficientsFromGraph(myGraph, filtOrder);
                gen   = TimeVaryingVARSTFunctionGenerator;
                gen.N = noOfNodes;
                gen.Q = filtOrder;
                gen.T = noOfObservations;
                gen.sigma = syntheticSigma;
                [myCoefficients, t_A, transitionCurve] = gen.smoothTransitionModel(myCoefficients1,myCoefficients2, ...
                    transitionPar, transitionTime, noOfObservations); %#ok<ASGLU>
                gen.t_coefficients = myCoefficients;
                for i_r = 1:n_regPars
                    % TIRSO
                    tirsos{i_r} = RecomidOGLE;
                    tirsos{i_r}.b_selectComid = 0;
                    tirsos{i_r}.b_adaptiveStepsize = 1;
                    tirsos{i_r}.sigma = initialSigma;
                    tirsos{i_r}.b_diminishingStepsize = 0;
                    tirsos{i_r}.ss_factor = 8;
                    tirsos{i_r}.forgettingFactor = f_factor;
                    tirsos{i_r}.mulFactor=1-f_factor;
                    %
                    tirsos{i_r}.regPar = v_regPars(i_r);
                    %Calling Simulator
                    [sq_error(exp,:,i_r), predict_error(exp,:,i_r), ...
                        power(exp,:,i_r)] =  OMVARSimulator(gen, tirsos{i_r}...
                        , filtOrder, noOfObservations, myCoefficients,0.15,1);
                    
                    % TISO, instantaneous stepsize
                    tisos{i_r} = RecomidOGLE;
                    tisos{i_r}.b_selectComid = 1;
                    tisos{i_r}.forgettingFactor = 0;
                    tisos{i_r}.b_maxInComid = 0;
                    tisos{i_r}.b_adaptiveStepsize = 1;
                    tisos{i_r}.sigma = initialSigma;
                    tisos{i_r}.b_diminishingStepsize = 0;
                    tisos{i_r}.ss_factor = 8;
                    %
                    tisos{i_r}.regPar = v_regPars(i_r);
                    %Calling Simulator
                    [sq_error(exp,:,n_regPars+i_r), predict_error(...
                        exp,:,n_regPars+i_r),power(exp,:,n_regPars+i_r)] ...
                        = OMVARSimulator(gen, tisos{i_r}, filtOrder, ...
                        noOfObservations, myCoefficients,0.15,1);
                end
                
                % Genie-aided
                [~, predict_error(exp,:,2*n_regPars+1), power(exp,:,2*n_regPars+1)] =...
                    OMVARSimulator(gen, 'genie', filtOrder, noOfObservations, myCoefficients,0.15,1);
                
                ltc.go(exp)
            end
            
            %% Calculating NMSD and plotting
            m_NMSD = zeros(2*n_regPars, size(sq_error, 2));
            normOfCoeff=zeros(1,noOfObservations);
            for t=1:noOfObservations
                normOfCoeff(t)= norm(vec(myCoefficients (:,:,:,t)))^2;
            end
            for i = 1:2*n_regPars
                devNumerator=mean(sq_error(:,:,i),1);
                m_NMSD(i,:)=devNumerator./normOfCoeff;
            end
            m_NMSE = zeros(2*n_regPars+1, size(predict_error,2));
            m_filteredNMSE = zeros(2*n_regPars+1, size(predict_error, 2));
            avgLength = 100;
            for k = 1:2*n_regPars+1
                m_NMSE(k,:)=mean(predict_error(:,:,k),1) ./ mean( power(:,:,k),1);
                m_filteredNMSE(k,:) = filter(ones(1, avgLength)/avgLength, 1, m_NMSE(k,:));
            end
            
            M(1,1)=GFigure();
            M(1,1).b_logy=1;
            for i = 1:n_regPars
                M(1,1).c_legend{i} = sprintf(...
                    'TIRSO, \\lambda=%0.5g', tirsos{i}.regPar);
                M(1,1).c_styles{i} = '-';
                M(1,1).c_legend{n_regPars+i} = sprintf(...
                    'TISO, \\lambda=%0.5g',  tisos{i}.regPar);
                M(1,1).c_styles{n_regPars+i} = '--';
            end
            M(1,1).colorPeriod = n_regPars;
            M(1,1).m_Y= m_NMSD;
            M(1,1).ch_ylabel='NMSD';
            M(1,1).ch_xlabel='Time';
            
            M(2,1)=GFigure();
            M(2,1).b_logy=1;
            M(2,1).c_legend=M(1,1).c_legend;
            M(2,1).c_legend{2*n_regPars+1} = 'Genie';
            M(2,1).colorPeriod = n_regPars;
            M(2,1).c_styles = M(1,1).c_styles;
            M(2,1).c_styles{2*n_regPars+1} = '-.';
            M(2,1).m_Y= m_filteredNMSE;
            M(2,1).ch_ylabel='NMSE';
            M(2,1).ch_xlabel='Time';
            F=GFigure('m_multiplot',M);
            
        end
        
        %% Experiments with real data from Lundin Norway
        % Figure 9
        function F = experiment_9(obj, niter)
            % shows the prediction error vs. pred. horizon 
            % for all 24 variables of the Lundin_System20 dataset obtained from 
            % Lundin Norway.
            nSensors      = 24;
            nHours        = 24;
            nTimeInstants = 6*60*nHours;            
            n_toPredict   = 9;
            
            regPar        = 1e-5;
            fFactor       = 0.9; 
            sFactor       = 2;
            fOrder        = 8;

            load('m_X_total.mat','m_X_totalData', 'samplingInterval', 'c_tagNames'); 
            
            m_X = m_X_totalData(1:nSensors, 1:nTimeInstants);
            %nSensors = size(m_X,1);
                        
            %initialization
            t_pE_indiv     = nan([nSensors, nTimeInstants, n_toPredict]);
            v_pE           = zeros([nTimeInstants, n_toPredict]);
            v_sampPow      = zeros(nTimeInstants, 1);
            m_NMSE         = zeros(size(v_pE));
            m_filteredNMSE = zeros(size(v_pE));
            avgLength      = 20;
            
            %% TIRSO configuration
            tirso = RecomidOGLE;
            tirso.b_adaptiveStepsize = 1;
            tirso.b_diminishingStepsize = 0;
            tirso.sigma = 0.01;
            
            tirso.ss_factor = sFactor;
            tirso.forgettingFactor = fFactor;
            tirso.mulFactor = 1-fFactor;
            tirso.regPar = regPar;
            if fFactor==0
                tirso.b_selectComid = 1;
%                 ch_legendEntry = 'TISO';
            else
                tirso.b_selectComid = 0;
%                 ch_legendEntry = sprintf(...
%                     'TIRSO, \\\\gamma=%0.5g', fFactor);
            end
            
            %% Estimation
            [v_pE(:,:), v_sampPow(:), ...
                ~, ~, t_pE_indiv(:,:,:), ...
                m_sampPow_indiv(:,:)] = ...
                obj.analyzeRealData(m_X, tirso, fOrder, n_toPredict);
            m_NMSE(:,:)=v_pE(:,:)./...
                mean( v_sampPow(:));
            t_NMSE_indiv(:,:,:) = ...
                t_pE_indiv(:,:,:)./...
                mean(m_sampPow_indiv(:,:), 2);
            m_filteredNMSE(:,:) = filter(ones(1, avgLength)/...
                avgLength, 1, m_NMSE(:,:));

            %% Plotting
            c_legend = c_tagNames;
            m_NMSE_indiv_vs_horizon = squeeze(mean(t_NMSE_indiv, 2, 'omitnan'));
            
            M(1,1)=GFigure('b_logy', 1, 'ch_title', 'NMSE');
            M(1,1).ch_xlabel = 'Prediction horizon (s)';
            M(1,1).c_legend = c_legend;
            M(1,1).m_Y = m_NMSE_indiv_vs_horizon;
            M(1,1).m_X = samplingInterval*(1:n_toPredict);
            
            F = M;
        end
        
        % Figure 10
        function F = experiment_10(obj, niter)
            % Plotting the recovered graph for the real data from the
            % Lundin_System20 dataset obtained from Lundin Norway.
            % 12 hours of operation
            
            nSensors      = 24; 
            nHours        = 12; 
            nTimeInstants = 6*60*nHours;
            n_toPredict   = 1; %!
            
            regPar        = 2e-2;
            fFactor       = 0.9; 
            sFactor       = 2; %scale factor of step size
            fOrder        = 8;
                        
            load('m_X_total.mat','m_X_totalData', 'c_tagNames'); 
            m_X = m_X_totalData(1:nSensors, 1:nTimeInstants);

            %% TIRSO configuration
            tirso = RecomidOGLE;
            tirso.b_adaptiveStepsize = 1;
            tirso.b_diminishingStepsize = 0;
            tirso.sigma = 0.01;
            
            tirso.ss_factor = sFactor;
            tirso.forgettingFactor = fFactor;
            tirso.mulFactor = 1-fFactor;
            tirso.regPar = regPar;
            if fFactor==0
                tirso.b_selectComid = 1;
            else
                tirso.b_selectComid = 0;
            end
            
            %%  Estimation
            [~, ~, t_Atirso] = obj.analyzeRealData(m_X, tirso, fOrder, n_toPredict);
                      
            %% Forming the graph
            t_norms = squeeze(norms(t_Atirso, 2, 3));
            v_avgNorms= mean(t_norms, 3, 'omitNaN');            
            v_avgNorms(1:(nSensors+1):end)=0;  % avoiding selfloops to plot
            [~, order] = sort(v_avgNorms(:), 'descend');
            n_toPreserve = 2.5*nSensors; % displaying the edges in the order of the number of nodes
            threshold = v_avgNorms(order(n_toPreserve));
            AdjacencyMatrix = ((v_avgNorms>threshold).*v_avgNorms)';% transposed
            nodes_labels = c_tagNames;
            G = digraph(AdjacencyMatrix, nodes_labels); 
            Weights=G.Edges.Weight;
            
            %% Plotting
            plot(G,'layout','circle','EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
            colorbar
            %title('Average estimated graph')
            set(gca,'XTick',[],'YTick',[])
            
            F(1)=GFigure.captureCurrentFigure();
        end
    end
    
    methods (Static)
        %
        function [m_pred_error, v_samplePower, t4_estimatedCoefficients,...
                t_A, t_pred_error_indiv, m_samplePower_indiv] = ...
                analyzeRealData(m_X, eTemplate, in_fOrder, n_toPredict)
            % This routine takes a matrix containing a stream of real data
            % and analyzes it using the method defined by the object
            % provided as an input in eTemplate, for example a Tirso/Tiso,
            % etc.
            if not(exist('n_toPredict', 'var'))
                n_toPredict = 1;
                warning('n_toPredict not provided. Assigned to 1');
            end
            [nSensors, nTimeInstants] = size(m_X);            
            
            assert(isa(eTemplate, 'OnlineGroupLassoEstimator'));
            
            ove = OnlineModularVAREstimator;
            ove.nTimeSeries = nSensors;
            ove.filterOrder = in_fOrder;
            ove.b_leastSquareSelfLoops = 0;
            createOgles(ove, eTemplate);
            inheritProperties(ove, eTemplate);
            initializeOgles(ove);
            
            % initialize arrays with metric outputs
            m_pred_error = nan(nTimeInstants, n_toPredict);
            t_pred_error_indiv = nan(nSensors, nTimeInstants, n_toPredict);
            v_samplePower= zeros(nTimeInstants, 1);
            m_samplePower_indiv = zeros(nSensors, nTimeInstants);
            t4_estimatedCoefficients = zeros(nSensors, ...
                nSensors, in_fOrder, nTimeInstants);
            t_A = zeros(nSensors, nSensors, nTimeInstants);
            
            % Main loop
            ltc = LoopTimeControl(nTimeInstants);
            for t = 1:nTimeInstants
                v_samplePower(t) = norm(m_X(:,t))^2;
                m_samplePower_indiv(:, t) = m_X(:,t).^2;
                aux_n_toPredict = min(n_toPredict, nTimeInstants-t+1);
                setTimeIndex(ove, t);
                if not(isempty(ove.t_estimatedCoefficients))
                    m_pX= ove.predict_static(ove.t_estimatedCoefficients...
                        , m_X(:, (t-in_fOrder):(t-1)),n_toPredict);
                    m_pred_errors_now = (m_pX(:,1:aux_n_toPredict) ...
                        - m_X(:,t-1+(1:aux_n_toPredict))).^2;
                    v_pred_errors = sum(m_pred_errors_now);
                    
                    t4_estimatedCoefficients(:,:,:,t) = ...
                        ove.t_estimatedCoefficients;
                else
                    v_pred_errors = nan(n_toPredict, 1);
                    m_pred_errors_now = nan(nSensors, n_toPredict);
                    %this line is probably not necessary
                end
                receiveSingleSample(ove, m_X(:,t));
                t_A(:,:,t) = ove.estimatedAdjacencyMatrix();
                
                for tau = 1:aux_n_toPredict
                    m_pred_error(t-1+tau, tau) = v_pred_errors(tau);
                    t_pred_error_indiv(:, t-1+tau, tau) = m_pred_errors_now(:, tau);
                end
                
                ltc.go(t);
            end
        end
        function [P_MD,P_FA,EIER]= computeEierMissdetectionFalsealarm(threshold,noOfNodes,trueWeightedAdjacency,estAdj)
            % used in experiment 2 to compute graph detection metrics
            noOfTrueEdges = nnz(trueWeightedAdjacency);
            noOfPossibleEdges= noOfNodes.^2 - noOfNodes;
            P_MD = nnz(trueWeightedAdjacency.*not(estAdj> threshold))/noOfTrueEdges;
            P_FA  = nnz(not(trueWeightedAdjacency).*(estAdj> threshold))/(noOfNodes.^2 - noOfTrueEdges);
            EIER = (nnz(trueWeightedAdjacency - (estAdj>threshold)))/noOfPossibleEdges;
        end
        
    end
end