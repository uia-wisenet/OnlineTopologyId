function [sq_errorTirso,sq_errorTiso,predict_errorTirso, predict_errorTiso,power,t_Atirso,t_Atiso]= ...
                    computationsForTirsoTiso(noOfNodes, n_regPars,v_regPars, f_factor, gen, grelos, glos, filtOrder, ...
                    noOfObservations, myCoefficients,exp,initialSigma)
                sq_errorTirso= zeros(noOfObservations,n_regPars);
                predict_errorTirso = zeros(noOfObservations,n_regPars);
                sq_errorTiso= zeros(noOfObservations,n_regPars);
                predict_errorTiso = zeros(noOfObservations,n_regPars);
                t_Atirso= zeros(noOfNodes,noOfNodes, noOfObservations,n_regPars);
                t_Atiso= zeros(noOfNodes,noOfNodes, noOfObservations,n_regPars);

                power = zeros(1,noOfObservations);
                for i_r = 1:n_regPars
                    %TIRSO COMID
                    grelos{i_r} = RecomidOGLE;
                    grelos{i_r}.b_selectComid = 0;
                    grelos{i_r}.b_adaptiveStepsize = 1;
                    grelos{i_r}.sigma = initialSigma;
                    grelos{i_r}.b_diminishingStepsize = 0;
                    grelos{i_r}.ss_factor = 8;
                    grelos{i_r}.forgettingFactor = f_factor;
                    grelos{i_r}.mulFactor=1-f_factor;
                    %
                    grelos{i_r}.regPar = v_regPars(i_r);
                    %Calling Simulator
                    [sq_errorTirso(:,i_r), predict_errorTirso(:,i_r), ...
                        power, t_Atirso(:,:,:,i_r), ~, ~...
                        , ~] =  OMVARSimulator(gen, grelos{i_r}...
                        , filtOrder, noOfObservations, myCoefficients);

                    % TISO COMID, instantaneous stepsize
                    glos{i_r} = RecomidOGLE;
                    glos{i_r}.b_selectComid = 1;
                    glos{i_r}.forgettingFactor = 0;
                    glos{i_r}.b_maxInComid = 0;
                    glos{i_r}.b_adaptiveStepsize = 1;
                    glos{i_r}.sigma = initialSigma;
                    glos{i_r}.b_diminishingStepsize = 0;
                    glos{i_r}.ss_factor = 8;
                    %
                    glos{i_r}.regPar = v_regPars(i_r);
                    %Calling Simulator
                    [sq_errorTiso(:,i_r), predict_errorTiso(:,i_r),~, ...
                        t_Atiso(:,:,:,i_r), ~,~,~] ...
                        = OMVARSimulator(gen, glos{i_r}, filtOrder, ...
                        noOfObservations, myCoefficients);
                   
                end
end