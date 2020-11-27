function F = plotFigures4and5_TIRSO (ch_algorithms, m_VARcoefs, T, P, ...
    markerDistance, c_EIER, c_NMSD, c_t_A_estimates, ch_styles)
% Plotting for Figs. 4 and 5
filterOrder = P;
noOfTimeInstants = T;
n_algs = length(ch_algorithms);
for k = 1:3
    v_figures(k) = figure(100-k); clf
end
figure(v_figures(3));
subplot(1, n_algs+1, 1);
imagesc(m_VARcoefs)
caxis_ref = caxis;
title('True VAR coefficients')
            
v_idx_Xaxis = 20:noOfTimeInstants;
v_timeAxis=filterOrder+1:noOfTimeInstants;
for idx_alg = 1:n_algs
    markerOffset = idx_alg*markerDistance;
    markerIndices = filterOrder + markerOffset ...
        + 0:(markerDistance*n_algs):noOfTimeInstants;
    figure(v_figures(1));
    semilogy(v_idx_Xaxis, c_EIER{idx_alg}(v_idx_Xaxis), ...
        ch_styles{idx_alg},'MarkerIndices',...
        markerIndices(markerIndices<length(v_idx_Xaxis)),...
        'LineWidth', 1.5);
    hold on
    figure(v_figures(2));
    semilogy(v_timeAxis, c_NMSD{idx_alg}(v_timeAxis),...
        ch_styles{idx_alg}, 'MarkerIndices',...
        markerIndices(markerIndices<length(v_timeAxis)),...
        'LineWidth', 1.5)
    hold on
    %title(['NMSD vs Time instants for ', num2str(noOfMCitr),' experiments'])
    figure(v_figures(3)); % plotting the estimates at T via imagesc
    subplot(1,n_algs+1, idx_alg+1)
    imagesc(c_t_A_estimates{idx_alg}(:,:,end,end))
    caxis(caxis_ref);
    title(ch_algorithms(idx_alg));
end
ch_ylabels = {'EIER', 'NMSD'};
for k = 1:2
    figure(v_figures(k));
    v_legends(k)=legend(ch_algorithms);
    set(v_legends(k),'FontSize',8);
    xlabel ('Time, t (samples)')
    ylabel (ch_ylabels{k})
    F(k) = GFigure.captureCurrentFigure();
end
figure(v_figures(3));
F(3) = GFigure.captureCurrentFigure();