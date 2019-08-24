function plot_Schmidts( tdvpOutput )
    %% Plot Schmidt values of a tdvpIntegrator samples output
    nSamples = length(tdvpOutput);
    N = length(tdvpOutput{1,1}) - 1;
    T = tdvpOutput{end,2};
  
  
    %% Plot Schmidt coeffs
    figure
    progressWindow = waitbar(0, 'plotting Schmidt coeffs...');
    for n=1:N-1
        subplot(2,ceil((N-1)/2),n)
        xlim([0,T])
        grid on
        hold on
        % first sample
        samplek = tdvpOutput{1,1};
        lambdak = sqrt(samplek{n+1}{2});
        scatter(zeros([1, length(diag(lambdak))]),diag(lambdak),'.');

        for k=2:nSamples
            % find nth sites Schmidt coeffs and plot
            samplek = tdvpOutput{k,1};
            tk = tdvpOutput{k,2};
            lambdak = sqrt(samplek{n+1}{2});
            scatter(tk*ones([1, length(diag(lambdak))]),diag(lambdak),'.');        

            waitbar((n*nSamples + k)/(N*nSamples))
        end
    end
    close(progressWindow)

end

