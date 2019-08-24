function mutual_info_series = plot_mutual_info(mps_series,time_series,show_plots,show_waitbar)
    %% Plot entropy along 4th cut and mutual information for 1-nth sites, N=9
    
    if nargin==2, show_plots=true; show_waitbar=false; end
    if nargin==3, show_waitbar=false; end
    
    N = length(mps_series{1}) - 1;
    
    nSamples = length(time_series);
    
    mutual_info_series = zeros([N,N,nSamples]);
    
    % get data
    if show_waitbar, h = waitbar(0,'Calculation progress...'); end
    for iSample = 1:nSamples
        mps = mps_series{iSample};
        
        % mutual info
        for n=1:N
            for m=1:N
                mutual_info_series(n,m,iSample) = mutual_information_mps(mps,n,m);
            end
        end
        if show_waitbar
            txt = ['Progress ',num2str(iSample),'/',num2str(nSamples),' samples'];
            waitbar(iSample/nSamples,h,txt)
        end
    end
    if show_waitbar, close(h), end
    
    %% plot
    
    if show_plots
        figure('Position',[100 300 1472 642])
        max_val = max(mutual_info_series(:));
        for n=1:N
            for m=n+1:N
                subplot(N,N,N*(n-1)+m-1)
                hold on
                plot(time_series,squeeze(real(mutual_info_series(n,m,:))))
                [max_val_n,max_val_idx] = max(real(mutual_info_series(n,m,:)));
                plot(time_series,squeeze(real(mutual_info_series(n,m,:)))/max_val_n,'r:')
                xlabel('t')
                xlim([min(time_series) max(time_series)])
                title(sprintf('I(%i,%i)',n,m))
                ylim([0 1])
            end
        end    
    end
end