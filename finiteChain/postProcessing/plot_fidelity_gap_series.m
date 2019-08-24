function fidelity_gap_series = plot_fidelity_gap_series(mps_series_1,mps_series_2,time_series,make_plots)
    %% plot series of fidelity gaps := 1 - |<1|2>|
    
    assert(length(mps_series_1)==length(mps_series_2))
    assert(length(mps_series_1)==length(time_series))
    
    if nargin==3
        make_plots=true;
    end
    
    tFinal=time_series(end);
    nSamples = length(time_series);
    
    fidelity_gap_series = zeros([1 length(time_series)]);
    for iSample=1:nSamples
        fidelity_gap_series(iSample) = 1 - fidelity_mps(mps_series_1{iSample}, mps_series_2{iSample});
    end
    
    if make_plots
        figure
        plot(time_series,fidelity_gap_series)
        ylabel('$1-|\langle\psi_1(t)|\psi_2(t)\rangle|$','interpreter','latex')
        xlabel('t')
        xlim([0 tFinal])
    end
    
    
    if nargout==0
        fidelity_gap_series=NaN;
    end
end

