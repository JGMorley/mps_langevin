function plot_mps_series_data(spinDimList,time_series,entropy_series,local_XYZ_series,...
                              H_given,energy_series,gs_fidelity_series,energy_gs,mps_gs,...
                              spin_operators)

    if nargin==4, H_given=false; end
                          
    N = length(spinDimList);
                        
    tFinal = time_series(end);
    tStart = time_series(1);
    % XYZ series
    figure('Position',[500 400 1000 400])

    nRows = floor(sqrt(N));
    nCols = ceil(N/nRows);
    for n=1:N
        subplot(nRows,nCols,n)
        hold on

        title(strcat('n=',string(n)))

        nth_idcs = (3*(n-1)+1):(3*(n-1)+3);
        plot(time_series,real(local_XYZ_series(nth_idcs(1),:)),'.b') % <Xn>
        plot(time_series,real(local_XYZ_series(nth_idcs(2),:)),'.r') % <Yn>
        plot(time_series,real(local_XYZ_series(nth_idcs(3),:)),'.',...
                                 'Color',[0.8 0.6 0.2]) % <Zn>

        if H_given && ~isequal(gs_fidelity_series,'none')
            X = spin_operators{n}{1};
            Y = spin_operators{n}{2};
            Z = spin_operators{n}{3};

            gs_X = real(singleSiteExpectationLCF(X,n,mps_gs));
            gs_Y = real(singleSiteExpectationLCF(Y,n,mps_gs));
            gs_Z = real(singleSiteExpectationLCF(Z,n,mps_gs));

            plot([tStart abs(tFinal)],gs_X*[1 1],'--b')
            plot([tStart abs(tFinal)],gs_Y*[1 1],'--r')
            plot([tStart abs(tFinal)],gs_Z*[1 1],'--','Color',[0.8 0.6 0.2])
        end

        xlabel('t')
        xlim([tStart abs(tFinal)])
        legend('<S_x>','<S_y>','<S_z>')    
    end 

    % (optional) energy and ground state fidelity
    if H_given && ~isequal(gs_fidelity_series,'none')
        figure

        subplot(1,2,1)
        hold on
        scatter(time_series,energy_series,'.')
        plot([tStart abs(tFinal)], energy_gs*[1 1],'--m');
        xlabel('t')
        xlim([0 abs(tFinal)])
        ylabel('<H>')

        subplot(1,2,2)
        scatter(time_series,gs_fidelity_series,'.')
        xlabel('t')
        xlim([tStart abs(tFinal)])
        ylabel('1-|\langle\psi|g.s.\rangle|')
        ylim([0 1])
    elseif H_given
        figure
        plot(time_series,real(energy_series),'.')
        xlabel('t')
        xlim([0 abs(tFinal)]);
        ylabel('<H>')
    end

    if N==1
        % there is no entropy
        return
    end

    % entropy series

    figure('Position',[500 400 1000 400])

    nRows = floor(sqrt(N-1));
    nCols = ceil((N-1)/nRows);

    Dmax = maxBondDim(spinDimList);
    Smax = log2(Dmax);
    for k=1:N-1
        subplot(nRows,nCols,k)
        hold on

        if k==1,     appstr='st cut'; 
        elseif k==2, appstr='nd cut';
        elseif k==3, appstr='rd cut';
        else,        appstr='th cut';
        end
        title(strcat('S(',string(k),appstr,')'))

        scatter(time_series,entropy_series(k,:),'.')

        if H_given && ~isequal(gs_fidelity_series,'none')
            kth_schmidts_gs = diag(mps_gs{k+1}{2});
            kth_entropy_gs = - sum(kth_schmidts_gs.*log2(kth_schmidts_gs));
            plot([tStart abs(tFinal)], kth_entropy_gs*[1 1], '--m')
        end

        kDmax = maxBondDim(spinDimList,k);
        kSmax = log2(kDmax);
        plot([tStart abs(tFinal)], [kSmax kSmax],'-.');
        xlabel('t')
        xlim([tStart abs(tFinal)])
        ylim([0 1.1*Smax]);
    end   
end

