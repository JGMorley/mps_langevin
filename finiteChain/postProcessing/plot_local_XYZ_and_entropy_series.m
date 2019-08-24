function plot_local_XYZ_and_entropy_series(spinDimList,time_series,local_XYZ_series,entropy_series)
  
    N = size(local_XYZ_series,1)/3;
    tFinal = time_series(end);

    %% plot!

    % XYZ series
    figure('Position',[500 400 1000 400])

    nRows = floor(sqrt(N));
    nCols = ceil(N/nRows);
    for n=1:N
        subplot(nRows,nCols,n)
        hold on

        title(strcat('n=',string(n)))

        nth_idcs = (3*(n-1)+1):(3*(n-1)+3);
        scatter(time_series,real(local_XYZ_series(nth_idcs(1),:)),'.b') % <Xn>
        scatter(time_series,real(local_XYZ_series(nth_idcs(2),:)),'.r') % <Yn>
        scatter(time_series,real(local_XYZ_series(nth_idcs(3),:)),'.',...
                                 'MarkerFaceColor',[0.2 0.3 0.1]) % <Zn>
                             
        xlabel('t')
        xlim([0 abs(tFinal)])
        legend('<S_x>','<S_y>','<S_z>')    
    end 
    
    if N==1
        % there is no entropy
        return
    end

    if nargin==4
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

            kDmax = maxBondDim(spinDimList,k);
            kSmax = log2(kDmax);
            plot([0 abs(tFinal)], [kSmax kSmax],'-.');
            xlabel('t')
            xlim([0 abs(tFinal)])
            ylim([0 1.1*Smax]);
        end   
    end
end

