function plot_tdvp_output(output,H)
    % plot various things about the mps_series
    if nargin==1
        H_given = false;
    else
        H_given = true;
    end
    
    mps_series = {output{:,1}}.';
    time_series = [output{:,2}];
    
    A0 = mps_series{1};
    tFinal = time_series(end);
    
    [~,D,d] = size(A0);
    [Sx,Sy,Sz] = spinMatrices((d-1)/2);

    % calculate expectations
    
    nSamples = length(mps_series);
    local_XYZ_series = zeros([3 nSamples]);  
    entropy_series = zeros([1, nSamples]);
    D_series = zeros([1 nSamples]);
    schmidtCoeffs = zeros([D nSamples]);
    if H_given
        energy_series = zeros([1 nSamples]);
    end

    for n = 1:nSamples
        % Find A, R, L
        An = mps_series{n};
        [~,Rn,Ln] = normalizeMPS(An);
        
        D = size(An,1);
        D_series(n) = D;

        % calculate expectations
        rho = ncon({An, Rn, conj(An), Ln}, {[1 2 -1], [2 3], [4 3 -2], [1 4]});
        energy_series(n) = ncon({An,An,Rn,conj(An),conj(An),Ln,H},...
                         {[1 2 5],[2 3 4],[3 6],[8 6 7],[9 8 10],[1 9],[10 7 5 4]});

        expSx = ncon({An,Rn,conj(An),Ln,Sx},{[2 3 6],[3 4],[1 4 5],[2 1],[5 6]});
        expSy = ncon({An,Rn,conj(An),Ln,Sy},{[2 3 6],[3 4],[1 4 5],[2 1],[5 6]});
        expSz = ncon({An,Rn,conj(An),Ln,Sz},{[2 3 6],[3 4],[1 4 5],[2 1],[5 6]});
        local_XYZ_series(:,n) = [expSx expSy expSz].';

        evalues = eig(rho);
        S = 0;
        for k = 1:size(evalues,1)
            S = S - evalues(k)*log(evalues(k));
        end
        entropy_series(n) = S;  

        mpsn = canonicalForm({1,An,0});
        schmidts = diag(mpsn{3});
        schmidtCoeffs(1:length(schmidts),n) = schmidts;
    end

    %% Check imaginary parts
    lgst_imag_XYZ = max(max(abs(imag(local_XYZ_series))));
    if lgst_imag_XYZ > 1000*eps
        warning('local_XYZ_series has imaginary part')
    end
    local_XYZ_series = real(local_XYZ_series);
    
    lgst_imag_entropy = max(abs(imag(entropy_series)));
    if lgst_imag_entropy > 1000*eps
        warning('entropy_series has imaginary part')
    end
    entropy_series = real(entropy_series);
    
    if H_given
        lgst_energy = max(abs(imag(energy_series)));
        if lgst_energy > 1000*eps
            warning('energy_series has imaginary part')
        end
        energy_series = real(energy_series);
    end
    
    %% plot!

    % XYZ series
    figure('Position',[500 400 1000 400])

    % X
    subplot(2,3,1)
    scatter(time_series,local_XYZ_series(1,:),'.')
    xlabel('t')
    xlim([0 abs(tFinal)])
    ylabel('<Sx>')    
    ylim([-0.5, 0.5])
    
    % Y
    subplot(2,3,2)
    scatter(time_series,local_XYZ_series(2,:),'.')
    xlabel('t')
    xlim([0 abs(tFinal)])
    ylabel('<Sy>')   
    ylim([-0.5, 0.5])
    
    % Z
    subplot(2,3,3)
    scatter(time_series,local_XYZ_series(3,:),'.')
    xlabel('t')
    xlim([0 abs(tFinal)])
    ylabel('<Sz>') 
    ylim([-0.5, 0.5])
    
    % S
    subplot(2,3,4)
    scatter(time_series,entropy_series,'.')
    xlabel('t')
    xlim([0 abs(tFinal)])
    ylabel('entropy')    
    
    % D
    subplot(2,3,5)
    plot(time_series,D_series,'.-')
    xlabel('t')
    xlim([0 abs(tFinal)])
    ylabel('D')    
    
    % Schmidt vals
    subplot(2,3,6)
    hold on
    for k=1:size(schmidtCoeffs,1)
        scatter(time_series,schmidtCoeffs(k,:),'.')
    end
    xlabel('t')
    xlim([0 abs(tFinal)])
    ylabel('Schmidt Coefficients')
    
    % (optional) E
    if H_given
        figure
        scatter(time_series,energy_series,'.')
        xlabel('t')
        xlim([0 abs(tFinal)])
        ylabel('<H>')
    end