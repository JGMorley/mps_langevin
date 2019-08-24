function ddAncl12_H = calculate_ddAncl12_H( n, N, H, mps_LCF, mps_RCF )
    %% Calculate d/d(An.l^12) <H>

    if N==1
        A = mps_LCF{2}{1};
        ddAncl12_H = (H*A).';
        return
    end
    
    % 1. Make tensorList =
    %    {AL1,...,ALN,H,conj(AL1),...,conj(AL[n-1]),conj(AR[n+1]),...conj(ARN)}
    
    tensorList = cell([1 2*N+1]);
    for k=1:n-1
        tensorList{k} = mps_LCF{k+1}{1};           % AL_k
        tensorList{N+1+k} = conj(mps_LCF{k+1}{1}); % conj(AL_k)
    end
    
    tensorList{n} = mps_LCF{n+1}{1};
    tensorList{N+1} = H;
    
    for k=n+1:N
        tensorList{k} = mps_LCF{k+1}{1};       % AL_k
        tensorList{N+1+k} = conj(mps_RCF{k+1}{1}); % conj(AR_k)
    end
    
    tensorList(N+1+n) = [];
    
    % 2. Use calculateC_FC to get legLinksList
    
    [~,~,legLinksList] = calculateC_FC(n, mps_LCF, H, N, 'True');
    
    % 3. Contract and output
    
    ddAncl12_H = ncon(tensorList, legLinksList);
end