function ddAnl12_Fmk = calculate_ddAnl12_Fmk( n, m, Fmk, mps_LCF)
    %% Calculate d/d(An.l^12) <F_{m,k}> (kth operator on mth site) 
    
    % 0. Handle N=1 separately
    N = length(mps_LCF) - 1;
    if N==1
        A = mps_LCF{2}{1};
        ddAnl12_Fmk = A'*Fmk;
        return
    end
    
    
    % 1. If m > n, set to zero (since using LCF and left gauge for d/dtAn)
    if m > n
        ddAnl12_Fmk = {};
        return
    end
    
    % 2. Extract mth and nth MPS
    ALm = mps_LCF{m+1}{1};
    ALn = mps_LCF{n+1}{1};  
    Ln12 = sqrt(mps_LCF{n+1}{2});
    
    % 3. Calculate the transfer matrix for sites p: m < p < n
    E = ncon({eye(size(ALm,2)), eye(size(ALm,2))}, {[-1 -2], [-3 -4]});
    for p = (m+1):(n-1)
        ALp = mps_LCF{p+1}{1};
        E = ncon({E,ALp,conj(ALp)},{[-1 1 -3 2],[1 -2 3],[2 -4 3]});
    end
    
    % 4. Perform final contractions, special cases for end sites
    
    if m==1
        if n==1
            ddAnl12_Fmk = ncon({Fmk,conj(ALm),conj(Ln12)},...
                            {[1 -1],[1 2],[-2 2]});
        elseif n<N
            ddAnl12_Fmk = ncon({ALm,Fmk,conj(ALm),E,conj(ALn),conj(Ln12)},...
                            {[1 2],[3 1],[3 4],[2 -1 4 5],[5 6 -3],[6 -2]});
        else
            ddAnl12_Fmk = ncon({ALm,Fmk,conj(ALm),E,conj(ALn)},...
                            {[1 2],[3 1],[3 4],[2 -1 4 5],[5 -2]});
        end
    elseif m<N
        if n==m
            ddAnl12_Fmk = ncon({Fmk,conj(ALn),Ln12},{[1 -3],[-1 2 1],[2 -2]});
        elseif n<N
            ddAnl12_Fmk = ncon({ALm,Fmk,conj(ALm),E,conj(ALn),Ln12},...
                            {[1 2 4],[5 4],[1 3 5],[2 -1 3 6],[6 7 -3],[7 -2]});
        else % n==N
            ddAnl12_Fmk = ncon({ALm,Fmk,conj(ALm),E,conj(ALn)*Ln12},...
                             {[1 4 3],[2 3],[1 5 2],[4 -1 5 6],[6 -2]});
        end 
    else
        % m = n = N
        ddAnl12_Fmk = ncon({Fmk,conj(ALn)},{[1 -2],[-1 1]});
    end

end

%% As a diagram:
%                            ________________
%                 ----ALm---|                |---[-1]      [-2]------
%                |    _|_   |                |        [-3]           |
%                |   |   |  |                |          |            |
%  ddAn_Fmk =    |   |Fmk|  | E(m+1,...,n-1) |          |            |
%                |   |___|  |  (TM of ALs)   |          |            |
%                |     |    |                |          |            |
%                 ----ALm---|________________|---------ALn--(Ln^12)--
%                 (mth site)                       (nth site)    
%
% where E is the transfer operator for sites (m+1) to (n-1) with index
% ordering as:
%                  ________________
%           [-1]--|                |--[-2]
%                 |                |
%                 |                |
%  E =            | E(m+1,...,n-1) |
%                 |                |
%                 |                |
%           [-3]--|________________|--[4]