function [left_env_n, right_env_n] = compute_env( mpsIn, n, RorL )
    % function for computing the left and right environment for the nth
    % site, making no assumptions about the gauge of mpsIn
    
    if nargin==3 && ~(strcmp(RorL,'R') || strcmp(RorL,'L') || strcmp(RorL,'B'))
        error('invalid RorL, use R or L or B(oth)')
    elseif nargin==2
        RorL = 'B'; % both
    end
    
    doL = false; doR = false;
    
    if strcmp(RorL,'L') || strcmp(RorL,'B'), doL = true; end
    if strcmp(RorL,'R') || strcmp(RorL,'B'), doR = true; end
    
    N = length(mpsIn) - 1;
    
    left_env_n = false; right_env_n = false;
    
    if doL
        % sweep left to right to find left_env_n
        for m=1:n-1
            Am = mpsIn{m+1}{1};
            if m==1
                L = Am'*Am;
            else
                L = ncon({conj(Am),L,Am},{[1 -1 3],[1 2],[2 -2 3]});
            end
        end
        
        if n==1, L = 1.; end
        left_env_n = L;
    end
    
    if doR
        % sweep right to left to find right_env_n
        for m=N:-1:n+1
            Am = mpsIn{m+1}{1};
            if m==N
                R = Am*Am';
            else
                R = ncon({Am,R,conj(Am)},{[-1 1 3],[1 2],[-2 2 3]});
            end
        end
        
        if n==N, R = 1.; end
        right_env_n = R;
    end
end