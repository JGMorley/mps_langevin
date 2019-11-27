function left_envts = find_left_environments_cell(mpsLCF,mpsRCF)
    % compute environments left of n, as <LCF(j)|RCF(j')>
    
    N = length(mpsLCF)-1;
    left_envts = cell([1 N]);
    left_envts{1} = 1.;

    for k=1:(N-1)
        ARk = mpsRCF{k+1}{1}; 
        ALk = mpsLCF{k+1}{1};   
        if k==1
            left = ALk'*ARk;
        else
            left = ncon({conj(ALk),left,ARk},{[1 -1 3],[1 2],[2 -2 3]});
        end
        left_envts{k+1} = left;
    end
end