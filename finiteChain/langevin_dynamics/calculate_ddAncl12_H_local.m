function ddAncl12_H = calculate_ddAncl12_H_local( ...
                                    n, N, Hcell, mps_LCF, mps_RCF, left_envts )
    %% Calculate *relevant* d/d(An.l^12) <H> using cell format for H 
    %  NB *relevant* here means terms that don't disappear on contraction
    %  with null vectors V
    if nargin==5
        left_envts = find_left_environments_cell(mps_LCF,mps_RCF); 
    end
    
    if N==1
        A = mps_RCF{2}{1};
        H = Hcell{1}{1};
        ddAncl12_H = (H*A).';
        return
    end
    
    if N==2
        if n==1
            O1  = Hcell{1}{1};
            O12 = Hcell{2}{1};
            AR1 = mps_RCF{2}{1};
            AR2 = mps_RCF{3}{1};
            
            singlesite = O1*AR1;% + AR1*AR2*Hcell{1}{2}*AR2'; 
            % ^^ second term dissapears with VV' contraction
            doublesite = ncon({O12,AR1,AR2,conj(AR2)},...
                              {[-1 4 1 3],[1 2],[2 3],[-2 4]});
        elseif n==2
            O1  = Hcell{1}{1};
            O2  = Hcell{1}{2};
            O12 = Hcell{2}{1};
            AL1 = mps_LCF{2}{1};
            AR1 = mps_RCF{2}{1};
            AR2 = mps_RCF{3}{:};
            
            singlesite = AL1'*O1*AR1*AR2 + AL1'*AR1*AR2*O2.';
            doublesite = ncon({O12,conj(AL1),AR1,AR2},...
                              {[1 -2 2 4],[1 -1],[2 3],[3 4]});
        end
        
        ddAncl12_H = singlesite + doublesite;
        return
    end

    % N>2:
    singlesite = Ngeq3_singlesite(Hcell,mps_LCF,mps_RCF,left_envts,n); 
    doublesite = Ngeq3_doublesite(Hcell,mps_LCF,mps_RCF,left_envts,n);
    ddAncl12_H = singlesite + doublesite;
end


function single_site = Ngeq3_singlesite(Hcell,mpsLCF,mpsRCF,left_envts,n)
    % calculate single-site H terms for mps with N>=3
    N = length(mpsLCF)-1;
    
    % n=1 as a special case
    if n==1
        O1 = Hcell{1}{1};
        AR1 = mpsRCF{2}{1};
        
        single_site = O1*AR1;
        return
    end
    
    % O1,...,O_{n-1} terms
    single_site = zeros(size(mpsLCF{n+1}{1}));
    for m=1:n-1
        Om = Hcell{1}{m};
        
        % construct left end
        ARm = mpsRCF{m+1}{1};
        ALm = mpsLCF{m+1}{1};
        if m==1
            left_end = ALm'*Om*ARm;
        else
            left_end = ncon({conj(ALm),Om,left_envts{m},ARm},...
                            {[1 -1 2],[2 3],[1 4],[4 -2 3]});
        end
        
        % contract left to right 
        for k=(m+1):(n-1)
            ARk = mpsRCF{k+1}{1};
            ALk = mpsLCF{k+1}{1};
            left_end = ncon({left_end,conj(ALk),ARk},...
                            {[1 3],[1 -1 2],[3 -2 2]});
        end
        
        % contract with right end
        ARn = mpsRCF{n+1}{1};
        if n==N
            single_site = single_site + left_end*ARn;
        else
            single_site = single_site + ncon({left_end,ARn},{[-1 1],[1 -2 -3]});
        end
    end
    
    % On term
    ARn  = mpsRCF{n+1}{1};
    On   = Hcell{1}{n};
    if n==N
        single_site = single_site + left_envts{n}*ARn*On.';
    else
        single_site = single_site + ncon({left_envts{n},ARn,On},...
                                         {[-1 1],[1 -2 2],[-3 2]});
    end
end


function double_site = Ngeq3_doublesite(Hcell,mpsLCF,mpsRCF,left_envts,n)
    % calculate double-site H terms for mps with N>=3
    N = length(mpsLCF)-1;
    
    double_site = zeros(size(mpsLCF{n+1}{1}));
    
    % n=1 special case
    if n==1
        O12 = Hcell{2}{1};
        AR1 = mpsRCF{2}{1};
        AR2 = mpsRCF{3}{1};
        AL2 = mpsLCF{3}{1};
        double_site = ncon({AR1,AR2,O12,conj(AR2)},...
                           {[1 2],[2 4 3],[-1 5 1 3],[-2 4 5]});
        return
    end 
    
    % n>1
    for m=1:(n-2)
        % construct left end
        Ommp = Hcell{2}{m};
        ARm = mpsRCF{m+1}{1};
        ARmp = mpsRCF{m+2}{1};
        ALm = mpsLCF{m+1}{1};
        ALmp = mpsLCF{m+2}{1};
        if m==1
            left_end = ncon({ARm,ARmp,Ommp,conj(ALm),conj(ALmp)},...
                            {[4 5],[5 -2 6],[1 3 4 6],[1 2],[2 -1 3]});
        else
            left_end = ncon({left_envts{m},ARm,ARmp,Ommp,conj(ALm),conj(ALmp)},...
                            {[2 5],[5 7 6],[7 -2 8],[4 3 6 8],[2 1 4],[1 -1 3]});
        end
        
        % contract left to right
        for k=(m+2):(n-1)
            ARk = mpsRCF{k+1}{1};
            ALk = mpsLCF{k+1}{1};
            left_end = ncon({left_end,ARk,conj(ALk)},{[1 2],[2 -2 3],[1 -1 3]});
        end
        
        % contract with right end
        ARn = mpsRCF{n+1}{1};
        if n==N
            double_site = double_site + left_end*ARn;
        else
            double_site = double_site + ncon({left_end,ARn},{[-1 1],[1 -2 -3]});
        end
    end
    
    %
    m=n-1;
    ARm  = mpsRCF{m+1}{1};
    ARmp = mpsRCF{m+2}{1};
    ALm  = mpsLCF{m+1}{1};
    Ommp = Hcell{2}{m};
    if m==1
        double_site = double_site + ...
                      ncon({conj(ALm),Ommp,ARm,ARmp},...
                           {[1 -1],[1 -3 2 4],[2 3],[3 -2 4]});
    elseif m==(N-1)
        double_site = double_site + ...
                      ncon({conj(ALm),left_envts{m},Ommp,ARm,ARmp},...
                           {[1 -1 2],[1 3],[2 -2 4 6],[3 5 4],[5 6]});
    else
        double_site = double_site + ...
                      ncon({conj(ALm),left_envts{m},Ommp,ARm,ARmp},...
                           {[1 -1 2],[1 3],[2 -3 4 6],[3 5 4],[5 -2 6]});
    end
    
    if n==N, return, end
    %
    m=n;
    ARm  = mpsRCF{m+1}{1};
    ARmp = mpsRCF{m+2}{1};
    Ommp = Hcell{2}{m};
    if m==(N-1)
        double_site = double_site + ...
                      ncon({left_envts{m},ARm,ARmp,Ommp,conj(ARmp)},...
                           {[-1 1],[1 3 2],[3 4],[-3 5 2 4],[-2 5]});
    else
        double_site = double_site + ...
                      ncon({left_envts{m},ARm,ARmp,Ommp,conj(ARmp)},...
                           {[-1 1],[1 3 2],[3 5 4],[-3 6 2 4],[-2 5 6]});
    end
end