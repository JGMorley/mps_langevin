function [ sizes, keys ] = sizeGenericTransferMatrix( mpsIn, r )
    % Calculate size (largest abs eigenvalue) of finite-chain transfer
    % matrix used in calculating correlation functions: C(n,n+r) = v.M.w where
    % v and w are specific to the operators of interest and M is
    % state-specific.
    %   [sizes, key] = sizeGenericTransferMatrix(mpsIn) gives all sizes
    %   [sizes, key] = sizeGenericTransferMatrix(mpsIn) gives sizes at length r
    %
    % 'key' is a cell array of strings labelling the dimensions of the
    % sizes array. eg. for N=3 might have key={'1,2','2,3','1,3'}.
    %
    % NB mpsIn must be in LCF!
    
    N = length(mpsIn) - 1;
    
    if nargin==1 
        r = 'all'; 
    else
        if ~any(r==1:N)
           'r = ',r
           error('invalid value for r')
        end
    end
    
    genkeys = nargout==2;
    
    if strcmp(r,'all')
        nmax = N-1;
        sizes = zeros([1 N*(N-1)/2]); % length = sum_{nmax=1}^{N-1} nmax
    else
        nmax = N-r;
        sizes = zeros([1 nmax]);
    end
    
    if genkeys
        keys = cell([1 length(sizes)]);
    end
    
    running_idx = 1;
    
    for n1 = 1:nmax
        % find D_n1, R_n1
        D_n1 = size(mpsIn{n1+1}{1},2);
        R_n1 = mpsIn{n1+1}{2};
        R_n1 = reshape(R_n1.', [D_n1^2 1]);
        
        % create T
        T = ncon({eye(D_n1),eye(D_n1)},{[-2 -4],[-1 -3]});
        T = reshape(T, D_n1^2*[1 1]);
        
        for n2 = (n1+1):N
            % find L_n2
            An2 = mpsIn{n2+1}{1};
            [D_n2m1,D_n2,~] = size(An2);
            L_n2 = eye(D_n2m1);
            L_n2 = reshape(L_n2, [D_n2m1^2 1]);
            
            % find M
            M = T - R_n1*L_n2.';
            
            % calc max(abs(eig(M)) and write to sizes
            singular_values = svd(M);
            if strcmp(r,'all')
                sizes(running_idx) = singular_values(1);
                if genkeys
                    keys{running_idx} = strcat(num2str(n1),',',num2str(n2));
                end
                running_idx = running_idx + 1;
            elseif abs(r - (n2 - n1)) < 0.001
                sizes(running_idx) = singular_values(1);
                if genkeys
                    keys{running_idx} = strcat(num2str(n1),',',num2str(n2));
                end
                running_idx = running_idx + 1;
            end
            
            % write to keys if requested
            
            
            % postmultiply T to get new T
            Tn2 = ncon({An2,conj(An2)},{[-2 -4 1],[-1 -3 1]});
            Tn2 = reshape(Tn2, [D_n2m1^2 D_n2^2]);
            T = T*Tn2;
        end
    end
end

% % code to test:
% mps = randmps(6,8,[2 2 2 2 2 2]);
% mps = canonicalFormFC(mps,'LCF',true);
% 
% R1 = mps{2}{2};
% A2 = mps{3}{1}; A3 = mps{4}{1}; A4 = mps{5}{1};
% L5 = eye(size(A4,2));
% 
% T234 = ncon({A2,A3,A4,conj(A2),conj(A3),conj(A4)},...
%             {[-2 1 3],[1 2 4],[2 -4 5],[-1 6 3],[6 7 4],[7 -3 5]});
% R1L5 = ncon({R1,L5},{[-2 -1],[-3 -4]});
% 
% M = T234 - R1L5;
% M = reshape(M,[size(R1,1)^2 size(L5,1)^2]);
% 
% svals = svd(M);
% 
% fnout = sizeGenericTransferMatrix(mps,4);
% 
% diff = svals(1) - fnout(1)