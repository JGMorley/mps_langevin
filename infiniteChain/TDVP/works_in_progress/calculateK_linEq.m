function K=calculateK(C,A,E,right,left)
    %% Calculate K, a [DxD] matrix
    % As shown in e.g. 'Summary' in Haegeman et al, PRB 88, 075133 (2013)
    %
    %% Inputs
    % C: [DxDxpxp] tensor
    %      Two-site Hamiltonian tensor contracted with two A-tensors
    % A: [DxDxp] tensor
    %      MPS tensor normalized st dominant eigenvalue of transfer matrix = 1
    % E: [D^2xD^2] matrix
    %      Transfer matrix
    % left: [DxD] matrix
    %      Dominant left eigenvector of the transfer matrix E
    %%
    
    D = size(A,1);
    
    % reshape right, left into vectors
    leftmat = left;
    
    right = reshape(right.',[],1);
    left  = reshape(left.' ,[],1);
    
    % Check MPS is normalized
    ERminusR = E*right - right;
    LEminusL = left.'*E - left.';
    
    right_condition = isequal(round(ERminusR,12), zeros([D^2,1]));
    left_condition  = isequal(round(LEminusL,12), zeros([1,D^2]));
    if ~right_condition || ~left_condition
        error('MPS appears to be unnormalized')
    end
    
    % Construct M and b in M.x = b
    S = right*left.';
    ImQEQ = E - S;
    I = eye(D^2);
    Q = I - S;
    Q = reshape(Q, [D D D D]);
    Q = permute(Q, [2 4 1 3]);
    
    M = ImQEQ.';
    b = ncon({leftmat, C, conj(A), conj(A), Q},....
             {[1 2], [1 7 3 5], [2 4 3], [4 6 5], [7 -1 6 -2]});
    b = reshape(b.', [D^2,1]);
    
    x = linsolve(M,b);
    
    K = reshape(x, [D D]).';
end

%% Another method (using iterative solver): - will be good for checking!
% K = bicgstab(M,b) solves M*K=b for K
% M = diag(ones(1,D)) - E + right*left
% b = (ncon({C,A',A',left},{[1 -1 2 3],[4 5 2],[5 -2 3],[1 4]})) + ...
% ((ncon({C,A',A',left,right},{[1 2 3 4],[5 6 3],[6 7 4],[1 5],[2 7]}))*left)

%% As a picture:
%  /  |--      /  |--|     |--|     |--
% |   |       |   |  |  C  |  |     |
% | K |     = | l |  |_____|  | E_P |  
% |   |       |   |   |   |   |     |
%  \  |--      \  |-- A - A --|_____|--