function [ A_new, right_new, U, Lambda ] = tangentBasis_unitaryGauge( A, right )
    %%TANGENTBASIS_UNITARYGAUGE puts left canonical form A in the unitary gauge
    %
    % > U(:,:,:,1) = A_new
    % > U(:,:,:,j!=1) gives left-gauge tangent vectors
    % > U[(sigma i),(rho j)] is unitary
    
    % find physical and bond dimensions d, D
    [~,D,d] = size(A);
 
    % reshape A[i, j, sigma] to A[(sigma, i), j]
    Areshaped = permute(A,         [1 3 2]); % A[i, sigma, j]
    Areshaped = reshape(Areshaped, d*D, D);  % A[(sigma, i), j]
    
    % SVD on bond index 'j'
    [U,Lambda,V] = svd(Areshaped); % U[(sigma i),(rho j)], V[i,j]
    V = V';
    
    % Gauge transform using V
    A_new = ncon({V,A,V'}, {[-1 1], [1 2 -3], [2 -2]});
    right_new = V*right*V';

    % Pre-multiply U by V and reshape to a rank-4 tensor
    U = reshape(U, [D d D d]); % U[(sigma i),(rho j)] -> U[i sigma j rho]

    U = ncon({V, U}, {[-1 1], [1 -2 -3 -4]}); % U[i sigma j rho]
    U = permute(U, [1 3 2 4]);                % U[i j sigma rho]
end

% Check U(:,:,:,1) with 
% U(:,:,:,1) - Anew

% Check tangent space properties with
% ncon({U(:,:,:,2:end), conj(Anew)}, {[1 -1 2 -2], [1 -3 2]})

% Check U is unitary with
% U = permute(U,[1 3 2 4]);
% U = reshape(U, [d*D d*D]);
% U*U' - eye(d*D)
% U'*U - eye(d*D)