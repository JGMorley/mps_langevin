function VR = rightGaugeTangentBasis( sqrtright, A, p, D, TOL )
    %% Find V_R, the null space of the A'*lambda tensor
    % As shown in e.g. 'Summary' in Haegeman et al, PRB 88, 075133 (2013)
    %
    %% Inputs
    % right: [DxD] matrix
    %      Dominant left eigenvector of the transfer matrix E
    % A: [DxDxp] tensor
    %      MPS tensor normalized st dominant eigenvalue of transfer matrix = 1
    % p: integer
    %      (optional) physical index dimension 
    % D: integer
    %      (optional) bond index dimension
    %%
    
    % check number of inputs
    switch nargin
        case 5
        case 4
            TOL = 6;
        case 3
            [~,D,~] = size(A);
            TOL = 6;
        case 2
            [~,D,p] = size(A);
            TOL = 6;
        otherwise
            err.identifier = 'leftGaugeTangentBasis:WrongNoInputs';
            err.message = 'Wrong number of input arguments';
            error(err);
    end
    
    % calculate R-tensor   
    R = ncon({A,sqrtright}, {[-1 1 -3], [1 -2]});
    R= permute(R,[1 3 2]); 
    R = reshape(R, [D,p*D]);
    
    % find null vectors
    VR = null(R);

    dimNullSpace = size(VR,2);
    
    VR = reshape(VR, [p, D, dimNullSpace]); % VL_{s, i, k}

    VR = permute(conj(VR),[3,2,1]);         % VL_{i, k, s}
    

end