function VL = leftGaugeTangentBasis( sqrtleft, A, p, D, TOL )
    %% Find V_L, the null space of the A'*lambda tensor
    % As shown in e.g. 'Summary' in Haegeman et al, PRB 88, 075133 (2013)
    %
    %% Inputs
    % left: [DxD] matrix
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
    
    % calculate L-tensor   
    L = ncon({conj(A),sqrtleft}, {[1 -3 -2], [-1 1]});  %L_{i, s, j}
    L = permute(L, [3 2 1]); % L_{j, s, i}
    L = reshape(L, [D,p*D]); % L_{j,(i,s)}
    %L = round(L,TOL); % why is this here?
    
    % find null vectors
    VL = null(L);  % VL_{(i,s), k} (kth column of VL is kth null vector of L)
    VLmat = VL;    % for checking unitarity later
    dimNullSpace = size(VL,2); % = pD - rank(L)
    
    % reshape the left index back out into a bond and a spin index
    try
        VL = reshape(VL, [p, D, dimNullSpace]); % VL_{s, i, k}
    catch ME
        if strcmp(ME.identifier,'MATLAB:getReshapeDims:notSameNumel')
            msg = ['Incorrect number of elements when reshaping VL.',...
                   sprintf(' size(VL) = %s.',mat2str(size(VL))),...
                   sprintf(' Attempted reshape(VL,[%d,%d,%d]).',p,D,D*(p-1)),...
                   ' size(A) = ',mat2str(size(A)),', size(left = ',...
                   mat2str(size(left))];
            causeException = MException('MATLAB:leftGaugeTangentBasis:numel',...
                                        msg);
            ME = addCause(ME,causeException);
        end
        rethrow(ME)
    end
    VL = permute(VL, [2 3 1]);         % VL_{i, k, s}
    
    %% check that VL contracted with lambda, A gives zeros as expected
    should_be_zeros = ncon({conj(A), sqrtleft, VL}, {[1 -1 2], [3 1], [3 -2 2]});
    if ~isequal(round(should_be_zeros, TOL),zeros([D dimNullSpace]))
        err.identifier = 'leftGaugeTangentBasis:outputError1';
        err.message = 'Output tensor is incorrect to within TOL';
        error(err)
    end
    
    %% check that VL*VL' gives a full-rank identity matrix within TOL digits
    iden_nullSpace = round(VLmat'*VLmat, TOL);
    if ~isequal(iden_nullSpace, eye(dimNullSpace))
        err.identifier = 'leftGaugeTangentBasis:outputError2';
        err.message = 'Output tensor is not correctly unitary to within TOL';
        error(err)
    end
end