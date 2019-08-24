function V = tangentSpace(n, A, left, right, varargin)
    %% Finds tangent space of MPS tensor A (left gauge)
    %
    %% Inputs
    %  A: tensor
    %    MPS tensor
    %  left: matrix
    %    left-environment of A
    %  right: matrix
    %    right-environment of A
    %  method: optional, str
    %    'null' for Haegeman construction
    %    'svd' for unitary gauge contruction
    %
    %% Outputs
    %  V: tensor
    %    tangent space basis of A, arranged in columns
    %%
    
    TOL = 8;
    
    %% Parse input arguments    
    p = inputParser;
    addOptional(p,'method','null',@(x)any(validatestring(x,{'null','svd'})))
    parse(p,varargin{:})
    
    l12 = sqrtm(left);
    
    if p.Results.method == 'null'
        % Haegeman construction using null()
        if n==1
            % A = A[\sigma, i]
            V = null(A'); % V[\sigma,i]
            return
        elseif isequal(size(size(A)), [1 2])
            % n == N
            [D,d] = size(A);
            L = l12*conj(A);  % L[i,\sigma]
            L = reshape(L.',[1 d*D]); % L[(i,\sigma)]
            
            V = null(L);             % V[(i,\sigma),k]
            V = reshape(V,d,D,[]);   % V[\sigma,i,k)
            V = permute(V, [2 3 1]); % V[i,k,\sigma]
            
            %% check that VL contracted with lambda, A gives zeros as expected
            should_be_zeros = ncon({conj(A),sqrtm(left),V},...
                                   {[1 2],[3 1],[3 -1 2]});
            if ~isequal(round(should_be_zeros, TOL),zeros([size(V,2) 1]))
                err.identifier = 'finiteChain:tangentSpace:outputError1';
                err.message = 'Output tensor is incorrect to within TOL';
                error(err)
            end
            
            %% check that V*V' gives a full-rank identity matrix within TOL digits
            iden_nullSpace = round(ncon({V, conj(V)}, {[1 -1 2],[1 -2 2]}), TOL);
            if ~isequal(iden_nullSpace, eye(size(V,2)))
                err.identifier = 'leftGaugeTangentBasis:outputError2';
                err.message = 'Output tensor is not correctly unitary to within TOL';
                error(err)
            end
            return
        end
        
        [D1,D2,d] = size(A);
        
        L = ncon({conj(A),sqrtm(left)},{[1 -1 -3],[-2 1]}); % L[i,j,\sigma]
        L = reshape(L,[D2 d*D1]); % L[i,(\sigma j)]
        
        V = null(L); % V[(\sigma j),k], where dimension(k) = min(d*D1-D2,0)
        Vmat = V; % for checking!
        dimNullSpace = size(V,2);
        if dimNullSpace ~= d*D1 - D2
            err.identifier = 'finiteChain:tangentSpace:wrongNumberNullDims';
            err.message = ['Unexpected number of null dimensions found',...
                           ', suggests not everything is full rank'];
            error(err);
        end
        V = reshape(V,D1,d,dimNullSpace); % V[j,\sigma,k]
        V = permute(V, [1 3 2]); % V[j,k,\sigma]
    else
        % Unitary gauge transformation 
        % NB this changes the whole state, which will need re-canonicalizing!
    end
    
    %% check that VL contracted with lambda, A gives zeros as expected
    should_be_zeros = ncon({conj(A), sqrtm(left), V},...
                           {[1 -1 2], [3 1], [3 -2 2]});
    if ~isequal(round(should_be_zeros, 10),zeros([D2 dimNullSpace]))
        err.identifier = 'finiteChain:tangentSpace:outputError1';
        err.message = 'Output tensor is incorrect to within TOL';
        error(err)
    end
     
    %% check that V*V' gives a full-rank identity matrix within TOL digits
    iden_nullSpace = round(Vmat'*Vmat, 10);
    if ~isequal(iden_nullSpace, eye(dimNullSpace))
        err.identifier = 'leftGaugeTangentBasis:outputError2';
        err.message = 'Output tensor is not correctly unitary to within TOL';
        error(err)
    end
end