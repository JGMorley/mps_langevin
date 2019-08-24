function X = squareDecomposition(M,TOL,TYPE)
    %Decompose a Hermitian, positive definite, matrix M into it's squares X, X'
    %  X = squareDecomposition(M, TOL) is returned such that M = X*X'.
    %    M is checked for Hermiticity and positive definiteness to within a
    %    tolerance of TOL digits
    %
    %  if TYPE==1, find M s.t. M = X*X'
    %  if TYPE==2, find M s.t. M = X*X (only possible for Hermitian M !)
    if nargin<2
        TOL = 8;
        TYPE = 1;
    elseif nargin<3 || (TYPE~=1 && TYPE~=2)
        % TYPE is missing or an invalid value
        TYPE = 1;
    end
    
    if TYPE==2
        X = sqrtm(M);
        return
    end
    
    dim = size(M,1); % bond index dimension
    
    % Scale tolerance by largest element of M
    max_abs_el = abs(max(max(M)));
    TOL = TOL - int64(ceil(log10(max_abs_el)));
    UVTOL = TOL - 3; % U==V check needs less precision than 
                     % final decomposition check
    
    % First check H is Hermitian and positive definite
    diff = round(M - M', TOL);
    HisHermitian = isequal(diff, zeros([dim,dim]));
    if not(HisHermitian) 
        errstr = strcat('Input matrix is not Hermitian',...
                        sprintf('\nInput matrix M = '),mat2str(M),...
                        sprintf('\nM - hc(M) = '),mat2str(M - M'));
        err.message = errstr;
        err.identifier = 'squareDecomposition:notHermitian';
        warning(err.message);
        %error(err);
    end
    
    e = eig(M);
    clear skip_UV_check
    if not(all(e + abs(e))) && TYPE==1  % if not all positive
        warning('Input matrix not positive definite')
        skip_UV_check = 1; % in this case U!=V is ok: X*X' - H still small
    end
    % Use svd() to decompose H
    [U,D,V] = svd(M); % s.t. H = U*D*V' with U = V
    
    % Check U and V are equal
    diff = round(U - V, UVTOL);
    if not(isequal(diff, zeros([dim,dim]))) && exist('skip_UV_check','var')
        warning('U and V not equal in svd - consider lowering TOL')
        [sprintf('\nH = '),mat2str(M),sprintf('\nU = '),mat2str(U),...
         sprintf('\nD = '),mat2str(D),sprintf('\nV = '),mat2str(V),...
         sprintf('\nTOL = '),num2str(UVTOL)]
        %error('stopping');
    end
    
    % Absorb square roots of the diagonal D into U, V
    if TYPE==1
        X = U*sqrt(D);
        diff = round(X*X' - M, TOL);
    else
        X = zeros(dim);
        for k=1:dim
            X = X + sqrt(D(k,k)) * U(:,k) * V(:,k)';
        end
        diff = round(X*X - M, TOL);
    end
    
    
    % Check X*X' = H (or X*X = H as appropriate)
    if not(isequal(diff, zeros([dim,dim])))
        warning('final decomposition incorrect')
        %error('final decomposition incorrect')
    end