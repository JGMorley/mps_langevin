function mps = randmps( N, DMAX, spinDimList )
    %% randmps Returns MPS with elements populated by rand()
    %    mps = randmps(N,DMAX,spinDimList)
    %        Randomized N-system MPS with maximum bond dimension DMAX

    ROUND_DIGITS = 5;
    
    if nargin == 3
        if ~isequal(size(spinDimList), [1 N])
            err.identifier = 'finiteChain:randmps:invalidSpinDimList';
            err.message = 'spinDimList has wrong dimensions';
            error(err)
        end
    else
        spinDimList = randi([2 10],[1 N]);
    end
    
    mps = cell(1,N);
    mps{1} = rand() + 1i*rand();
    
    d = spinDimList(1);
    D = min([spinDimList(1) prod(spinDimList(2:N)) DMAX]);
    dims_1 = [d D];
    mps{2}{1} = round(rand(dims_1) + 1i*rand(dims_1), ROUND_DIGITS);
    mps{2}{2} = round(rand() + 1i*rand(), ROUND_DIGITS);
    
    if N == 1
        return
    end
    
    for k=2:N-1
        d = spinDimList(k);
        D_left  = min([prod(spinDimList(1:k-1)) prod(spinDimList(k:N)) DMAX]);
        D_right = min([prod(spinDimList(1:k)) prod(spinDimList(k+1:N)) DMAX]);
        dims_k = [D_left D_right d];
        mps{k+1}{1} = round(rand(dims_k) + 1i*rand(dims_k),ROUND_DIGITS);
        mps{k+1}{2} = round(rand() + 1i*rand(),ROUND_DIGITS);
    end
    
    d = spinDimList(N);
    D = min([prod(spinDimList(1:N-1)) spinDimList(N) DMAX]);
    dims_N = [D d];
    mps{N+1}{1} = round(rand(dims_N) + 1i*rand(dims_N),ROUND_DIGITS);  
    mps{N+1}{2} = round(rand() + 1i*rand(),ROUND_DIGITS);
end