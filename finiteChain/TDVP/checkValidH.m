function checkValidH( H, mpsIn )
    %% Raises errors if H doesn't match mps
    %
    %% Inputs
    % H: tensor
    %   Hamiltonian tensor H[i1,...,iN,j1,....,jN] = <i1...|H|j1...>
    %    or, cell { {O1,O2,...,ON},{O12,O23,...}}
    % mpsIn: mps cell array
    %   state corresponding to H
    %%
    
     N = length(mpsIn) - 1;
     
     if N==1
         d = max(size(mpsIn{2}{1}));
         if iscell(H)
             assert(isempty(H{2}),'H{2} should be empty')
             assert(isequal(size(H{1}{1}),[d d]),'O1 has wrong dimensions')
         else
             assert(isequal(size(H),[d d]),'H has wrong dimensions')
         end
         return
     end
    
    if iscell(H)
        for n=1:N
            % check singlesite terms
            On = H{1}{n};
            if n==1
                dn = size(mpsIn{n+1}{1},1);
            elseif n==N
                dn = size(mpsIn{n+1}{1},2);
            else
                dn = size(mpsIn{n+1}{1},3);
            end
            errstr = sprintf('On has wrong dimensions, n=%i',n);
            assert(isequal(size(On),[dn dn]),errstr);
        end
        
        for n=1:(N-1)
            % check doublesite terms
            Onnp = H{2}{n};
            if n==1
                dn = size(mpsIn{n+1}{1},1);
            else
                dn = size(mpsIn{n+1}{1},3);
            end
            if (n+1)==N
                dnp = size(mpsIn{n+2}{1},2);
            else
                dnp = size(mpsIn{n+2}{1},3);
            end
            errstr = sprintf('Onnp has wrong dimensions, n=%i, np=%i',n,n+1);
            assert(isequal(size(Onnp),[dn dnp dn dnp]),errstr);
        end
    else
        % Check H has correct rank and dimensions
        if ~isequal(size(size(H)), [1 2*N])
            err.identifier = 'finiteChain:checkValidH:rankH';
            err.message = 'H has incorrect rank';
            error(err)
        end

        physDimList = zeros([1 N]);
        physDimList(1) = size(mpsIn{2}{1},1);
        for n=2:N
            dims = size(mpsIn{n+1}{1});
            physDimList(n) = dims(end);    
        end

        if ~isequal(size(H), [physDimList physDimList])
            err.identifier = 'finiteChain:checkValidH:dimH';
            err.message = 'H has incorrect dimensions';
            error(err)
        end    
    end
end