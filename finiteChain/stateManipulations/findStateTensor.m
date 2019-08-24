function Psi = findStateTensor( mpsIn )
    %% findStateTensor Contract MPS tensors to find state tensor 
    %    Psi = findStateTensor(mps)
    %        Finds Psi, whose elements are given by Psi(s1,...) = <s1,...|psi>

    checkValidMPS(mpsIn);
    
    N = length(mpsIn) - 1;
    
    if N==1
        Psi = mpsIn{2}{1}; 
        return
    end
    
    tensorList = cell([1,N]);
    legLinksList = cell([1,N]);
    
    tensorList{1} = mpsIn{2}{1};
    tensorList{N} = mpsIn{N+1}{1};
    
    legLinksList{1} = [-1 1];
    legLinksList{N} = [N-1, -N];
    
    for n=2:N-1
        legLinksList{n} = [n-1, n, -n];
        tensorList{n} = mpsIn{n+1}{1};
    end
    
    Psi = ncon(tensorList, legLinksList);
end

%% As a picture:
%
%  ____________________________       
% |            Psi             |     
% |____________________________|  =   A--A--A--  ...    ... --A--A--A 
%  |  |  |  ...     ... |  |  |       |  |  |                 |  |  |