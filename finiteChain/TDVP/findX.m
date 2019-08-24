function X = findX(n,N,C,left,right,V)
    %% Contract to find X, which parameterizes the change in MPS tensor dA
    % NB Currently written for the left-tangentSpace form
    %
    %% Inputs
    %  n : integer
    %    site index
    %  N:  integer
    %    total number sites
    %  C: tensor
    %    conj(A) away from a full <psi|H|psi> contraction
    %  left: matrix
    %  right: matrix
    %    left- and right- environments of MPS tensor
    %  V: tensor
    %    left tangent space
    %
    %% Outputs
    %  X: matrix
    %    matrix that parameterizes dA. Has dimensions [dimNullSpace Dright]
    %%
    
    lm12c = conj(inv(sqrtm(left)));
    rm12c = conj(inv(sqrtm(right)));
    
    if N==1
        X = ncon({conj(V),rm12c,C.'},{[1 -1],[-2 2],[1 2]});
        return
    end
    
    if n==1
        X = ncon({conj(V),rm12c,C},{[1 -1],[-2 2],[1 2]});
    elseif n == N
        X = ncon({conj(V),C,lm12c},{[1 -1 2],[3 2],[3 1]});
    else
        X = ncon({conj(V),lm12c,C,rm12c},{[1 -1 2],[3 1],[3 4 2],[4 -2]});
    end
end

%% As a picture:
%
%  ==X--  =  
%             O-- ... --------O--------------O--- ... --O
%             |_______________|______________|__________|
%            |___________________________________________|
%             |               |              |          |
%             O-- ... --lm12--V==   --rm12---O--- ... --O
%
%         =
%               --------C-------------
%              |        |             |           
%              |        |             |
%               --lm12--V==   --rm12--
%
% where lm12 = left^(-1/2), rm12 = right^(-1/2), and all tensors but C are
% elementwise complex conjugated
%%