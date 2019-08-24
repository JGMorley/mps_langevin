function contracted = nsp_exptn_contraction( n, N, ddAnc_O, nsp )
    %% Perform an often-needed contraction between null space projector and
    %  d/dconj(An) <O>, where O is some (Hermitian) operator
    
    if N==1
        contracted = (ddAnc_O*nsp).';
        return
    end
    
    if n==1
        contracted = ncon({ddAnc_O,nsp},{[1 -2],[1 -1]});
    elseif n==N
        contracted = ncon({ddAnc_O,nsp},{[1 2],[1 2 -1 -2]});
    else
        contracted = ncon({ddAnc_O,nsp},{[1 -2 2],[1 2 -1 -3]});
    end
    
end

%% As a picture
%                ____________________________
%               |                            |
%               |       d/dconj(An) <O>      |
% contracted =  |                            |      ,
%               |     __________________     |
%               |    |  _|__            |    |
%               |    |-|   |         ---|    |
%               |____| |nsp|        |   |____|
%                      |   |        |
%                (-1)--|___|         -----(-2
%                        |
%                      (-3)
%
% where
%
%           (-2)                    (-2)
%          ___|___                    |
%   (-1)--|        |        (-1)-conj(VL)--|
%         |  nsp   |   =                   |
%   (-3)--|________|        (-3)------VL---|
%             |                       |
%            (-4)                    (-4)            ,
%
% and index ordering for d/dconj(An)<O> follows from order for An.