function mpsOut = projectOutNullDims( mpsIn )
    %% projectOutNullDims Truncate null bond dimensions of mpsIn
    %    mpsOUt = projectOutNullDims(mpsIn)
    %        Gauge transformed MPS with null bond dimensions removed

    mpsOut = mpsIn;
    N = length(mpsOut) - 1;
    
    for n=1:N-1
        Ln = mpsOut{n+1}{2};
        D_old = size(Ln,1);
        D_new = rank(Ln);
        if D_new~=size(Ln,1)
            if ~issorted(-diag(Ln))
                err.identifier = 'finiteChain:projectOutNullDims:notSorted';
                err.message = 'Singular values not sorted';
                error(err);
            end
            P = [eye(D_new), zeros(D_new,(D_old - D_new))];
            
            An = mpsOut{n+1}{1};
            Anp = mpsOut{n+2}{1};
            mpsOut{n+1}{1} = ncon({An,P'},{[-1 1 -3], [1 -2]});
            mpsOut{n+2}{1} = ncon({P,Anp},{[-1 1], [1 -2 -3]});
            mpsOut{n+1}{2} = P*Ln*P';
        end
    end
end

