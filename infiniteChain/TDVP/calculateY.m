function Y = calculateY(A,sqrtl,sqrtr,VL,VR,H)
    % This function calculates Y which is used to increase
    % the bond dimension of an MPS in inverse-free TDVP.
    %% Diagram:
    %             5
    %   2 __--A-------A--__ 8
    %    |    |4     6|    |
    %  sqrtl  ====H====   sqrtr
    %    |_   |3     7|   _|
    %    1 --VL-_   _-VR-- 9
    %            | |
    %           _| |_
    %         (-1)  (-2)
    %    

    Y=ncon({sqrtl,A,      A,      H,        sqrtr,conj(VL),conj(VR)},...
           {[1 2],[2 5 4],[5 8 6],[3 7 4 6],[8 9],[1 -1 3],[-2 9 7]});
end