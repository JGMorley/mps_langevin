function fidelity = fidelity_mps( mps1, mps2, abs_value )
    %% calculate the fidelity between two mps objects
    %    fidelity between state |1> and |2> defined as |<1|2>|
    %  optional argument abs_value determines whether |.| is taken
    
    if nargin==2, abs_value = true; end
    
    N = length(mps1)-1;
    if length(mps2)~=(N+1)
        error('input mps of different lengths')
    end
    
    if N==1
        A1 = mps1{2}{1};
        A2 = mps2{2}{1};
        if ~isequal(size(A1),size(A2)), A2 = A2.'; end
        if size(A1,1)==1
            fidelity = abs(A2*A1');
        else
            fidelity = abs(A2'*A1);
        end
        return
    end

    % contract from right
    for n=N:-1:2
            A1m = mps1{n+1}{1};
            A2m = mps2{n+1}{1};
            if n==N
                R = A1m*A2m';
            else
                R = ncon({A1m,R,conj(A2m)},{[-1 1 3],[1 2],[-2 2 3]});
            end
    end
    
    % first site
    A1 = mps1{2}{1};
    A2 = mps2{2}{1};
    fidelity = trace(A1*R*A2');
    if abs_value
        fidelity = abs(fidelity);
    end
end