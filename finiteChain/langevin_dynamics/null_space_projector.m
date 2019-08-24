function nsp = null_space_projector( An, n, N, direct)
    %% Construct null space projector - see below for diagram (NB LCF assumed)
    %  direct = true:  Use nsp = II - Aconj(A)
    %  direct = false: Calculate VL first
    
    if nargin==3
        direct = true;
    end
    
    if N==1
        % We don't bother with the null space projector for N=1
        nsp = eye(size(An,1));
        return
    end
    
    if n==1
        if direct
            nsp = eye(size(An,1)) - conj(An)*An.';
        else
            VL = leftTangentSpace(n,An,1);
            nsp = conj(VL)*VL.';
        end
    elseif size(size(An),2) ~= 2
        % 1 < n < N
        if direct
            nsp = ncon({eye(size(An,1)), eye(size(An,3))},{[-1 -3],[-2 -4]});
            nsp = nsp - ncon({conj(An), An},{[-1 1 -2],[-3 1 -4]});
        else
            VL = leftTangentSpace(n,An,eye(size(An,1)));
            nsp = ncon({conj(VL), VL},{[-1 1 -2],[-3 1 -4]});
        end
    else
        % n = N
        if direct
            nsp = ncon({eye(size(An,1)), eye(size(An,2))},{[-1 -3],[-2 -4]});
            nsp = nsp - ncon({conj(An),An},{[-1 -2],[-3 -4]});
        else
            VL = leftTangentSpace(n,An,eye(size(An,1)));
            nsp = ncon({conj(VL),VL},{[-1 1 -2],[-3 1 -4]});            
        end
    end
end

%% As a diagram:
%
% nsp(i,sigma,j,rho) = sum_k VL(j,rho,k) conj(VL(i,sigma,k))
%
%                   ->       (-2)
%                              |
%                    (-1)-conj(VL)--|
%                                   |
%                    (-3)------VL---|
%                              |
%                            (-4)
%
%                    = I(i,j)I(sigma,rho) 
%                                        - \sum_k An(j,rho,k) conj(An)(i,isigma)