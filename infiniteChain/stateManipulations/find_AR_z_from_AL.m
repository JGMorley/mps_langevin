function [AR, z] = find_AR_z_from_AL(AL,z_input,d,D)
    %% find_AR_z_from_AL Finds AR, z from AL via successive RQ decompositions
    % z_input is an initial z which is used to start the process. 
    % Perhaps the z from the previous timestep should be used.
    lambda_approx = z_input;
    
    if isequal(lambda_approx,'none')
        lambda_approx = eye(D)/sqrt(D);
    end
        
    nLoops = 1e4;
    CONVERGENCE_THRESHOLD = 100*eps;
    lambda_n=cell(nLoops,1);
    AR_n=cell(nLoops,1);
    for k=1:nLoops
        
        % Form AC
        AC = ncon({AL, lambda_approx}, {[-1 1 -2],[1 -3]}); % AC[i,\sigma,j]
        AC=reshape(AC,[D,D,d]);        AC = reshape(AC, [D, d*D]); % AC[i,(j,\sigma)]

        [r,q] = rq(AC);
        % Fix the signs of r and q so lambda is positive
        fix1=zeros(D,1);
        for n=1:D
            fix1(n)=sign(sign(r(n,n))+0.1);
        end
        
        fix = diag([fix1;1*ones([(d-1)*D,1])]);
        r=r*fix;
        q=fix*q;
        
        % Project out redundant dimensions
        lambda_new=r(1:D,1:D);
        Q=q';
        AR=reshape(q(1:D,:),[D,d,D]);
        AR=permute(AR,[1,3,2]);
        lambda_n{k}=lambda_new;
        AR_n{k}=AR;
        
        % Test for convergence
        diff = lambda_new - lambda_approx;
        if max(abs(diff(:))) < CONVERGENCE_THRESHOLD
            z = lambda_new/sqrt(trace(lambda_new*lambda_new'));           
            return
        elseif k>2 && trace(lambda_new)>0 
            diff2=lambda_new-lambda_n{k-2};
            if max(abs(diff2(:))) < CONVERGENCE_THRESHOLD
                z = lambda_new/sqrt(trace(lambda_new*lambda_new'));
                return
            end          
        end
        
        % intialize next loop
        lambda_approx = lambda_new;
    end
    error(['Failed to converge within ',num2str(nLoops),' loops.'])
end


function [R Q]=rq(A)
    % function [R Q]=rq(A)
    % A [m x n] with m<=n
    % return R [m x n] triangular and
    % Q [n x n] unitary (i.e., Q'*Q=I)
    % such that A=R*Q
    % Author: Bruno Luong
    % Last Update: 04/Oct/2008

    [m n]=size(A);
    if m>n
        error('RQ: Number of rows must be smaller than column');
    end

    [Q R]=qr(flipud(A).');
    R=flipud(R.');
    R(:,1:m)=R(:,m:-1:1);
    Q=Q.';
    Q(1:m,:)=Q(m:-1:1,:);

end