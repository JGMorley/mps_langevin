function [A_L_new,z_new,Dnew,comp_err] = decrease_bond_dim(A,comp_err_thresh)
    % This algorithm finds a new AL and AR and z 
    % which will be used in the time evolution algorithm
    % to find AL(t+dt) etc. 

    % Compression is done by effectively neglecting the smallest Schmidt
    % values.
    
    if nargin==1, comp_err_thresh=1e-8; end
    
    % Take first chi_new singular values, this can be done by pre- and
    % post-multiplying by [eye(chi_new) zeros([chi_new, chi_old-chi_new])]
    D = size(A,1);
    Dnew = D-1;
    
    mps = canonicalForm({1,A,'none'});
    z   = sqrt(mps{3});
    A   = mps{2};
    
    I0 = [eye(Dnew) zeros([Dnew, D-Dnew])];
    
    AL_comp = ncon({I0,A,I0.'},{[-1 1],[1 2 -3],[2 -2]});
    
    % For this approx'n err = || |new> - |old> ||^2 = sum of squared
    % neglected Schmidt values
    diagz = diag(z);
    comp_err = sum(diagz(Dnew+1:D).^2);
    if comp_err > comp_err_thresh
        A_L_new = A;
        Dnew = D;
        z_new = z;
    else
        A_L_new = AL_comp;
        z_new = diag(diagz(1:Dnew));
    end
end

