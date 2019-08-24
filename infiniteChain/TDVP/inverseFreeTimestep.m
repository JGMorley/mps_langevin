function dA = inverseFreeTimestep( AL, AR, z, dt, H, d, D, KWARGS )
    %%INVERSEFREETIMESTEP Take inverse free timestep
    % This function takes a left gauge AL, right gauge AR, generator of
    % the environment z, AL's nullspace VL and a Hamiltonian H and updates AL
    % from AL(t) to AL(t+dt). AR(t+dt) and z(t+dt) must then be determined from
    % this
    
    % First, verify gauge conditions
    verify_gauge_conditions(AL, z, D, KWARGS.CFORM_TOL)
    
    if KWARGS.RK4
        %dA = inverseFreeTimestep_RK4( AL , AR , H , z, D , d , dt );
        error('Inverse-free RK4 not working properly yet')
    else
        % 1. Find AC = AL.z = z.AR
        AC = ncon({AL, z}, {[-1 1 -3],[1 -2]}); 

        % 2. find d/dt(AL).z
        VL = leftGaugeTangentBasis(eye(D), AL, d, D, KWARGS.VL_TOL);
        ALdotz =find_Adotz(AL,AR,VL,H,z,D);        

        % 3. find AL.d/dt(z)
        zdot = find_zdot(AL,ALdotz,AR);
        ALzdot = ncon({AL, zdot}, {[-1 1 -3],[1 -2]}); 

        % 4. find new AC, z
        if KWARGS.IMAG_TIME
            AC_new = AC - dt*(ALdotz + ALzdot);
            z_new = z - dt*zdot;
        else
            AC_new = AC - 1i*dt*(ALdotz + ALzdot);
            z_new = z - 1i*dt*zdot;
        end

        % 5. decompose into new AL and z
        AL_new = minimize_ALz_minus_AC(z_new, AC_new,D,d);      
        dA = AL_new - AL;
    end
end


%% Subroutines


function verify_gauge_conditions(AL, z, D, TOL)
    %% Check gauge conditions are satisfied
    
    I = diag(ones([D,1])); 
    R = z*z';

    % conditions are: IE = I, ER = R
    % 1. IE = I
    IE = ncon({AL, conj(AL)}, {[1 -1 2], [1 -2 2]});
    IE = round(IE,TOL);

    IEisI = isequal(IE, I);

    % 2. ER = R
    ER = ncon({AL, conj(AL), R}, {[-1 1 2], [-2 3 2], [1 3]});
    ER = round(ER,TOL);
    R = round(R, TOL);
    ERisR = isequal(ER, R);

    % Raise errors if necessary
    if ~IEisI && ~ERisR
        absdiff1 = abs(IE - I); absdiff2 = abs(ER - R);
        size_err = max([absdiff1(:), absdiff2(:)]);
        err.message = ['Neither gauge condition for (AL, z) is ',...
                       'satisfied. Largest elementwise abs error = ',...
                       num2str(size_err)];
        err.idenfitifier = 'infiniteChain:TDVP:inverseFreeTimestep:NEITHER';
        error(err);
    elseif ~ERisR
            absdiff = abs(ER - R);
            size_err = max(absdiff(:));
            err.message = ['Right gauge condition for (AL, z) is not ',...
                           'satisfied. Largest elementwise abs error = ',...
                           num2str(size_err)];
            err.idenfitifier = 'infiniteChain:TDVP:inverseFreeTimestep:RIGHT';
            error(err);
    elseif ~IEisI
            absdiff = abs(IE - I);
            size_err = max(absdiff(:));
            err.message = ['Left gauge condition for (AL, z) is not ',...
                           'satisfied. Largest elementwise abs error = ',...
                           num2str(size_err)];
            err.idenfitifier = 'infiniteChain:TDVP:inverseFreeTimestep:LEFT';
            error(err);
    end
end


function [ AL_new ] = inverseFreeTimestep_RK4(AL,AR,H,z,D,d,dt)
    %% 4th order Runge Kutta - inverse free algorithm.
    % 
    % I want to find the 4th order RK vesion of A(t+dt) 
    % To do this, we can't naively find 
    % A(t+dt) = A(t) + dt/6 * (A_k1 + 2*A_k2 + 2*A_k3 + A_k4)
    % This isn't properly normalized for example. 
    % Instead we need to find d(z)/dt and d(AC)/dt and construct A(t+dt) from
    % that. 

    % Basically I want to find 
    % z_k1 = f(t , A(t) )
    % z_k2 = f(t + dt/2 , A_k1(t + dt/2) )
    % z_k3 = f(t + dt/2 , A_k2(t + dt/2) )
    % z_k4 = f(t + dt , A_k3(t + dt) )

    % z(t+dt) = z(t) + dt/6 * (z_k1 + 2*z_k2 + 2*z_k3 + z_k4)

    % Basically I want to find 
    % AC_k1 = g(t , A(t))
    % AC_k2 = g(t + dt/2 , A_k1(t + dt/2) )
    % AC_k3 = g(t + dt/2 , A_k2(t + dt/2) )
    % AC_k4 = g(t + dt , A_k3(t + dt) )

    % AC(t+dt) = AC(t) + dt/6 * (AC_k1 + 2*AC_k2 + 2*AC_k3 + AC_k4)

    % From this I will use minimize_ALz_minus_AC to find the 
    % 4th order RK version of A(t+dt).

    %% 1. Find z_k1 and AC_k1 and also the AL,AR,z which will be used
    % as input to find k2 components.

    VL = leftGaugeTangentBasis(eye(D), AL, d, D);

    ALdotz = find_Adotz(AL,AR,VL,H,z,D);

    zdot = find_zdot(AL,ALdotz,AR);
    ALzdot = ncon({AL, zdot}, {[-1 1 -3],[1 -2]});

    AC = ncon({AL, z}, {[-1 1 -3],[1 -2]});

    % construct ztilde and ACtilde, these are the ingredients
    % needed to determine the new AL

    % In this case we want a dt/2 timestep. 
    ztilde = z - 1i*(dt/2)*zdot;
    ACtilde = AC - 1i*(dt/2)*(ALdotz + ALzdot);

    % perform minimization to find the new AL
    AL_1 = minimize_ALz_minus_AC(ztilde, ACtilde,D,d);

    % z_k1 and AC_k1 will eventually be used in the RK calculation.
    z_k1= -1i*zdot;
    AC_k1= -1i*(ALdotz + ALzdot);

    [AR_1,z_1]=find_AR_z_from_AL(AL_1,d,D);
       
    %% 2. Find z_k2 and AC_k2 and also the AL,AR,z which will be used
    % as input to find k3 components.
        
    VL_1 = leftGaugeTangentBasis(eye(D), AL_1, d, D);

    ALdotz = find_Adotz(AL_1,AR_1,VL_1,H,z_1,D);

    zdot = find_zdot(AL_1,ALdotz,AR_1);
    ALzdot = ncon({AL_1, zdot}, {[-1 1 -3],[1 -2]});


    % construct ztilde and ACtilde, these are the ingredients
    % needed to determine the new AL

    % In this case we want a dt/2 timestep. 
    ztilde = z - 1i*(dt/2)*zdot;
    ACtilde = AC - 1i*(dt/2)*(ALdotz + ALzdot);

    % perform minimization to find the new AL
    AL_2 = minimize_ALz_minus_AC(ztilde, ACtilde,D,d);

    z_k2= -1i*zdot;
    AC_k2= -1i*(ALdotz + ALzdot);

    [AR_2,z_2]=find_AR_z_from_AL(AL_2,d,D);

    %% 3. Find z_k3 and AC_k3 and also the AL,AR,z which will be used
    % as input to find k4 components.
        
    VL_2 = leftGaugeTangentBasis(eye(D), AL_2, d, D);

    ALdotz = find_Adotz(AL_2,AR_2,VL_2,H,z_2,D);

    zdot = find_zdot(AL_2,ALdotz,AR_2);
    ALzdot = ncon({AL_2, zdot}, {[-1 1 -3],[1 -2]});

    % construct ztilde and ACtilde, these are the ingredients
    % needed to determine the new AL

    % In this case we want a dt timestep. 
    ztilde = z - 1i*(dt)*zdot;
    ACtilde = AC - 1i*(dt)*(ALdotz + ALzdot);

    % perform minimization to find the new AL
    AL_3 = minimize_ALz_minus_AC(ztilde, ACtilde,D,d);

    z_k3= -1i*zdot;
    AC_k3= -1i*(ALdotz + ALzdot);

    [AR_3,z_3]=find_AR_z_from_AL(AL_3,d,D);

    %% 4. Find z_k4 and AC_k4.
    VL_3 = leftGaugeTangentBasis(eye(D), AL_3, d, D);

    ALdotz = find_Adotz(AL_3,AR_3,VL_3,H,z_3,D);

    zdot = find_zdot(AL_3,ALdotz,AR_3);
    ALzdot = ncon({AL_3, zdot}, {[-1 1 -3],[1 -2]});

    z_k4= -1i*zdot;
    AC_k4= -1i*(ALdotz + ALzdot);

    %% 5. We can now construct our final z and AC from which we can 
    % construct the final A(t+dt).
    z_final = z + (dt/6) * ( z_k1 + 2*z_k2 + 2*z_k3 + z_k4);
    AC_final = AC + (dt/6) * ( AC_k1 + 2*AC_k2 + 2*AC_k3 + AC_k4);

    AL_new = minimize_ALz_minus_AC(z_final, AC_final,D,d);
end


function zdot = find_zdot(AL, ALdotz, AR)
    %% find_zdot Calculates d/dt(z) using algorithm from Vid Stojevic
    % M.zdot = v
    %%
    %  =M=  =       ___
    %         (i)--|   |--(j)
    %              | M | 
    %        (ip)--|___|--(jp)
    %  
    
    D = size(AL,1);     
       
    % 1. Find M
    I = eye(D);
    M_term1 = ncon({I, I}, {[-2, -4], [-1, -3]});
    M_term2 = ncon({AL, conj(AR)}, {[-1, -3, 1],[-2, -4, 1]});
    
    M = M_term1 - M_term2; % M[i,j,ip,jp]
    M = reshape(permute(M, [3 1 4 2]),[D^2, D^2]); % M[(ip,i), (jp,j)]
    
    % 2. Find v
    
    v = ncon({ALdotz, conj(AR)}, {[-1 1 2],[-2 1 2]});  % v[jp,j]
    v = reshape(v.', [D^2,1]); % v[(jp,j)]
    
    % 3. Invert M to find zdot. 
    % We are assuming using pinv is OK.
    % If there are issues, this could be the cause.
   invM = pinv(M);
    
    zdot = invM*v; % zdot[(jp,j)]
   %Azdot=linsolve(M,v);
   
    zdot = reshape(zdot, [D, D]).'; % zdot[jp,j]
end  


function Adotz = find_Adotz(A_L,A_R,V_L,H,z,chi)
    %% Adotz Find d/dt(AL).z    
    
    % This finds the inverse of the transfer operator.
    E_L = ncon({A_L,conj(A_L)},{[-2 -4 1],[-1 -3 1]}); 
    E_L=reshape(E_L,chi^2,chi^2);
    left=eye(chi);
    right=z*z';
    
    % We construct a matrix Q which is a projector on to the 
    % subspace of eigenvalues < 1 for the transfer operator
    Q=eye(chi^2)-reshape(right.',[chi^2,1])*reshape(left.',[1,chi^2]);
    
    % Because E_L has its eigenvalue = 1 part projected out this should
    % be well behaved. (The eigenvalue = 1 is zero due to the constuction
    % of VL)
    
    E_N=reshape(Q*inv(eye(chi^2)-Q*E_L*Q)*Q,chi,chi,chi,chi);
    
    
    ncon({A_L,A_L,H,conj(A_L),conj(A_L)},...
            {[1 4 2],[4 -1 5],[3 7 2 5],[1 6 3],[6 -2 7]});
    x1=ncon({A_L,A_L,H,conj(A_L),conj(A_L),E_N,A_L,conj(V_L),z},...
            {[1 4 2],[4 8 5],[3 7 2 5],[1 6 3],[6 9 7],...
             [9 8 11 10],[10 13 12],[11 -1 12],[13 -2]});

    x2=ncon({A_L,z,A_R,H,conj(V_L), conj(A_R)},{[1 4 2],[4 5], [5 8 6],[3 7 2 6],...
                                                     [1 -1 3],[-2 8 7]});

    x3=ncon({A_L,z,A_R,H,conj(A_L),conj(V_L)},{[1 4 2], [4 5],[5 -2 6],[3 8 2 6],...
                                                     [1 7 3],[7 -1 8]});

    x=x1+x2+x3;
    Adotz=ncon({V_L,x},{[-1 1 -3],[1 -2]});
    
    % Perhaps the relevant factor of -1i should be introduced here.
end


function A_L=minimize_ALz_minus_AC(ztilde,ACtilde,chi,d)
    % Take A_C (i , j, sigma) and combine in to
    % ac ( (sigma,i), j)
    ac=reshape(permute(ACtilde,[3,1,2]),[d*chi,chi]);

    % qr on ztilde and ac;
    [q,r]=qr(ztilde);
    [Q,R]=qr(ac);

    % You want R and r to be positive along the diagonal
    % This will flip the signs of q/Q and r/R when necessary
    fix1=diag(sign(diag(r)));
    fix2=zeros(chi,1);
    
    % This construction works for chi=1, others don't.
    for n=1:chi
        fix2(n)=sign(sign(R(n,n))+0.1);
    end

    fix2 = diag([fix2;1*ones([(d-1)*chi,1])]);
    q=q*fix1;
    r=fix1*r;
    Q=Q*fix2;
    R=fix2*R;

    % Take the first chi components of Q((sigma,i),(delta,j))
    % Reshape to (i,j,sigma)
    A_L_q=permute(reshape(Q(:,1:chi),[d,chi,chi]),[2,3,1]);
    if chi==1
        A_L=A_L_q;
    else
        % Contract q' on to finally find A_L.
        A_L=ncon({A_L_q,q'},{[-1 1 -3],[1 -2]});
    end
end