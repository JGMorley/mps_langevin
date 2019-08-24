function [A_L_new,A_R_new,z_new,chi_new] = increase_bond_dim(AL,AR,z,H,chi,D,d)
    % This algorithm finds a new AL and AR and z 
    % which will be used in the time evolution algorithm
    % to find AL(t+dt) etc. 

    % chi_new is the new bond dimension, it will be either
    % chi (if Y=zeros), d*chi_old (in most cases) 
    % or D (if the necessary bond dimension has been reached).

    if ischar(AR) && ischar(z)
        if strcmpi(AR,'none') && strcmpi(z,'none')
           if chi==1
            AR=AL;
            z=1;
           else

            [AR,z]=find_AR_z_from_AL(AL,z,d,chi);
           end
        end
    end
    
    
    if chi<D
        % VL and VR form the basis from which the 
        % update will be determined and Y are the 
        % coefficients of these basis elements.
        VL = leftGaugeTangentBasis(eye(chi), AL, d, chi);
        VR=rightGaugeTangentBasis(z,AL,d,chi);

        Y=calculateY(AL,eye(chi),z,VL,VR,H);

        % Round Y so bond dimension doesnt increase
        % when it shouldn't.
        Y=round(Y,14);

        [u,s,v]=svd(Y);
        v=v';
        Z12=u*sqrt(s);
        Z21=sqrt(s)*v;
        B12=ncon({VL,Z12},{[-1 1 -3],[1 -2]});
        B21=ncon({Z21,VR,z},{[-1 1],[1 2 -3],[2 -2]});

        % Deals with Y being zero by not increasing
        % the bond dimension.
        if Y==zeros(size(VL,1),size(VR,1))
            A_L_new=AL;
            A_R_new=AR;
            z_new=z;
            chi_new=chi;

        else
            % Distinguish between chi*d being > D and < D.
            % If chi*d > D then we want to truncate the matrices
            % we are using so the bond dimension is correct.
            % The truncation is OK due to things being ordered by
            % size of singular value.
            if D-chi < (d-1)*chi
                delta_D=D-chi;
                A_L_new=zeros(D,D,d);
                A_R_new=zeros(D,D,d);
                z_new=zeros(D);

                % Creating the new AL with the old AL and B12
                % which pads it out. 
                A_L_new(1:chi,1:chi,:)=AL;
                A_L_new(1:chi,(chi+1):(chi+delta_D),:)=B12(1:chi,1:delta_D,:);

                % Creating the new AR with the old AL and B21
                % which pads it out.
                A_R_new(1:chi,1:chi,:)=AR;
                A_R_new((chi+1):(chi+delta_D),1:chi,:)=B21(1:delta_D,1:chi,:);

                % The new z is the old z with nothing extra added.
                z_new(1:chi,1:chi)=z;
                chi_new=D;

            else
                A_L_new=zeros(d*chi,d*chi,d);
                A_R_new=zeros(d*chi,d*chi,d);
                z_new=zeros(d*chi);

                % Creating the new AL with the old AL and B12
                % which pads it out.                   
                A_L_new(1:chi,1:chi,:)=AL;
                A_L_new(1:chi,(chi+1):(chi+(d-1)*chi),:)=B12;

                % Creating the new AR with the old AL and B21
                % which pads it out.

                A_R_new(1:chi,1:chi,:)=AR;
                A_R_new((chi+1):(chi+(d-1)*chi),1:chi,:)=B21;

                % The new z is the old z with nothing extra added.
                z_new(1:chi,1:chi)=z;
                chi_new=d*chi;
            end
        end
    else
        % If this error happens then fix the integrator so it doesn't!
        error('Dont want to increase the bond dimension, it is big enough!')
    end
end