function [ Sx, Sy, Sz, Ssq ] = spinMatrices( S )
    %spinMatrices Returns spin matrices in Zeeman basis with spin q number S
    %
    %Expressions taken from http://easyspin.org/documentation/spinoperators.html
    
    % Round S to nearest integer
    twiceS = double(int64(2*S));
    S = twiceS / 2;
    dim = twiceS + 1;
    SSp1 = S*(S+1);
    mValues = S:-1:-S;

    % Populate matrices element-wise, iterate over the diagonal
    
    Sx = zeros(dim);
    Sy = zeros(dim);
    Sz = zeros(dim);
    Ssq = zeros(dim);
    
    % 1st row
    Sx(1,2) = sqrt( SSp1 - S*(S-1) )/2;
    Sy(1,2) = sqrt( SSp1 - S*(S-1) )/(2*1i);
    Sz(1,1) = S;
    Ssq(1,1) = SSp1;
    
    % middle rows
    for iRow=2:dim-1
        mm = mValues(iRow-1);
        m0 = mValues(iRow);
        mp = mValues(iRow+1);
        
        % Sx
        Sx(iRow, iRow+1) = sqrt( SSp1 - m0*mp ) / 2;
        Sx(iRow, iRow-1) = sqrt( SSp1 - m0*mm ) / 2;
        
        % Sy
        Sy(iRow, iRow+1) = sqrt( SSp1 - m0*mp ) / ( 2*1i );
        Sy(iRow, iRow-1) = -sqrt( SSp1 - m0*mm ) / ( 2*1i );
        
        % Sz, Ssq
        Sz(iRow, iRow) = m0;
        Ssq(iRow, iRow) = SSp1;
    end
    
    % last row
    Sx(end, end-1) = sqrt( SSp1 - (-S)*(-S+1) ) / 2;
    Sy(end, end-1) = - sqrt( SSp1 - (-S)*(-S+1) ) / ( 2*1i );
    Sz(end) = -S;
    Ssq(end) = SSp1;    
end