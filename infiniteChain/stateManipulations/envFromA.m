function [Lambda, E] = envFromA( A, FORM, fast )
    % Finds environment Lambda from A, FORM, given A is correctly normalized
    %
    % Checks that A is normalized unless fast is set to 1.
    % fast takes values 0 or 1.
    if nargin < 3
        fast = 0;
    else
        if (fast~=0 && fast~=1)
            err.identifier = 'envFromA:InvalidInput:fast';
            err.message = 'Invalid value of fast. Should be 0 or 1.';
            error(err);
        end
    end
    
    % check A is normalized if fast=0
    D = size(A,1);
    if fast == 0
        norm = ncon({A, conj(A)}, {[1 2 3], [1 2 3]}) / D;
        if round(norm,8)~=1
            err.identifier = 'envFromA:InvalidInput:A';
            err.message = 'Input tensor A is not normalized';
            error(err);
        end
    end
    
    % Find Lambda  
    E = ncon({A,conj(A)},{[-2 -4 1],[-1 -3 1]});
    E = reshape(E,[D^2,D^2]);
    
    if FORM==1
        % left canonical form
        [Lambda,~] = eigs(E,1);
    elseif FORM==2
        % right canonical form
        [Lambda,~] = eigs(E.',1);
    end
    
    % reshape as a matrix
    Lambda = reshape(Lambda,[D,D]);
    
    % normalize Lambda
    phase = sign(Lambda(1,1)); % phase factor of first element
    norm = trace(Lambda / phase);
    Lambda = Lambda / (norm*phase);
end

