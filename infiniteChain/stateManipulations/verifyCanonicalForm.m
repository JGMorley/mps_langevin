function canform = verifyCanonicalForm(mpsIn, TOL)
    %Verifies that mpsIn is in canonical form up to a tolerance of TOL digits
    %  canform = verifyCanonicalForm(mpsIn, TOL) returns string 'TRUE' if
    %    the conditions of eqs. 4,5 of arXiv:0711.3960, or their equivalents, 
    %    are satisfied up to a tolerance of TOL digits, and the string 'FALSE' 
    %    if they aren't.
    
    if nargin < 2
        TOL = 10;
    end
    
    FORM = mpsIn{1};
    
    switch FORM
        case 1 % left canonical form
            canform = verifyCanonicalFormA(mpsIn, TOL);
        case 2 % right canonical form
            canform = verifyCanonicalFormA(mpsIn, TOL);
        case 3 % Vidal form
            canform = verifyCanonicalFormGL(mpsIn, TOL);
    end
end

function canform = verifyCanonicalFormGL(mpsIn, TOL)
    % verifies that an mpsIn in Vidal form is in canonical form  
    % Assume FORM==3
    gammaIn  = mpsIn{2};
    lambdaIn = mpsIn{3};
    
    if nargin < 2
        TOL = 10;
    end
    
    % Extract bond order D
    [D,~,~] = size(gammaIn);
    
    identity = diag(ones([D,1]));
    
    tensors = {gammaIn, lambdaIn, conj(gammaIn), conj(lambdaIn)};
    legLinksRI = {[-1 2 3], [2 4], [-2 5 3], [5 4]};
    legLinksIL = {[2 -1 4], [1 2], [3 -2 4], [1 3]};

    canform = 'FALSE';
 
    RI = ncon(tensors, legLinksRI);
    RI = RI/RI(1,1);  % normalize by 1st element to eleminate eta
    RI = round(RI,TOL);
    
    IL = ncon(tensors, legLinksIL);
    IL = IL/IL(1,1);  % normalize by 1st element to elminate eta
    IL = round(IL,TOL);
    
    RIisI = isequal(RI,identity);
    ILisI = isequal(IL,identity);
    
    if RIisI && ILisI
        canform = 'TRUE';
    end
end

function canform = verifyCanonicalFormA(mpsIn, TOL)
    % verifies that an mpsIn given in terms of A-tensors is in canonical form
    
    canform = 'FALSE';
    
    if nargin < 2
        TOL = 10; % This may need tuning lower than in canonicalForm() 
    end
    
    [FORM, A, L] = mpsIn{:};
    
    % Extract bond order D to calculate identity
    [D,~,~] = size(A);    
    I = diag(ones([D,1])); 
    
    switch FORM
        case 0
            error('function called with FORM = 0');
        case 1
            % conditions are: IE = I, EL = L
            % 1. IE = I
            IE = ncon({A, conj(A)}, {[1 -1 2], [1 -2 2]});
            IE = IE/IE(1,1); % normalize
            IE = round(IE,TOL);
            
            IEisI = isequal(IE, I);
            
            % 2. EL = L
            EL = ncon({A, conj(A), L}, {[-1 1 2], [-2 3 2], [1 3]});
            EL = EL*L(1,1)/EL(1,1);
            EL = round(EL,TOL);
            
            L = round(L, TOL);
            ELisL = isequal(EL, L);
            
            % Now check that both IEisI and ELisL
            if IEisI && ELisL
                canform = 'TRUE';
            end
        case 2
            % conditions are: LE = L, EI = I
            % 1. LE = L
            LE = ncon({L, A, conj(A)}, {[1 2], [1 -1 3], [2 -2 3]});
            LE = LE*(L(1,1)/LE(1,1)); % normalize with respect to L
            LE = round(LE,TOL);
            
            L = round(L, TOL);
            LEisL = isequal(LE, L);
            
            % 2. EI = I
            EI = ncon({A, conj(A)}, {[-1 1 2], [-2 1 2]});
            EI = EI/EI(1,1);
            EI = round(EI,TOL);
            
            EIisI = isequal(EI, I);
            
            % Now check that both IEisI and ELisL
            if LEisL && EIisI
                canform = 'TRUE';
            end
    end
end