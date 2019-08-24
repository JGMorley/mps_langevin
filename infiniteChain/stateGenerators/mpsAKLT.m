function state = mpsAKLT(FORM)
    %Function to return a MPS describing the ground state of the 1d AKLT model 
    %  state = mpsAKLT(FORM) gives
    %    for FORM = 3, AKLT state in Vidal form
    %    for FORM = 1, AKLT state in left-canonical form, and
    %    for FORM = 2, AKLT state in right-canonical form.
    %
    %  see Orus' 2013 review at arXiv:1306.2164
    
    % The matrix associated with each physical index value is a Pauli
    % matrix
    state{1} = FORM; % NB state{2} and state{3} are the same in any form
    state{2} = cat(3, [0 1;1 0], [0 -1i;1i 0], [1 0;0 -1]);
    state{3} = [1 0;0 1];
    
    %NB p = 3, D = 2
end