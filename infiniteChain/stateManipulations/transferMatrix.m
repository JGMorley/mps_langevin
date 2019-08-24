function [E, VL, VR, eValue] = transferMatrix(mpsIn)
%Find the transfer matrix and dominant left/right eigenvectors of mpsIn
%
%   E: 
%      the transfer matrix as a rank-4 tensor with indices labelled in order
%      (top-left,top-right,bottom-left,bottom-right) when drawn in the usual way
%   VL:
%      dominant left-eigenvector as a DxD matrix normalized by first
%      element phase
%   VR:
%      dominant right-eigenvector as a DxD matrix normalized by first
%      element phase
%   eValue:
%      dominant eigenvalue. By construction E will be Hermitian so this is
%      unambiguous
%
%   NB won't work with mpsIn in Vidal form
    
    if mpsIn{1} == 3
        err.message = 'function called with incompatible Vidal form MPS';
        err.identifier = 'transferMatrix:VidalFormIncompatible';
        error(err);        
    end
    
    A = mpsIn{2};
    D = size(A,1);
    
    % Construct transfer matrix E and reshape appropriately
    E = ncon({A,conj(A)},{[-1 -2 1],[-3 -4 1]});
    E = permute(E,[3,1,4,2]);
    E = reshape(E,[D^2,D^2]);
    
    % Finding the dominant eigenvectors of R and L, V_R and V_L respectively,
    % and the corresponding eigenvalues.
    [VR,Eta_R]=eigs(E,1);
    [VL,Eta_L]=eigs(E.',1);

    if abs(Eta_R - Eta_L) > 1e-10
        err.message = 'Dominant left and right eigenvalues unequal';
        err.identifier = 'transferMatrix:evaluesNotEqual';
        error(err);
    end
    
    eValue = Eta_R;

    % Reshape V_R and V_L from D^2 element vectors to DxD matrices
    VR=reshape(VR,[D,D]).';
    VL=reshape(VL,[D,D]);

    % normalize by dividing by normalized first element
    VR = VR/sign(VR(1,1)); % sign(z) = z./abs(z) for complex z
    VL = VL/sign(VL(1,1));
end

