d=2;
D=2;
% 1. create random mps
mps = mpsGenerator(d, D, 1);

% 2. put mps in LCF and get A
mps =  canonicalForm(mps);
A = mps{2};
right = mps{3};

% 3. run tangentBasis_unitaryGauge and check the outputs
[ A_new, right_new, U, Lambda ] = tangentBasis_unitaryGauge( A, right );

%%
% Check U(:,:,:,1) with 
U(:,:,:,1) - A_new

% Check tangent space properties with
ncon({U(:,:,:,2:end), conj(A_new)}, {[1 -1 2 -2], [1 -3 2]})

% Check U is unitary with (incorrect atm)
UU = permute(U,[1 3 2 4]);
UU = reshape(U, [d*D d*D]);
UU*UU' - eye(d*D)
UU'*UU - eye(d*D)