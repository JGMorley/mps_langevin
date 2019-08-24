% Test the FC TDVP code by comparing our results for rate function:
% l(t) = -ln(|<psi_t|psi_0>|)/N.
% against analytical results for a finite-chain XXZ model of spin-1/2s.
%
% H = -J\sum_i (Sz_i Sz_(i+1) + g Sx_i)
%
% We quench from a Neel state (ground state if J=h=0, h_st!=0) to a
% Hamiltonian for t>0 with h=h_st = 0.
%
%-%-%-%
% NB - Analytical solution only holds for N a multiple of 4
%-%-%-%

N = 4;
spinDimList = 2*ones([1 N]);

% Initial state is |+x>^N
psi0 = ones(spinDimList) / sqrt(prod(spinDimList));
[mps_LCF, mps_RCF] = MPSdecomposition(psi0,2); 

%% XXZ Hamiltonian

J = 1; g = 0;
[Sx,Sy,Sz] = spinMatrices(1/2);
X = 2*Sx; Y = 2*Sy; Z = 2*Sz;
I = eye(2);

% 1. Create length N cell array of {I,I,I,...,I}
all_Is = cell([1 N]);
for k=1:N
    all_Is{k} = I;
end

% 2. Add XX terms starting at site n and Z term on n for n=1:N-1
H = 0.;
for n=1:N-1
    ZZn = all_Is;
    ZZn{n} = Z; 
    ZZn{n+1} = Z;
    H = H + J*kronlist(ZZn);
    
    Xn = all_Is;
    Xn{n} = X;
    H = H + g*kronlist(Xn);
end

% also need 1st-nth site coupling!
ZZn = all_Is;
ZZn{N} = Z; 
ZZn{1} = Z;
H = H + J*kronlist(ZZn);

% 3. now add Z term for Nth site
Xn = all_Is;
Xn{N} = X;
H = H + g*kronlist(Xn);

% 4. reshape and permute
H = permute(reshape(H, [fliplr(spinDimList) fliplr(spinDimList)]), 2*N:-1:1);

%%
T = 6; dt = 0.1; number_samples = floor(T/dt);
output = tdvpIntegratorFC( mps_LCF, H, T, dt, number_samples, 'LCF', true );
%%

nSamples = floor(T/dt);%length(output);
times = zeros([1 nSamples]);
rate_fn = zeros([1 nSamples]);
absoverlaps = zeros([1 nSamples]);
analytical_rate_fn = zeros([1 nSamples]);

for k=1:nSamples
    t = output{k,2};
    times(k) = t;
    psik = findStateTensor(output{k,1});
    f = abs(ncon({psi0,psik},{1:N,1:N}));
    absoverlaps(k) = f;
    rate_fn(k) = -log(f^2)/N;
    
    analytical_rate_fn(k) = -(2/N)*log(cos(J*t)^N + sin(J*t)^N);
end

%%

figure

subplot(1,2,1)
hold on
scatter(times,rate_fn,'.')
scatter(times,absoverlaps,'.')
scatter(times,analytical_rate_fn,'.')
legend('-ln|<psi0|psit>| / N','|<psi0|psit>|','analytical rate fn')

subplot(1,2,2)
hold on
scatter(times,rate_fn - analytical_rate_fn,'.')
legend('numerical analytical difference')

if false
    plot_Schmidts(output)
end