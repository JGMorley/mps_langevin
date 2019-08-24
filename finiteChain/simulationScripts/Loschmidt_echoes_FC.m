% Test the FC TDVP code by comparing our results for rate function:
% l(t) = -ln(|<psi_t|psi_0>|)/N.
% against analytical results for a finite-chain XXZ model of spin-1/2s.
%
% H = J\sum_j[Sx_j Sx_(j+1) + Sy_j Sy_(j+1)] + \sum_n (h+h_dt) Sz_j
%
% We quench from a Neel state (ground state if J=h=0, h_st!=0) to a
% Hamiltonian for t>0 with h=h_st = 0.
%
%-%-%-%
% NB - Analytical solution only holds for N a multiple of 4
%-%-%-%

N = 8;
spinDimList = 2*ones([1 N]);

% Initial state is Neel state in z-direction
psi0 = zeros(spinDimList);
NeelComponent = cell([1 N]);
temp = 1;
for n=1:N
    if temp==1
        NeelComponent{n} = temp;
        temp = temp + 1;
    else
        NeelComponent{n} = temp;
        temp = temp - 1;
    end
end
psi0(NeelComponent{:}) = 1;
[mps_LCF, mps_RCF] = MPSdecomposition(psi0,2); 

%% XXZ Hamiltonian

J = 1; h = 0;
[Sx,Sy,Sz] = spinMatrices(1/2);
I = eye(2);

% 1. Create length N cell array of {I,I,I,...,I}
all_Is = {};
for k=1:N
    all_Is = {all_Is{:} , I};
end

% 2. Add XX terms starting at site n and Z term on n for n=1:N-1
H = 0.;
for n=1:N-1
    XXn = all_Is;
    XXn{n} = Sx; 
    XXn{n+1} = Sx;
    H = H + J*kronlist(XXn);
    
    YYn = all_Is;
    YYn{n} = Sy; 
    YYn{n+1} = Sy;
    H = H + J*kronlist(YYn);
    
    Zn = all_Is;
    Zn{n} = Sz;
    H = H + h*kronlist(Zn);
end

% 3. now add Z term for Nth site and Nth-1st site couplings
XXn = all_Is;
XXn{N} = Sx; 
XXn{1} = Sx;
H = H + J*kronlist(XXn);

YYn = all_Is;
YYn{N} = Sy; 
YYn{1} = Sy;
H = H+ J*kronlist(YYn);
    
Zn = all_Is;
Zn{N} = Sz;
H = H + h*kronlist(Zn);

% 4. reshape and permute
H = permute(reshape(H, [fliplr(spinDimList) fliplr(spinDimList)]), 2*N:-1:1);

%%
T = 1; dt = 0.01; number_samples = floor(T/dt);
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
    
    lt = 0;
    for n=1:floor(N/2)
        if mod(N, 2) == 0 % N even
            a = pi/N;
        else
            a = pi/N;
        end
        kn = 2*pi*n/N;
        ekn = -J*cos(kn + a);
        
        lt = lt - (1/N)*log(cos(t*ekn)^2);
    end        
    analytical_rate_fn(k) = lt;
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