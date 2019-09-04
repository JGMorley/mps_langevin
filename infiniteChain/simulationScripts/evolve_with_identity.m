%% Simple test script for evolving with identity

%% setup

d = 2; % physical dimension
D = 2; % bond dimension

A0 = rand([2 2 2]); % MPS with randomized elements
                    % NB not normalized etc
                    
H = eye(d^2);       % 2 site identity matrix
H = permute(reshape(H,[d d d d]),[2 1 4 3]);  % reshape into tensor

T  = 1;    % Final time of evolution
dt = 0.1;  % timestep

%% evolution

samples = tdvpIntegrator(A0, H, T, dt);

%% analysis

% various gauge transformations will be made, so may not directly see
% unchanging elements

nSamples = length(samples);
time_series = zeros([1 nSamples]);
fidelity_series = zeros([1 nSamples]);
for k=1:nSamples
    time_series(k) = samples{k,2};
    Ak = samples{k,1};
    fidelity_series(k) = find_fidelity(A0,Ak);
end

%% plot

figure
hold on
box on
plot(time_series,1 - fidelity_series,'*-')
xlabel('t')
ylabel('1 - F(A(0),A(t)')
title("MPS doesn't change under evolution with identity")
xlim([0 T])
ylim(max(10*eps,max(abs(1-fidelity_series)))*[-1 1])
set(gcf,'color','w')

%% function for finding fidelity

function abs_overlap = abso(A0, A1)
    % find |<A0|A1>| / N, via largest eigenvalue of transfer matrix
    D = size(A0,1);
    assert(isequal(size(A0),size(A1)),'MPS are different sizes')
    
    E = ncon({A0,conj(A1)},{[-2 -4 1],[-1 -3 1]});
    E = reshape(E,[D^2,D^2]);
    abs_overlap = abs(eigs(E, 1)); % absolute value of largest eigenvector
end


function F = find_fidelity(A0, A1)
    % find F = |<A0|A1>| /  sqrt( |<A0|A0> * <A1|A1> )
    
    F = abso(A0,A1) / sqrt(abso(A0,A0) * abso(A1,A1));
end