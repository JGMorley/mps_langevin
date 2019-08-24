%% Look for non-analyticities in Loschmidt echo

D = 8; % Bond dimension
J = 1.; g0 = 1.5; g1 = 0.2; % Initial, final quenching parameters

%% 1. Hamiltonian tensor
d = 2; % local spin dimension
[X,Y,Z] = spinMatrices((d-1)/2);
X = 2*X; Y = 2*Y; Z = 2*Z;
I = eye(d);

% H_{i,i+1} = -J*Z.Z -(Jg/2)(X.I + I.X)]

H = @(J,g) permute(...
                      reshape(...
             -J*kron(Z,Z) - (J*g)*kron(I,X),[d d d d]...
                              ),[2 1 4 3]...
                   );
               
%% 2. Imag time integrator to find ground state with J0, g0
[~,Arand] = canonicalForm({1,rand([D D d]) + 1i*rand([D D d]),eye(D)});
T = 5; dt = T/500;

samples_GS = tdvpIntegrator(Arand, H(J,g0), T, dt, 'IMAG_TIME',true,'INVERSE_FREE',true,...
     'SCHMIDT_TH',1e-8,'CFORM_TOL',8);

plot_energy_density(samples_GS,H(J,g0)); % for some reason getting -ve flip?

A0 = samples_GS{end,1}; % Approximate ground state

%% 3. Real time integrator starting at ground state of old Hamiltonian

T = 5; dt = T/2000;
samples = tdvpIntegrator(A0,H(J,g1), T, dt,'SCHMIDT_TH',1e-8,'RK4',true);

%% 4. Plot Loschmidts
% define exact rate function
epsilon = @(k, J, g) 2.*J.*sqrt((g - cos(k)).^2 + sin(k).^2);
theta = @(k, g) 0.5 .* atan(sin(k) ./ (g - cos(k)));
phi = @(k, g0, g1) theta(k,g0) - theta(k,g1);
integrand = @(z, k, J, g0, g1) log(cos(phi(k,g0,g1)).^2 + ...
                                sin(phi(k,g0,g1)).^2 .* exp(-2.*z.*epsilon(k,J,g1)));
f = @(g0, g1, J, z) -1*integral(@(k)integrand(z,k,J,g0,g1)./(2*pi),0,pi); 
l = @(g0, g1, J, t) f(g0,g1,J,1i*t) + f(g0,g1,J,-1i*t);
ratefnExact = @(t) l(g0,g1,J,t);

absRateFnData = zeros([2,length(samples)]); % exact along 2nd axis

progressWindow = waitbar(0, 'Calculating Loschmidt amplitudes');
% extract rate function from data
for n = 1:length(samples)
    % Find A, R, L
    [An,t] = samples{n,:};
    [~,Rn,Ln] = normalizeMPS(An);

    % find rateFnData(n)
    En0 = ncon({An,conj(A0)},{[-2 -4 1],[-1 -3 1]});
    En0 = reshape(En0,[D^2,D^2]);
    y = eigs(En0, 1);
    absRateFnData(1,n) = - 2*log(abs(y));
    %yData(n) = abs(y);

    absRateFnData(2,n) = ratefnExact(t);
    waitbar(n/length(samples))
end
close(progressWindow)
%%
figure, hold on
times = cell2mat(samples(:,2)).';
scatter(times,absRateFnData(1,:),'.')
scatter(times,absRateFnData(2,:),'.')
xlabel('t')
ylabel('|y(t)|')
legend('TDVP output', 'Exact')

%% Plotting functions
function plot_energy_density(samples,H)
    energy_density = zeros([1,length(samples)]);
    progressWindow = waitbar(0, 'Calculating energies');
    cform_warning = true;
    for k = 1:length(samples)
        A = samples{k,1};
        [A,r,l] = normalizeMPS(A);
        if strcmp(verifyCanonicalForm({1,A,r},8),'FALSE') && cform_warning
            warning('not canonical form!')
            cform_warning = false;
        end
        warning('off','ncon:suboptimalsequence');
        energy_density(k) = ncon({A,A,H,conj(A),conj(A),r,l},...
                           {[1 3 2],[3 5 4],[6 8 2 4],[10 7 6],[7 9 8],[5 9],[1 10]});
        warning('on','ncon:suboptimalsequence');
        
        waitbar(k/length(samples), progressWindow, 'Calculating energies');
    end
    close(progressWindow)
    
    figure
    times = cell2mat(samples(:,2)).';
    scatter(times,real(energy_density),'.');
    xlabel('t');
    ylabel('<H>');
end

function plot_Loschmidt_amplitude(samples,J,g0,g,A0,D)
    % define exact rate function
    epsilon = @(k, J, g) 2.*J.*sqrt((g - cos(k)).^2 + sin(k).^2);
    theta = @(k, g) 0.5 .* atan(sin(k) ./ (g - cos(k)));
    phi = @(k, g0, g1) theta(k,g0) - theta(k,g1);
    integrand = @(z, k, J, g0, g1) log(cos(phi(k,g0,g1)).^2 + ...
                                    sin(phi(k,g0,g1)).^2 .* exp(-2.*z.*epsilon(k,J,g1)));
    f = @(g0, g1, J, z) -1*integral(@(k)integrand(z,k,J,g0,g1)./(2*pi),0,pi); 
    l = @(g0, g1, J, t) f(g0,g1,J,1i*t) + f(g0,g1,J,-1i*t);
    ratefnExact = @(t) l(g0,g,J,t);
    
    absRateFnData = zeros([2,length(samples)]); % exact along 2nd axis
    
    progressWindow = waitbar(0, 'Calculating Loschmidt amplitudes');
    % extract rate function from data
    for n = 1:length(samples)
        % Find A, R, L
        [An,t] = samples{n,:};
        [~,Rn,Ln] = normalizeMPS(An);

        % find rateFnData(n)
        En0 = ncon({An,conj(A0)},{[-2 -4 1],[-1 -3 1]});
        En0 = reshape(En0,[D^2,D^2]);
        y = eigs(En0, 1);
        absRateFnData(1,n) = - 2*log(abs(y));
        yData(n) = abs(y);

        absRateFnData(2,n) = ratefnExact(t);
        waitbar(n/length(samples))
    end
    close(progressWindow)
    
    figure, hold on
    times = cell2mat(samples(:,2)).';
    scatter(times,absRateFnData(1,:),'.')
    scatter(times,absRateFnData(2,:),'.')
    xlabel('t')
    ylabel('|y(t)|')
    legend('TDVP output', 'Exact')
end