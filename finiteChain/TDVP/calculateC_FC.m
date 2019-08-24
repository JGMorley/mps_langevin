function [C,tensorList,legLinksList] = calculateC_FC(n, mpsIn, H, N, noC)
    %% Perform contractions to find C-tensor used in TDVP integration
    % NB Bottom right corner, ie. right of empty bit, needs fixing!
    %% Inputs
    % 
    % n: integer
    %   size index to be left open, 1 <= n <= N
    % mpsIn: MPS cell-array
    %   input state
    % H: tensor of rank 2N
    %   Hamiltonian tensor
    % N: integer 
    %   Number of sites in system
    %
    %% Outputs
    % C: tensor
    %   If n!=1 and n!=N, then C = C[i,j,sigma]
    %   If n==1 then C = C[sigma,i] (spin index first)
    %   If n==N then C = C[j,sigma] (bond index first)
    %   This follows the convention for the 1st and Nth MPS tensors.
    %
    %%
    if N==1
        A = mpsIn{2}{1};
        C = (H*A).';
        tensorList   = 'NoneApplicable';
        legLinksList = 'NoneApplicable';
        return
    end
    
    % Create cell-array of tensors {A1,...,AN,H,conj(A1),...,conj(AN)} 
    % (with conj(An) missing)
    tensorList = cell([1 2*N+1]);
    for k=1:N
        Ak = mpsIn{k+1}{1};
        tensorList{k} = Ak;
        tensorList{N+1+k} = conj(Ak);
    end
    tensorList{N+1} = H;
    tensorList(N+1+n) = [];
    
    % Create cell-array of leg links
    legLinksList = cell([1 2*N+1]);
    
    % 1st, Nth, N+2th and 2Nth entries (corners of diagram)
    legLinksList{1}     = [1 2];
    legLinksList{N}     = [(2*N-2) (2*N-1)];
    legLinksList{N+2}   = [2*N (2*N+1)];
    legLinksList{2*N+1} = [(4*N-6) (4*N-5)] + double(n==1)*[1 1];
    
    % Nth entry (for H)
    topRow = 2*(1:N) - ones([1 N]);
    bottomRow = topRow + (2*N-(n~=1))*ones([1 N]);
    bottomRow = [bottomRow(1:n-1) -3 bottomRow(n:N-1)];
    bottomRow(n+1:N) = bottomRow(n+1:N) - ones([1,N-n]);
    legLinksList{N+1} = [bottomRow topRow];
    
    for k=2:N-1
        % kth entry
        legLinksList{k} = 2*k*[1 1 1] + [-2 0 -1];
        % N+1+kth entry
        legLinksList{N+1+k} = 2*(N+k-1)*[1 1 1] + [-1 1 0];
    end
    
    % Overwrite bottom MPS indices to the right of the nth site
    for k=n+1:N-1
        shift = 3 - double(n==1);
        legLinksList{N+1+k} = legLinksList{N+1+k} - shift*ones([1 3]);
    end
    
    % Pop the nth MPS tensor on the bottom
    legLinksList(N+1+n) = [];
    
    % Overwrite with open bond indices to left and right of nth site
    
    if n==1
        legLinksList{N+1}(1) = -1;
        legLinksList{N+2}(1) = -2;
    elseif n==N
        legLinksList{N+1+n-1}(2) = -1;
        legLinksList{N+1}(N) = -2;
    else
        legLinksList{N+1+n-1}(2) = -1;
        legLinksList{N+1+n}(1)   = -2;
    end
    
    % Call ncon() to make C
    C = 0.;
    if noC ~= 'True'
        C = ncon(tensorList,legLinksList);
    end
end

%% As a diagram:
%
%  C[i,j,\sigma,n] = C[(*),n] = 
%
%           O-- ... --O--O--O-- ... --O
%           |_________|__|__|_________|
%          |___________________________|
%           |         |  |  |         |
%           O-- ... --O-(*)-O-- ... --O
%
% where n is the site index, i, j are the left, right open 
% bond indices and \sigma is the open physical index.
%
% 'O' = site-specific MPS tensor
% large rectangle = Hamiltonian tensor
%%
