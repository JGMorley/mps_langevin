classdef EnvironmentParams < matlab.mixin.Copyable % allows copying
    properties
        nSites
        CouplingOperators
        CouplingStrengths
        CouplingNoises
        CouplingTemperatures
        FixedT = []
        FixedCouplingStrength = []
    end
    methods
        function EP = EnvironmentParams(N) % constructor function
            EP.nSites = N;
            EP.CouplingOperators = num2cell(cell([1,N])); % (cell of cells)
            EP.CouplingStrengths = num2cell(cell([1,N]));
            EP.CouplingNoises = num2cell(cell([1,N]));
            EP.CouplingTemperatures = num2cell(cell([1,N]));
        end
        
        function add_Fn(EP,n,F,CouplingStrength,T)
            % add a coupling operator to the nth site, with optionally
            % specified CouplingStrength and T(emperature).
            
            % Add F to nth site
            if ~isequal(EP.CouplingOperators{n},{[]})       
                EP.CouplingOperators{n}{end+1} = {};
            end
            EP.CouplingOperators{n}{end} = F;
            
            % Initialize empty cells for CouplingStrength and T
            if ~isequal(EP.CouplingStrengths{n},{[]}) 
                EP.CouplingStrengths{n}{end+1} = {};
            end
            if ~isequal(EP.CouplingTemperatures{n},{[]})
                EP.CouplingTemperatures{n}{end+1} = {};
            end
            
            if nargin==3
                if ~isequal(EP.FixedT,[])
                    EP.CouplingTemperatures{n} = EP.FixedT;
                end
                if ~isequal(EP.FixedCouplingStrength,[])
                    EP.CouplingStrength{n}{end} = EP.FixedCouplingStrength;
                end
            elseif (nargin==4)||(nargin==5)
                % Specify CouplingStrength
                EP.CouplingStrengths{n}{end} = CouplingStrength;
                EP.FixedCouplingStrength = [];
            end
            if nargin==5
                % Specify T(emperature)
                EP.CouplingTemperatures{n}{end} = T;
                EP.FixedT = [];
            end
            
            if nargin~=3 && nargin~=4 && nargin~=5
                error('Invalid number input arguments')
            end
        end
        
        function set_FixedT(EP,T)
            for n=1:EP.nSites
                for k=1:length(EP.CouplingTemperatures{n})
                    if ~isequal(EP.CouplingTemperatures{n}{k},[])
                        EP.CouplingTemperatures{n}{k} = T;
                    end
                end
            end
            EP.FixedT = T;
        end
        
        function set_FixedCouplingStrength(EP,CouplingStrength)
            for n=1:EP.nSites
                for k=1:length(EP.CouplingStrengths{n})
                    if ~isequal(EP.CouplingStrengths{n}{k},[])
                        EP.CouplingStrengths{n}{k} = CouplingStrength;
                    end
                end
            end
            EP.FixedCouplingStrength = CouplingStrength;
        end
        
        function generate_noises(EP,dt,NOISE_WOUT_DISSIPATION)
            if nargin<3
                NOISE_WOUT_DISSIPATION = false;
            end
            % generate stochastic noise terms for each coupling operator
            EP.check_valid(false)
            for n=1:EP.nSites
                EP.CouplingNoises{n} = cell([1,length(EP.CouplingOperators{n})]);
                for k=1:length(EP.CouplingOperators{n})
                    gamma = EP.CouplingStrengths{n}{k};
                    temp = EP.CouplingTemperatures{n}{k};
                    if NOISE_WOUT_DISSIPATION
                        noise_nk = randn()*sqrt(2*temp*dt);
                    else
                        noise_nk = randn()*sqrt(2*gamma*temp*dt);
                    end
                    EP.CouplingNoises{n}{k} = noise_nk;
                end
            end
        end
        
        function check_valid(EP, check_noises_too)
            % Check Operators, Noises, Strengths and Temperatures are of
            % compatible dimensions and are of square numeric type
            if nargin==1
                check_noises_too = true;
            end
            
            for n = 1:EP.nSites
                op_length = length(EP.CouplingOperators{n});
                
                if isequal(EP.CouplingOperators{n},{[]})
                    % no coupling operator specified, so no need to check
                    continue
                end
                
                % 1. Check dimensions
                if op_length ~= length(EP.CouplingStrengths{n})
                    errmsg = ...
                       ['CouplingStrengths is wrong length.',newline,...
                        'length(CouplingOperators) = ',num2str(op_length),...
                        newline,'length(CouplingStrengths) = ',...
                        num2str(length(EP.CouplingStrengths{n}))','.'];
                    error(errmsg)
                elseif op_length ~= length(EP.CouplingTemperatures{n})
                    errmsg = ...
                       ['CouplingTemperatures is wrong length.',newline,...
                        'length(CouplingOperators) = ',num2str(op_length),...
                        newline,'length(CouplingTemperatures) = ',...
                        num2str(length(EP.CouplingTemperatures{n}))','.'];
                    error(errmsg)
                elseif check_noises_too
                    if op_length ~= length(EP.CouplingNoises{n})
                        errmsg = ...
                          ['CouplingNoises is wrong length.',newline,...
                           'length(CouplingOpeNoises) = ',...
                           num2str(length(EP.CouplingNoises{n}))','.'];
                        error(errmsg)        
                    end
                end
                for k=1:op_length
                    % 2. Check type
                    Fnk = EP.CouplingOperators{n}{k};
                    Tnk = EP.CouplingTemperatures{n}{k};
                    Gnk = EP.CouplingStrengths{n}{k};
                    
                    try 
                        Fnk_transpose = Fnk.';
                    catch 
                        errmsg = [num2str(k),'th coupling operator on ',...
                                  num2str(n),'th site is not rank 2.'];
                        error(errmsg)
                    end
                    
                    if ~isequal(size(Fnk), size(Fnk_transpose))
                        errmsg = [num2str(k),'th coupling operator on ',...
                                 num2str(n),'th site is not square.',newline,...
                                 'Its size is ',num2str(size(Fnk)),'.'];
                        error(errmsg)
                    elseif ~isnumeric(Tnk) || ~isequal(size(Tnk),[1 1])
                        errmsg = [num2str(k),'th coupling temperature on ',...
                                  num2str(n),'th site is invalid..'];
                        error(errmsg)
                    elseif ~isnumeric(Gnk) || ~isequal(size(Gnk),[1 1])
                        errmsg = [num2str(k),'th coupling strength on ',...
                                  num2str(n),'th site is invalid.'];
                        error(errmsg)
                    elseif check_noises_too
                        Enk = EP.CouplingNoises{n}{k};
                        if ~isnumeric(Enk) || ~isequal(size(Enk),[1 1])
                            errmsg = [num2str(k),'th coupling noise on ',...
                                      num2str(n),'th site is invalid.'];
                            error(errmsg);
                        end
                    end
                end
            end
        end
        
        function set_manual_noise(EP,manual_noise,idx)
            % set CouplingNoises from idx elements of manual_noise, which
            % should be a cell-array of the same structure.
            
            % check manual_noise has right structure
            assert(iscell(manual_noise),'manual_noise not a cell')
            assert(length(manual_noise)==EP.nSites,'manual_noise wrong length')
            for n=1:EP.nSites
                nOps = length(EP.CouplingOperators{n});
                if isequal(EP.CouplingOperators{n},{[]}), nOps=0; end
                assert(length(manual_noise{n})==nOps,'manual_noise wrong shape')
            end
            
            % set CouplingNoises
            for n=1:EP.nSites
                nOps = length(EP.CouplingOperators{n});
                if isequal(EP.CouplingOperators{n},{[]}), nOps=0; end
                for k=1:nOps
                    EP.CouplingNoises{n}{k} = manual_noise{n}{k}(idx);
                end
            end
        end
        
        function answer = is_frictionless(EP)
            % determine whether EP is frictionless by checking each 
            % coupling strength
            for n=1:EP.nSites
                nOps = length(EP.CouplingOperators{n});
                for k=1:nOps
                    gamma = EP.CouplingStrengths{n}{k};
                    frictionless = isequal(gamma,0) || isequal(gamma,[]);
                    if ~frictionless
                        answer = false;
                        return
                    end
                end
            end
            
            answer = true;
        end
        
        function answer = is_empty(EP)
            % determine whether there are any coupling operators
            for n=1:EP.nSites
                nOps = length(EP.CouplingOperators{n});
                for k=1:nOps
                    if ~isempty(EP.CouplingOperators{n}{k})
                        answer = false;
                        return
                    end
                end
            end
            
            answer = true;
        end
                        
                 
    end
end