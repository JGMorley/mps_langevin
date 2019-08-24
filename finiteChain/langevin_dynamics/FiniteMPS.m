classdef finiteMPS < matlab.mixin.Copyable % allows copying
    properties
        tensors
        nSites
        localDimList
        gauge
        norm
    end
    methods
        function fMPS = finiteMPS(localDimList) % constructor function
            fMPS.localDimList = localDimList;
            fMPS.nSites = length(localDimList);
        end 
        
        %~% state manipulations
        function generate_rand(fMPS)
            % generate tensors with elements given by rand()
        end
        
        function generate_maxCatState(fMPS)
            % generate CAT state with max. possible Schmidt values
        end
        
        function normalize(fMPS)
            % divide first tensor by state norm
        end
        
        function check_valid_MPS(fMPS)
            % check bond dimensions are compatible
        end
        
        function compute_env(fMPS,n,RorL)
            % compute R, L or both environments
        end
        
        function variational_compression(fMPS,Dmax)
            % variationally compress bond dimension to Dmax
        end
        
        function get_gauge(fMPS)
            % work out gauge, LCF or RCF or MCFn=...
        end
        
        function canonicalize(fMPS,GAUGE,n)
            % put in right, left or mixed canonical form
        end
        
        %~% evolution
        function evolve_closed(H,EnvParams,tFinal,nSteps,imag_time)
            % evolve, set properties time_series and mps_series
        end
        
        function evolve_langevin(H,EnvParams,tFinal,nSteps,imag_time)
            % evolve, set properties time_series and mps_series
        end
    end
end