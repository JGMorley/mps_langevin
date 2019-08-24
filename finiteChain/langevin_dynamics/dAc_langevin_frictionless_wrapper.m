function [dAc,h,qs,debug_structs] = dAc_langevin_frictionless_wrapper( n, mps_mixed, H, dt, ...
                                                 EnvParams, scheme, ...
                                                 varargin )                                 
    %% Find update dAc using scheme
    
    % check valid scheme
    valid_scheme_vals = {'Euler','RK2','RK4'};
    scheme_match = zeros([1 length(valid_scheme_vals)]);
    for k=1:length(valid_scheme_vals)
        scheme_match(k) = isequal(scheme,valid_scheme_vals{k});
    end
    assert(any(scheme_match),'invalid value for scheme')

    % find first-order update
    [dA0, h_k1, qs, debug_k1] = dAc_langevin_frictionless(n,mps_mixed,H,dt,EnvParams,varargin);
    debug_structs.k1 = debug_k1;
    
    dAc_noise = dA0 - h_k1*dt;
    
    if isequal(scheme,'Euler')
        h = h_k1;
        dAc = h_k1*dt + dAc_noise;
        return
    end
%     mps_k1 = mps_mixed;
%     mps_k1{n+1}{1} = mps_k1{n+1}{1} + 0.5*h_k1*dt;
%     
%     % get k2
%     [~,h_k2,~,debug_k2] = dAc_langevin_Euler( n, mps_k1, H, dt, EnvParams, varargin );
%     debug_structs.k2 = debug_k2;
%     if isequal(scheme,'RK2')
%         h = h_k2;
%         dAc = h_k2*dt + dAc_noise;
%         return
%     end
% 
%     mps_k2 = mps_mixed;
%     mps_k2{n+1}{1} = mps_k2{n+1}{1} + 0.5*h_k2*dt;
%     
%     % get k3
%     [~,h_k3,~,debug_k3] = dAc_langevin_Euler( n, mps_k2, H, dt, EnvParams, varargin );
%     debug_structs.k3 = debug_k3;
% 
%     mps_k3 = mps_mixed;
%     mps_k3{n+1}{1} = mps_k3{n+1}{1} + h_k3*dt;
%     
%     % get k4
%     [~,h_k4,~,debug_k4] = dAc_langevin_Euler( n, mps_k3, H, dt, EnvParams, varargin );
%     debug_structs.k4 = debug_k4;
%     
%     % output
%     h = (1/6)*(h_k1 + 2*h_k2 + 2*h_k3 + h_k4);
%     dAc = h*dt + dAc_noise;
%     qs = NaN;
end