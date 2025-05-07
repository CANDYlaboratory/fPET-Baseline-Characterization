% fPETsimulateTAC - Simulates a TAC using the two-compartment model of tracer
% dynamics.
%
% Parameters:
%   kin_consts_in - 1x4 vector of kinetic constants [K1, k2, k3, k4]
%   AIF           - Arterial Input Function (array representing input concentrations)
%   n_timesteps   - Number of time steps in simulation (scalar)
%   T_step        - Time step for simulation (minutes)
%   paradigm      - Array that defines the task paradigm; percent change in k3
%   varargin      - Optional parameters:
%                   - 'activation': Flag indicating a task paradigm; requires
%                     a vector of percent changes to k3 as the following argument
%
% Returns:
%   TAC      

function TAC = fPETsimulateTAC(kin_consts_in, AIF, n_timesteps, T_step, varargin)

    paradigm = zeros(n_timesteps, 1);

    for i = 1:length(varargin)
        if strcmpi(varargin{i}, 'activation')
            paradigm = varargin{i+1};
            i = i+1;
        end
    end

    TAC = zeros(length(paradigm), 1);

    K1 = kin_consts_in(1)*ones(length(paradigm), 1);
    k2 = kin_consts_in(2)*ones(length(paradigm), 1);
    k3 = kin_consts_in(3)*ones(length(paradigm), 1) + kin_consts_in(3)*paradigm;
    k4 = kin_consts_in(4)*ones(length(paradigm), 1);

    TAC = simulateWithImplicitEuler([K1 k2 k3 k4], AIF, n_timesteps, T_step);
end

%% Implicit Euler Implementation
% % written by Grant Hartung

% Numerically integrates the tissue compartment model using the implicit
% Euler method.

function Csig = simulateWithImplicitEuler(coeffs, inlet_func, T, dt)

    k = coeffs;
    if size(k, 1) > 1
        for i = 1:size(k, 1)
            AMx(:, :, i) = getAMx(k(i, :));
        end
    else
        AMx = getAMx(k);
    end
    A0 = [0 0 0]';
    
    C1BC = inlet_func;   
    
    C = integrateIntime(AMx,T,dt,A0,C1BC,false);
    
    Csig = (C(2,:) + C(3,:))';

end

function C = integrateIntime(AMx,T,dt,A0,C1BC,ploton)

    nTimesteps = T;
    identityMx = eye(3);
    identityMx(1,:) = 0;  %zero out first row, bc it is a BC anyway
    
    firstRow = AMx(1,:,1);
    x = A0;
    xDynamic = zeros(length(x),nTimesteps);
    
    % % implicit euler is (eye-dt*A)*xNew = xOld
    if size(AMx, 3) == 1
        dynamicAMx = identityMx - dt.*AMx;
    end
    dynamicAMx(1,:) = firstRow;
    for i = 1:nTimesteps
        if size(AMx, 3) > 1
            dynamicAMx = identityMx - dt.*AMx(:,:,i);
            dynamicAMx(1,:) = firstRow;
        end
        x(1) = C1BC(i);
        xNew = dynamicAMx\x;
        xDynamic(:,i) = xNew;
        x = xNew;
    end

    if ploton
        clf
        plot((1:T)*dt,xDynamic')
        legend('c_a_r_t','c_t_i_s_s_F_r_e_e','c_t_i_s_s_B_o_u_n_d')
        c4 = xDynamic(2,:) + xDynamic(3,:);
        hold on, plot((1:T)*dt,c4);  
        legend('c_a_r_t','c_t_i_s_s_F_r_e_e','c_t_i_s_s_B_o_u_n_d','c_s_i_g_n_a_l')
    end
    C = xDynamic;

end

function A = getAMx(k)
    % algebraic-differential system
    % Row1 is algebraic BC for input
    % Row 2 and 3 are dynamic
    A = [1          0       0;
        k(1)    -k(2)-k(3)  0;
        0           k(3)    0]; 
end
