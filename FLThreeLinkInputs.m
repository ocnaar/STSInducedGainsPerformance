%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez-Aroche                                                  %
% Berkeley Center for Control and Identification                          %
% Spring 2017                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Retrieve input applied to a three-link robot by a controller based on   %
% feedback linearization, and control allocation from the closed-loop     %
% simulation of the system.                                               %
%                                                                         %
% Input                                                                   %
%                                                                         %
% x: n by 6 array with state trajectory from closed-loop simulation.      % 
% 	x(:,1): angular position of link 1 relative to the horizontal in [rad]%
% 	x(:,2): angular position of link 2 relative to link 1 in [rad].       % 
% 	x(:,3): angular position of link 3 relative to link 2 in [rad].       %
% 	x(:,4): angular velocity of link 1 in [rad/s].                        %
% 	x(:,5): angular velocity of link 2 in [rad/s].                        %
% 	x(:,6): angular velocity of link 3 in [rad/s].                        %
% xd: n by 6 array with reference state trajectory.                       % 
% par: structure with nominal parameters of the system, and controller    %
% settings.                                                               %
% 	par.lb: lower bounds for system input in control allocation.          % 
% 	par.ub: upper bounds for system input in control allocation.          %
% 	par.m1: mass of link 1 in [kg].                                       % 
% 	par.m2: mass of link 2 in [kg].                                       %
% 	par.m3: mass of link 3 in [kg].                                       %
% 	par.I1: moment of inertia of link 1 about its Center of Mass (CoM)    %
%           in [kg.m^2].                                                  %
% 	par.I2: moment of inertia of link 2 about its CoM in [kg.m^2].        %
% 	par.I3: moment of inertia of link 3 about its CoM in [kg.m^2].        %
% 	par.l1: length of link 1 in [m].                                      %
% 	par.l2: length of link 2 in [m].                                      %
% 	par.l3: length of link 3 in [m].                                      %
% 	par.lc1: distance from ankle joint to CoM of link 1 in [m].           % 
% 	par.lc2: distance from knee joint to CoM of link 2 in [m].            %
% 	par.lc3: distance from hip joint to CoM of link 3 in [m].             %
% 	par.g: acceleration of gravity [m/s^2].                               %
% 	par.K: feedback gain in synthetic input.                              %
% 	par.T: total simulation time.                                         %
% 	par.uconfig: Integer defining the configuration of input u.           %
% 	  Torques are in [N.m].                                               %
% 	  The forces are applied at the shoulder joint, and are in [N].       %
% 	  1: u=[ankle torque; shoulder torque; Fx; Fy].                       %
% 	  2: u=[knee torque; shoulder torque; Fx; Fy].                        %
% 	  3: u=[hip torque; shoulder torque; Fx; Fy].                         %
% 	  4: u=[ankle torque; knee torque; shoulder torque; Fx; Fy].          %
% 	  5: u=[ankle torque; hip torque; shoulder torque; Fx; Fy].           %
% 	  6: u=[knee torque; hip torque; shoulder torque; Fx; Fy].            %
% 	  7: u=[ankle torque;knee torque;hip torque; shoulder torque; Fx; Fy].%
% 	  8: u=[ankle torque; knee torque; hip torque; shoulder torque].      %
% 	  9: u=[ankle torque; knee torque; hip torque].                       %
%   par.wa: weight for ankle torque in control allocation.                %
%   par.wk: weight for knee torque in control allocation.                 %
%   par.wh: weight for hip torque in control allocation.                  %
%   par.ws: weight for shoulder torque in control allocation.             %
%   par.wFx: weight for horizontal force at the shoulder.                 %
%   par.wFy: weight for vertical force at the shoulder.                   %
%                                                                         %
% Output                                                                  %
%                                                                         %
% U: n by nu array with input applied by the controller.                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function U = FLThreeLinkInputs(x,xd,par)

% Parameters of the system.
lb = par.lb;     % Lower bounds for input.
ub = par.ub;     % Upper bounds for input.
m1 = par.m1;     % Mass of link 1 in [kg].
m2 = par.m2;     % Mass of link 2 in [kg].
m3 = par.m3;     % Mass of link 2 in [kg].
I1 = par.I1;     % Moment of inertia of link 1 about its CoM in [kg.m^2].
I2 = par.I2;     % Moment of inertia of link 2 about its CoM in [kg.m^2].
I3 = par.I3;     % Moment of inertia of link 3 about its CoM in [kg.m^2].
l1 = par.l1;     % Length of link 1 in [m].
l2 = par.l2;     % Length of link 2 in [m].
l3 = par.l3;     % Length of link 3 in [m].
lc1 = par.lc1;   % Distance from ankle joint to CoM of link 1 in [m].
lc2 = par.lc2;   % Distance from knee joint to CoM of link 2 in [m].
lc3 = par.lc3;   % Distance from hip joint to CoM of link 3 in [m].
g = par.g;       % Acceleration of gravity in [m/s^2].
K = par.K;       % Feedback gain in synthetic input.
uconfig = par.uconfig; % Configuration of input.

% Weights for control allocation problem.
wa = par.wa;
wk = par.wk;
wh = par.wh;
ws = par.ws;
wFx = par.wFx;
wFy = par.wFy;

% Terms of generalized force matrix, control allocation weights, 
% and bounds depending on input configuration.
switch uconfig
    case 1
        % Upper body loads, and actuation at the ankle.
        Loads = 1;                 % Presence of upper body loads.
        lbcfg = [lb(1);lb(4:end)]; % Lower bounds on controlled input.
        ubcfg = [ub(1);ub(4:end)]; % Upper bounds on controlled input.
        Atau1 = [1;0;0];           % Generalized forces term.
        W = diag([wa,ws,wFx,wFy]); % Control allocation weights.
    case 2
        % Upper body loads, and actuation at the knee.
        Loads = 1;                 % Presence of upper body loads.
        lbcfg = [lb(2);lb(4:end)]; % Lower bounds on controlled input.
        ubcfg = [ub(2);ub(4:end)]; % Upper bounds on controlled input.
        Atau1 = [0;1;0];           % Generalized force term.
        W = diag([wk,ws,wFx,wFy]); % Control allocation weights.
    case 3
        % Upper body loads, and actuation at the hip.
        Loads = 1;                 % Presence of upper body loads.
        lbcfg = [lb(3);lb(4:end)]; % Lower bounds on controlled input.
        ubcfg = [ub(3);ub(4:end)]; % Upper bounds on controlled input.
        Atau1 = [0;0;1];           % Generalized force term.
        W = diag([wh,ws,wFx,wFy]); % Control allocation weights.
    case 4
        % Upper body loads, and actuation at the ankle, and knee.
        Loads = 1;                    % Presence of upper body loads.
        lbcfg = [lb(1:2);lb(4:end)];  % Lower bounds on controlled input.
        ubcfg = [ub(1:2);ub(4:end)];  % Upper bounds on controlled input.
        Atau1 = [1 0;0 1;0 0];        % Generalized force term.
        W = diag([wa,wk,ws,wFx,wFy]); % Control allocation weights.
    case 5
        % Upper body loads, and actuation at the ankle, and hip.
        Loads = 1;                    % Presence of upper body loads.
        lbcfg = [lb(1);lb(3:end)];    % Lower bounds on controlled input.
        ubcfg = [ub(1);ub(3:end)];    % Upper bounds on controlled input.
        Atau1 = [1 0; 0 0; 0 1];      % Generalized force term.
        W = diag([wa,wh,ws,wFx,wFy]); % Control allocation weights.
    case 6
        % Upper body loads, and actuation at the knee, and hip.
        Loads = 1;                    % Presence of upper body loads.
        lbcfg = lb(2:end);            % Lower bounds on controlled input.
        ubcfg = ub(2:end);            % Upper bounds on controlled input.
        Atau1 = [0 0; 1 0; 0 1];      % Generalized force term.
        W = diag([wk,wh,ws,wFx,wFy]); % Control allocation weights.
    case 7
        % Upper body loads, and actuation at the ankle, knee, and hip.
        Loads = 1;                     % Presence of upper body loads.
        lbcfg = lb;                    % Lower bounds on controlled input.
        ubcfg = ub;                    % Upper bounds on controlled input.
        Atau1 = eye(3);                % Generalized force term.
        W = diag([wa,wk,wh,ws,wFx,wFy]); % Control allocation weights.
    case 8
        % Actuation at the ankle, knee, hip, and shoulder.
        Loads = 0;                   % Presence of upper body loads.
        lbcfg = lb(1:4);             % Lower bounds on controlled input.
        ubcfg = ub(1:4);             % Upper bounds on controlled input.
        Atau1 = [eye(3),-ones(3,1)]; % Generalized force matrix.
        W = diag([wa,wk,wh,ws]);     % Control allocation weights.
    case 9
        % Actuation at the ankle, knee, and hip.
        Loads = 0;            % Presence of upper body loads.
        lbcfg = lb(1:3);      % Lower bounds on controlled input.
        ubcfg = ub(1:3);      % Upper bounds on controlled input.
        Atau1 = eye(3);       % Generalized force matrix.
        W = diag([wa,wk,wh]); % Control allocation weights.
    otherwise
        error('Unrecognized configuration of input. Options range from 1 to 9.')
end

% Number of time samples in simulation data.
n = size(x,1);

% Array for storing input.
U = zeros(size(x,1),size(W,1));

for i=1:n
    % State of the system.
    th = x(i,1:3);
    om = x(i,4:6);
    omsq = [om(1);om(1)+om(2);om(1)+om(2)+om(3)].^2;
    
    % Trigonometric functions.
    s2 = sin(th(2));
    s3 = sin(th(3));
    s23 = sin(th(2)+th(3));
    c2 = cos(th(2));
    c3 = cos(th(3));
    c23 = cos(th(2)+th(3));
    sv = sin([th(1);th(1)+th(2);th(1)+th(2)+th(3)]);
    cv = cos([th(1);th(1)+th(2);th(1)+th(2)+th(3)]);
    
    % Mass matrix.
    M11 = I1 + I2 + I3 + lc1^2*m1 + m2*(l1^2 + 2*l1*lc2*c2 + lc2^2) + m3*(l1^2 + 2*l1*l2*c2 + 2*l1*lc3*c23 + l2^2 + 2*l2*lc3*c3 + lc3^2);
    M12 = I2 + I3 + lc2*m2*(l1*c2 + lc2) + m3*(l1*l2*c2 + l1*lc3*c23 + l2^2 + 2*l2*lc3*c3 + lc3^2);
    M13 = I3 + lc3*m3*(l1*c23 + l2*c3 + lc3);
    M22 = I2 + I3 + lc2^2*m2 + m3*(l2^2 + 2*l2*lc3*c3 + lc3^2);
    M23 = I3 + lc3*m3*(l2*c3 + lc3);
    M33 = I3 + lc3^2*m3;
    M = [M11,M12,M13;M12,M22,M23;M13,M23,M33];
    
    % Constant terms.
    % k0 = 1/(m1+m2+m3);
    k1 = lc1*m1+l1*m2+l1*m3;
    k2 = lc2*m2+l2*m3;
    k3 = lc3*m3;
    
    % Vector of gravity effects.
    fg = -g*[k1,k2,k3;0,k2,k3;0,0,k3]*cv;
    
    % Coriolis effects matrix.
    k4 = l1*(k2*s2+k3*s23);
    k5 = k3*l2*s3;
    C = [k4,-k2*l1*s2+k3*l2*s3,-k3*l1*s23-k3*l2*s3;k4,k5,-k5;l1*k3*s23,k5,0];
    
    % Gravity, and Coriolis effects.
    f = fg-C*omsq;
    
    % Matrix of generalized force.
    if Loads
        Atau = [Atau1,-ones(3,1),-[l1,l2,l3;0,l2,l3;0 0 l3]*sv,[l1,l2,l3;0,l2,l3;0 0 l3]*cv];
    else
        Atau = Atau1;
    end
    
    % Synthetic input.
    v = -K*(x(i,:)'-xd(i,1:6)');
    
    % Feedback linearization of the system.
    b = M*(v+xd(i,7:end)')-f;
    
    % Solve constrained linear least-squares problem for control allocation.
    options = optimoptions('lsqlin','Algorithm','active-set','Display','off');
    u = lsqlin(W,zeros(size(W,1),1),[],[],Atau,b,lbcfg,ubcfg,[],options);
    
    % Store input.
    U(i,:) = u';
end