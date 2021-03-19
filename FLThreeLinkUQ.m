%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez-Aroche                                                  %
% Berkeley Center for Control and Identification                          %
% Spring 2017                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Closed-loop dynamics of a three-link planar robot controlled with       %
% feedback linearization, and control allocation. The reference trajectory%
% is built with a cubic polynomial satisfying zero-slope conditions at the%
% initial, and final configurations desired for the robot in the space of % 
% z: zi, and zf, respectively.                                            %
%                                                                         %
% Input                                                                   %
%                                                                         %
% t: time in [s].                                                         %
% x: state of the system.                                                 %
% 	x(1): angular position of link 1 relative to the horizontal in [rad]. %
% 	x(2): angular position of link 2 relative to link 1 in [rad].         % 
% 	x(3): angular position of link 3 relative to link 2 in [rad].         %
% 	x(4): angular velocity of link 1 in [rad/s].                          %
% 	x(5): angular velocity of link 2 in [rad/s].                          %
% 	x(6): angular velocity of link 3 in [rad/s].                          %
% zi: initial configuration of the three-link robot in the space of z.    %
% 	zi(1): angular position of link 2 relative to link 1 in [rad].        %
% 	zi(2): x coordinate of the position of the Center of Mass (CoM) in [m]%
% 	zi(3): y coordinate of the position of the CoM for the robot in [m].  %
% 	zi(4): angular velocity of link 2 in [rad/s].                         %
% 	zi(5): x coordinate of the velocity of the CoM in [m/s].              %
% 	zi(6): y coordinate of the velocity of the CoM in [m/s].              %
% 	zi(7): angular acceleration of link 2 in [rad/s^2].                   %
% 	zi(8): x coordinate of the acceleration of the CoM in [m/s^2].        %
% 	zi(9): y coordinate of the acceleration of the CoM in [m/s^2].        %
% zf: final configuration of the three-link robot in the space of z.      %
% 	zf(1): angular position of link 2 relative to link 1 in [rad].        %
% 	zf(2): x coordinate of the position of the CoM for the robot in [m].  %
% 	zf(3): y coordinate of the position of the CoM for the robot in [m].  %
% 	zf(4): angular velocity of link 2 in [rad/s].                         %
% 	zf(5): x coordinate of the velocity of the CoM in [m/s].              %
% 	zf(6): y coordinate of the velocity of the CoM in [m/s].              %
% 	zf(7): angular acceleration of link 2 in [rad/s^2].                   %
% 	zf(8): x coordinate of the acceleration of the CoM in [m/s^2].        %
% 	zf(9): y coordinate of the acceleration of the CoM in [m/s^2].        %
% par: structure with nominal parameters of the system, and controller    %
% settings.                                                               %
% 	par.P: coefficients of cubic polynomials for reference trajectories.  % 
% 	par.lb: lower bounds for system input in control allocation.          % 
% 	par.ub: upper bounds for system input in control allocation.          %
% 	par.m1: mass of link 1 in [kg].                                       % 
% 	par.m2: mass of link 2 in [kg].                                       %
% 	par.m3: mass of link 3 in [kg].                                       %
% 	par.I1: moment of inertia of link 1 about its CoM in [kg.m^2].        %
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
% UQpar: structure with values for the parameter of the system which are  %
% unknown to the feedback linearization controller.                       %
% 	UQpar.m1: mass of link 1 in [kg].                                     %
% 	UQpar.m2: mass of link 2 in [kg].                                     %
% 	UQpar.m3: mass of link 3 in [kg].                                     %
% 	UQpar.I1: moment of inertia of link 1 about its CoM in [kg.m^2].      %
% 	UQpar.I2: moment of inertia of link 2 about its CoM in [kg.m^2].      %
% 	UQpar.I3: moment of inertia of link 3 about its CoM in [kg.m^2].      %
% 	UQpar.l1: length of link 1 in [m].                                    %
% 	UQpar.l2: length of link 2 in [m].                                    %
% 	UQpar.l3: length of link 3 in [m].                                    %
% 	UQpar.lc1: distance from ankle joint to CoM of link 1 in [m].         %
% 	UQpar.lc2: distance from knee joint to CoM of link 2 in [m].          %
% 	UQpar.lc3: distance from hip joint to CoM of link 3 in [m].           %
%                                                                         %
% Output                                                                  %
%                                                                         %
% xdot: time derivative of the state.                                     %
% 	xdot(1): angular velocity of link 1 in [rad/s].                       %
%   xdot(2): angular velocity of link 2 in [rad/s].                       %
%   xdot(3): angular velocity of link 3 in [rad/s].                       %
%   xdot(4): angular acceleration of link 1 in [rad/s^2].                 %
%   xdot(5): angular acceleration of link 2 in [rad/s^2].                 %
%   xdot(6): angular acceleration of link 3 in [rad/s^2].                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xdot = FLThreeLinkUQ(t,x,zi,zf,par,UQpar)

% Nominal parameter values of the system used by the controller. 
P = par.P;       % Matrix with coefficients of polynomials for reference trajectories. 
lb = par.lb;     % Lower bound for system input. 
ub = par.ub;     % Upper bound for system input. 
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
T = par.T;       % Total simulation time in [s]. 
uconfig = par.uconfig; % Integer to specify input configuration.

% Weights for input in control allocation problem.
wa = par.wa;
wk = par.wk;
wh = par.wh;
ws = par.ws;
wFx = par.wFx;
wFy = par.wFy;

% Reference trajectory in the space of z.
z = desiredz(P,t,T,zi(1:3),zf(1:3));

% Reference trajectory in the space of theta. 
thd = z2theta3link(z,par);

% System state.
th = x(1:3);
om = x(4:6);

% Square of angular velocity terms.
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

% Terms of generalized force matrix, control allocation weights, and input
% bounds depending on configuration. 
switch uconfig
    case 1
        % Upper body loads, and actuation at the ankle.
        Loads = 1;                 % Presence of upper body loads. 
        lbcfg = [lb(1);lb(4:end)]; % Lower bounds on controlled inputs.
        ubcfg = [ub(1);ub(4:end)]; % Upper bounds on controlled inputs.
        Atau1 = [1;0;0];           % Generalized force term. 
        W = diag([wa,ws,wFx,wFy]); % Control allocation weights.
    case 2
        % Upper body loads, and actuation at the knee.
        Loads = 1;                 % Presence of upper body loads.
        lbcfg = [lb(2);lb(4:end)]; % Lower bounds on controlled inputs.
        ubcfg = [ub(2);ub(4:end)]; % Upper bounds on controlled inputs.
        Atau1 = [0;1;0];           % Generalized force term. 
        W = diag([wk,ws,wFx,wFy]); % Control allocation weights.
    case 3
        % Upper body loads and actuation at the hip.
        Loads = 1;                 % Presence of upper body loads.
        lbcfg = [lb(3);lb(4:end)]; % Lower bounds on controlled inputs.
        ubcfg = [ub(3);ub(4:end)]; % Upper bounds on controlled inputs.
        Atau1 = [0;0;1];           % Generalized force term. 
        W = diag([wh,ws,wFx,wFy]); % Control allocation weights.
    case 4
        % Upper body loads and actuation at the ankle and knee.
        Loads = 1;                    % Presence of upper body loads.
        lbcfg = [lb(1:2);lb(4:end)];  % Lower bounds on controlled inputs.
        ubcfg = [ub(1:2);ub(4:end)];  % Upper bounds on controlled inputs.
        Atau1 = [1 0;0 1;0 0];        % Generalized force term. 
        W = diag([wa,wk,ws,wFx,wFy]); % Control allocation weights.
    case 5
        % Upper body loads and actuation at the ankle and hip.
        % u = [ankle torque, hip torque, shoulder torque, Fx, Fy]'
        Loads = 1;                    % Presence of upper body loads.
        lbcfg = [lb(1);lb(3:end)];    % Lower bounds on controlled inputs.
        ubcfg = [ub(1);ub(3:end)];    % Upper bounds on controlled inputs.
        Atau1 = [1 0; 0 0; 0 1];      % Generalized forces term. 
        W = diag([wa,wh,ws,wFx,wFy]); % Control allocation weights.
    case 6
        % Upper body loads and actuation at the knee and hip.
        Loads = 1;                    % Presence of upper body loads.
        lbcfg = lb(2:end);            % Lower bounds on controlled inputs.
        ubcfg = ub(2:end);            % Upper bounds on controlled inputs.
        Atau1 = [0 0; 1 0; 0 1];      % Generalized force term.
        W = diag([wk,wh,ws,wFx,wFy]); % Control allocation weights.
    case 7
        % Upper body loads and actuation at the ankle, knee and hip.
        Loads = 1;                     % Presence of upper body loads.
        lbcfg = lb;                    % Lower bounds on controlled inputs.
        ubcfg = ub;                    % Upper bounds on controlled inputs.
        Atau1 = eye(3);                % Generalized force term. 
        W = diag([wa,wk,wh,ws,wFx,wFy]); % Control allocation weights.
    case 8
        % Actuation at the ankle, knee, hip, and shoulder.
        Loads = 0;                   % Presence of upper body loads.
        lbcfg = lb(1:4);             % Lower bounds on controlled inputs.
        ubcfg = ub(1:4);             % Upper bounds on controlled inputs.
        Atau1 = [eye(3),-ones(3,1)]; % Generalized force matrix.
        W = diag([wa,wk,wh,ws]);     % Control allocation weights.
    case 9
        % Actuation at the ankle, knee, and hip.
        Loads = 0;            % Presence of upper body loads.
        lbcfg = lb(1:3);      % Lower bounds on controlled inputs.
        ubcfg = ub(1:3);      % Upper bounds on controlled inputs.
        Atau1 = eye(3);       % Generalized force matrix.
        W = diag([wa,wk,wh]); % Control allocation weights.
    otherwise
        error('Unrecognized configuration of inputs. Options range from 1 to 9.')
end

% Matrix of generalized force.
if Loads
   Atau = [Atau1,-ones(3,1),-[l1,l2,l3;0,l2,l3;0 0 l3]*sv,[l1,l2,l3;0,l2,l3;0 0 l3]*cv];
else
   Atau = Atau1;  
end

% Synthetic input. 
v = -K*(x-thd(1:6));

% Feedback linearization of the system.
b = M*(v+thd(7:end))-f;

% Solve constrained linear least-squares problem for control allocation.
options = optimoptions('lsqlin','Algorithm','active-set','Display','off');
u = lsqlin(W,zeros(size(W,1),1),[],[],Atau,b,lbcfg,ubcfg,[],options);

% Values for the parameter of the system which are unknown to the controller.
m1 = UQpar.m1;     % Mass of link 1 in [kg]. 
m2 = UQpar.m2;     % Mass of link 2 in [kg]. 
m3 = UQpar.m3;     % Mass of link 3 in [kg]. 
I1 = UQpar.I1;     % Moment of inertia of link 1 about its CoM in [kg.m^2].
I2 = UQpar.I2;     % Moment of inertia of link 2 about its CoM in [kg.m^2].
I3 = UQpar.I3;     % Moment of inertia of link 3 about its CoM in [kg.m^2].
l1 = UQpar.l1;     % Length of link 1 in [m].
l2 = UQpar.l2;     % Length of link 2 in [m].
l3 = UQpar.l3;     % Length of link 3 in [m].
lc1 = UQpar.lc1;   % Distance from ankle joint to CoM of link 1 in [m]. 
lc2 = UQpar.lc2;   % Distance from knee joint to CoM of link 2 in [m].
lc3 = UQpar.lc3;   % Distance from hip joint to CoM of link 3 in [m].

% Mass matrix.
M11 = I1 + I2 + I3 + lc1^2*m1 + m2*(l1^2 + 2*l1*lc2*c2 + lc2^2) + m3*(l1^2 + 2*l1*l2*c2 + 2*l1*lc3*c23 + l2^2 + 2*l2*lc3*c3 + lc3^2);
M12 = I2 + I3 + lc2*m2*(l1*c2 + lc2) + m3*(l1*l2*c2 + l1*lc3*c23 + l2^2 + 2*l2*lc3*c3 + lc3^2);
M13 = I3 + lc3*m3*(l1*c23 + l2*c3 + lc3);
M22 = I2 + I3 + lc2^2*m2 + m3*(l2^2 + 2*l2*lc3*c3 + lc3^2);
M23 = I3 + lc3*m3*(l2*c3 + lc3);
M33 = I3 + lc3^2*m3;
M = [M11,M12,M13;M12,M22,M23;M13,M23,M33];
Minv = M\eye(3);

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

% Time derivatives of state vector.
xdot = zeros(6,1);
xdot(1:3) = om;
xdot(4:6) = Minv*(f+Atau*u);