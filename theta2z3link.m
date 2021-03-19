%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez-Aroche                                                  %
% Berkeley Center for Control and Identification                          %
% Spring 2017                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Map n vectors in the space of theta into the space of z for a three-link%
% planar robot with revolute joints.                                      %
%                                                                         %
% Input                                                                   %
%                                                                         %
% Theta: 9 by n array.                                                    %
% 	Theta(1,:): angular position of link 1 relative to the horizontal[rad]% 
% 	Theta(2,:): angular position of link 2 relative to link 1 in [rad].   %
% 	Theta(3,:): angular position of link 3 relative to link 2 in [rad].   %
% 	Theta(4,:): angular velocity of link 1 in [rad/s].                    %
% 	Theta(5,:): angular velocity of link 2 in [rad/s].                    %
% 	Theta(6,:): angular velocity of link 3 in [rad/s].                    %
% 	Theta(7,:): angular acceleration of link 1 in [rad/s^2].              %
% 	Theta(8,:): angular acceleration of link 2 in [rad/s^2].              %
% 	Theta(9,:): angular acceleration of link 3 in [rad/s^2].              %
%                                                                         %
% par: structure containing the parameters of the three-link robot.       %
% 	par.l1: length of link 1 in [m].                                      %
% 	par.lc1: distance from ankle to Center of Mass (CoM) of link 1 in [m].%
% 	par.l2: length of link 2 in [m].                                      %
% 	par.lc2: distance from knee joint to CoM of link 2 in [m].            %
% 	par.lc3: distance from hip joint to CoM of link 3 in [m].             %
% 	par.m1: mass of link 1 in [kg].                                       %
% 	par.m2: mass of link 2 in [kg].                                       %
% 	par.m3: mass of link 3 in [kg].                                       %
%                                                                         %
% Output                                                                  %
%                                                                         %
% Z: 9 by n array.                                                        %
% 	Z(1,:): angular position of link 2 relative to link 1 in [rad].       %
% 	Z(2,:): x coordinate of the position of the robot CoM in [m].         %
% 	Z(3,:): y coordinate of the position of the robot CoM in [m].         %
% 	Z(4,:): angular velocity of link 2 in [rad/s].                        %
% 	Z(5,:): x coordinate of the velocity of the robot CoM in [m/s].       %
% 	Z(6,:): y coordinate of the velocity of the robot CoM in [m/s].       %
% 	Z(7,:): angular acceleration of link 2 in [rad/s^2].                  %
% 	Z(8,:): x coordinate of the acceleration of the robot CoM in [m/s^2]. %
% 	Z(9,:): y coordinate of the acceleration of the robot CoM in [m/s^2]. %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Z = theta2z3link(Theta,par)

% Parameters of the three-link robot.
l1 = par.l1;        % Length of link 1 [m].
lc1 = par.lc1;      % Distance from ankle joint to CoM of link 1 [m].
l2 = par.l2;        % Length of link 2 [m].
lc2 = par.lc2;      % Distance from knee joint to CoM of link 2 [m].
lc3 = par.lc3;      % Distance from hip joint to CoM of link 3 [m].
m1 = par.m1;        % Mass of link 1 [kg].
m2 = par.m2;        % Mass of link 2 [kg].
m3 = par.m3;        % Mass of link 3 [kg].

% Constant terms.
k0 = 1/(m1+m2+m3);
k1 = lc1*m1+l1*m2+l1*m3;
k2 = lc2*m2+l2*m3;
k3 = lc3*m3;

% Number of vectors in the space of theta. 
n = size(Theta,2);

% Array for storing vectors in the space of z.
Z = zeros(9,n);

% Map vectors in the space of theta into the space of z.
for i=1:n
    % Angular positions of the links. 
    th1 = Theta(1,i);
    th2 = Theta(2,i);
    th3 = Theta(3,i);
    
    % Sine and cosine functions.
    s1 = sin(th1);
    s12 = sin(th1+th2);
    s123 = sin(th1+th2+th3);
    c1 = cos(th1);
    c12 = cos(th1+th2);
    c123 = cos(th1+th2+th3);
    
    % Calculate position coordinates of the CoM.
    x = k0*(k1*c1+k2*c12+k3*c123);
    y = k0*(k1*s1+k2*s12+k3*s123);
    
    % Angular velocities of the links. 
    om1 = Theta(4,i);
    om2 = Theta(5,i);
    om3 = Theta(6,i);
    
    % Calculate velocity of the CoM.
    vx = -om1*y-om2*k0*(k2*s12+k3*s123)-om3*k0*k3*s123;
    vy = om1*x+om2*k0*(k2*c12+k3*c123)+om3*k0*k3*c123;
    
    % Angular accelerations of the links.
    al1 = Theta(7,i);
    al2 = Theta(8,i);
    al3 = Theta(9,i);
    
    % Square of angular velocities.
    omsq = [om1, om2, om3].^2;
    
    % Calculate acceleration of the CoM. 
    ax = -omsq*[x;k0*(k2*c12+k3*c123);k0*k3*c123];
    ax = ax-2*om1*om2*k0*(k2*c12+k3*c123)-2*(om1+om2)*om3*k0*k3*c123;
    ax = ax-al1*y-al2*k0*(k2*s12+k3*s123)-al3*k0*k3*s123;
    
    ay = -omsq*[y;k0*(k2*s12+k3*s123);k0*k3*s123];
    ay = ay-2*om1*om2*k0*(k2*s12+k3*s123)-2*(om1+om2)*om3*k0*k3*s123;
    ay = ay+al1*x+al2*k0*(k2*c12+k3*c123)+al3*k0*k3*c123;
    
    % Store vector in the space of z.
    Z(:,i) = [th2;x;y;om2;vx;vy;al2;ax;ay];
end