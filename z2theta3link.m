%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez-Aroche                                                  %
% Berkeley Center for Control and Identification                          %
% Spring 2017                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Map vectors in the space of z into the space of theta for a three-link  %
% planar robot. For the complete algorithm please refer to:               %
%                                                                         %
% O.Narvaez-Aroche, A. Packard, and M. Arcak, “Motion planning of the     %
% sit-to-stand movement for powered lower limb orthoses”, ASME 2017       %
% Dynamic Systems & Control Conference, October 2017. Tysons, VA, USA.    % 
% http://dx.doi.org/10.1115/DSCC2017-5289                                 %
%                                                                         %
% Input                                                                   %
%                                                                         %
% Z: 9 by n array.                                                        %
% 	Z(1,:): angular position of link 2 relative to link 1 in [rad].       %
% 	Z(2,:): x coordinate of the position of the robot Center of Mass (CoM)%
%           in [m].                                                       %
% 	Z(3,:): y coordinate of the position of the robot CoM in [m].         %
% 	Z(4,:): angular velocity of link 2 in [rad/s].                        %
% 	Z(5,:): x coordinate of the velocity of the robot CoM in [m/s].       %
% 	Z(6,:): y coordinate of the velocity of the robot CoM in [m/s].       %
% 	Z(7,:): angular acceleration of link 2 in [rad/s^2].                  %
% 	Z(8,:): x coordinate of the acceleration of the robot CoM in [m/s^2]. %
% 	Z(9,:): y coordinate of the acceleration of the robot CoM in [m/s^2]. %
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Theta = z2theta3link(Z,par)

% Parameters of the three-link robot.
l1 = par.l1;        % Length of link 1 in [m].
lc1 = par.lc1;      % Distance from ankle joint to CoM of link 1 in [m].
l2 = par.l2;        % Length of link 2 in [m].
lc2 = par.lc2;      % Distance from knee joint to CoM of link 2 in [m].
lc3 = par.lc3;      % Distance from hip joint to CoM of link 3 in [m].
m1 = par.m1;        % Mass of link 1 in [kg].
m2 = par.m2;        % Mass of link 2 in [kg].
m3 = par.m3;        % Mass of link 2 in [kg].

% Constant terms.
k0 = 1/(m1+m2+m3);
k1 = lc1*m1+l1*m2+l1*m3;
k2 = lc2*m2+l2*m3;
k3 = lc3*m3;

% y coordinate of the CoM of the three-link robot in the vertical position.
yvert = k0*(k1+k2+k3);

% Number of vectors in the space of z. 
n = size(Z,2);

% Array for storing vectors in the theta space. 
Theta = zeros(9,n);

% Map vectors in the space of z into the space of theta. 
for i=1:n
    if Z(3,i)==yvert
        % Angles of the robot in the vertical position. 
        Theta(1,i) = pi/2;
        Theta(2:3,i) = 0;
        % Velocities and accelerations cannot be computed since matrix V
        % becomes singular. 
        Theta(4:end,i) = NaN;
    else
        % Compute angle th1.
        th2 = Z(1,i);
        alpha = pi+th2;
        if Z(2,i)>=0
            if Z(3,i)>=0
                beta = atan(Z(3,i)/Z(2,i));
            else
                beta = 2*pi+atan(Z(3,i)/Z(2,i));
            end
        else
            beta = pi+atan(Z(3,i)/Z(2,i));
        end
        % Squared norm of r1+r2. 
        c = k0^2*(k1^2+k2^2+2*k1*k2*cos(th2));
        varphi = asin(max(min(k0*k2*sin(alpha)/sqrt(c),1),-1));
		
        % Norm of the position of the CoM.
        a = norm(Z(2:3,i));
        phi = acos(max(min(((k0*k3)^2-c-a^2)/(-2*sqrt(c)*a),1),-1));
        th1 = beta-phi+varphi;
        
        % Compute angle th3.
        psi = asin(max(min(sqrt(c)*sin(phi)/(k0*k3),1),-1));
        th3 = psi-th1-th2+beta;
        
        % Compute angular velocities om1, and om3.
        om2 = Z(4,i);
        s12 = sin(th1+th2);
        s123 = sin(th1+th2+th3);
        c12 = cos(th1+th2);
        c123 = cos(th1+th2+th3);
        f1= Z(5:6,i)-om2*[-k0*(k2*s12+k3*s123); k0*(k2*c12+k3*c123)];
        V = [-Z(3,i), -k0*k3*s123; Z(2,i), k0*k3*c123];
        om = V\f1;
        
        % Compute angular accelerations al1, and al3.
        al2 = Z(7,i);
        q = [-k0*(k2*s12+k3*s123);k0*(k2*c12+k3*c123)]*al2;
        q = q-[Z(2,i),k0*(k2*c12+k3*c123),k0*k3*c123; Z(3,i),k0*(k2*s12+k3*s123),k0*k3*s123]*[om(1)^2;om2^2;om(2)^2];
        q = q-2*om(1)*om2*[k0*(k2*c12+k3*c123);k0*(k2*s12+k3*s123)];
        q = q-2*(om(1)+om2)*om(2)*[k0*k3*c123;k0*k3*s123];
        al = V\(Z(8:9,i)-q);
        
        % Vector in the space of theta. 
        Theta(:,i) = [th1;th2;th3;om(1);om2;om(2);al(1);al2;al(2)];
    end
end