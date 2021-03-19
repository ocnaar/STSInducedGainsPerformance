%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez-Aroche                                                  %
% Berkeley Center for Control and Identification                          %
% Summer 2017                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Obtain reference trajectories for a sit-to-stand movement in the space  %
% of z. The trajectories for the angular position of link 2 relative to   %
% link 1, and the position of the Center of Mass (CoM) of the robot are   %
% interpolated for the time values in array t using polynomials, whose    %
% coefficients are specified in array A. The time derivatives of the      %
% polynomials are used for the velocity, and acceleration trajectories.   %
%                                                                         %
% Input                                                                   %
%                                                                         %
% A: 3 by m array of polynomial coefficients ordered by descending powers.%
% 	A(1,:): polynomial coefficients for interpolation of the angular      %
% 	  position of link 2 relative to link 1.                              %
% 	A(2,:): polynomial coefficients for interpolation of the x coordinate %
% 	  of the position of the robot CoM.                                   %
% 	A(3,:): polynomial coefficients for interpolation of the y coordinate %
% 	  of the position of the robot CoM.                                   %
% t: n by 1 time array with values between 0 and T in [s].                %
% T: final simulation time in [s].                                        %
% zi: initial configuration of the three-link robot in the space of z.    %
% 	zi(1): angular position of link 2 relative to link 1 in [rad].        %
% 	zi(2): x coordinate of the position of the robot CoM in [m].          %
% 	zi(3): y coordinate of the position of the robot CoM in [m].          %
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
%                                                                         %
% Output                                                                  %
%                                                                         %
% Z: 9 by n array.                                                        %
% 	Z(1,:): angular position of link 2 relative to link 1 in [rad].       %
% 	Z(2,:): x coordinate of the position of the CoM for the robot in [m]. %
% 	Z(3,:): y coordinate of the position of the CoM for the robot in [m]. %
% 	Z(4,:): angular velocity of link 2 in [rad/s].                        %
% 	Z(5,:): x coordinate of the velocity of the CoM in [m/s].             %
% 	Z(6,:): y coordinate of the velocity of the CoM in [m/s].             %
% 	Z(7,:): angular acceleration of link 2 in [rad/s^2].                  %
% 	Z(8,:): x coordinate of the acceleration of the CoM in [m/s^2].       %
% 	Z(9,:): y coordinate of the acceleration of the CoM in [m/s^2].       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Z = desiredz(A,t,T,zi,zf)

% Number of time steps.
nt = length(t);

% Number of polynomials parameterizing reference trajectories.  
np = size(A,1);

% Array for storing desired trajectories.
Z = zeros(3*np,nt);

for j=1:nt
    % Normalize time.
    tn = t(j)/T; 
    for i=1:np
        % Vector of coefficients for trajectory.
        coef = A(i,:);
        m = length(coef);
        p = m-1:-1:0;
        b = (tn*ones(m,1)).^(p');
        phi = coef*b;
        
        % Vector of coefficients for first time derivative of trajectory.
        coefp = coef(1:end-1).*p(1:end-1);
        phip = coefp*b(2:end);
        
        % Vector of coefficients for second time derivative of trajectory.
        coefpp = coefp(1:end-1).*p(2:end-1);
        phipp = coefpp*b(3:end);
        
        % Reference trajectory in the space of z.
        Z(i,j) = zi(i)+(zf(i)-zi(i))*phi;
        Z(np+i,j) = (zf(i)-zi(i))*phip*(1/T);
        Z(2*np+i,j) = (zf(i)-zi(i))*phipp*(1/T^2);
    end
end