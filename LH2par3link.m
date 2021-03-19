%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez-Aroche                                                  %
% Berkeley Center for Control and Identification                          %
% Summer 2017                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Use an array of random samples to obtain cells containing structures    %
% with values for the parameters of a three-link planar robot. The values %
% for the parameters lie within a known interval.                         %
%                                                                         %
% Input                                                                   %
%                                                                         %
% LH: n by 12 array of random samples between 0 and 1.                    %
% parmin: 12 by 1 vector with lower bound for the parameter uncertainty.  % 
% 	parmin(1): Mass of link 1 in [kg].                                    %
%   parmin(2): Mass of link 2 in [kg].                                    %
%   parmin(3): Mass of link 3 in [kg].                                    %
%   parmin(4): Moment of inertia of link 1 about its Center of Mass (CoM) %
%   	in [kg.m^2].                                                      %
%   parmin(5): Moment of inertia of link 2 about its CoM in [kg.m^2].     %
%   parmin(6): Moment of inertia of link 3 about its CoM in [kg.m^2].     %
%   parmin(7): Length of link 1 in [m].                                   %
%   parmin(8): Length of link 2 in [m].                                   %
%   parmin(9): Length of link 3 in [m].                                   %
%   parmin(10): Distance from ankle joint to CoM of link 1 in [m].        %
%   parmin(11): Distance from knee joint to CoM of link 2 in [m].         %
%   parmin(12): Distance from hip joint to CoM of link 3 in [m].          %
% parmax: 12 by 1 vector with upper bound for the parameter uncertainty.  %
%                                                                         %
% Output                                                                  %
%                                                                         %
% UQpar: 1 by n cell array with structures containing the values for the  %
% parameters of the three-link robot.                                     %
% 	UQpar{i}.g: Acceleration of gravity [m/s^2].                          %
%   UQpar{i}.m1: Mass of link 1 in [kg].                                  %
%   UQpar{i}.m2: Mass of link 2 in [kg].                                  %
%   UQpar{i}.m3: Mass of link 3 in [kg].                                  %
%   UQpar{i}.I1: Moment of inertia of link 1 about its Center of Mass     %
%   	(CoM) in [kg.m^2].                                                %
%   UQpar{i}.I2: Moment of inertia of link 2 about its CoM in [kg.m^2].   %
%   UQpar{i}.I3: Moment of inertia of link 3 about its CoM in [kg.m^2].   %
%   UQpar{i}.l1: Length of link 1 in [m].                                 %
%   UQpar{i}.l2: Length of link 2 in [m].                                 %
%   UQpar{i}.l3: Length of link 3 in [m].                                 %
%   UQpar{i}.lc1: Distance from ankle joint to CoM of link 1 in [m].      %
%   UQpar{i}.lc2: Distance from knee joint to CoM of link 2 in [m].       %
%   UQpar{i}.lc3: Distance from hip joint to CoM of link 3 in [m].        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UQpar = LH2par3link(LH,parmin,parmax)

% Number of sampled values for system parameters. 
ne = size(LH,1);

% Cell for storing structures with parameters of the system.
UQpar = cell(1,ne);

for i=1:ne
    % Parameter values from Latin Hypercube Sampling and bounds.  
    realpar = parmin + (parmax-parmin).*LH(i,:)';
    UQpar{i}.g = 9.81;            % Acceleration of gravity in [m/s^2].
    UQpar{i}.m1 = realpar(1);     % Mass of link 1 in [kg].
    UQpar{i}.m2 = realpar(2);     % Mass of link 2 in [kg].
    UQpar{i}.m3 = realpar(3);     % Mass of link 3 in [kg].
    UQpar{i}.I1 = realpar(4);     % Moment of inertia of link 1 about its CoM in [kg.m^2].
    UQpar{i}.I2 = realpar(5);     % Moment of inertia of link 2 about its CoM in [kg.m^2].
    UQpar{i}.I3 = realpar(6);     % Moment of inertia of link 3 about its CoM in [kg.m^2].
    UQpar{i}.l1 = realpar(7);     % Length of link 1 in [m].
    UQpar{i}.l2 = realpar(8);     % Length of link 2 in [m].
    UQpar{i}.l3 = realpar(9);     % Length of link 3 in [m].
    UQpar{i}.lc1 = realpar(10);   % Distance from ankle joint to CoM of link 1 in [m].
    UQpar{i}.lc2 = realpar(11);   % Distance from knee joint to CoM of link 2 in [m].
    UQpar{i}.lc3 = realpar(12);   % Distance from hip joint to CoM of link 3 in [m].
end