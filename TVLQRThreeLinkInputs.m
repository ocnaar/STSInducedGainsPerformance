%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez-Aroche                                                  %
% Berkeley Center for Control and Identification                          %
% Spring 2017                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Compute input applied to a three-link robot by a finite time horizon    %
% LQR controller from its closed-loop simulation data.                    %
%                                                                         %
% Input                                                                   %
%                                                                         %
% T: n by 1 time array for closed-loop simulation.                        %
% X: n by 6 array with state trajectory from closed-loop simulation.      % 
% 	X(:,1): angular position of link 1 relative to the horizontal in [rad]%
% 	X(:,2): angular position of link 2 relative to link 1 in [rad].       % 
% 	X(:,3): angular position of link 3 relative to link 2 in [rad].       %
% 	X(:,4): angular velocity of link 1 in [rad/s].                        %
% 	X(:,5): angular velocity of link 2 in [rad/s].                        %
% 	X(:,6): angular velocity of link 3 in [rad/s].                        %
% Ktvmat: nu by 6 LTV object with LQR gain.                               %
% tgrid: m by 1 time array for samples in reference trajectories.         %
% xref: m by 6 array with samples from reference state trajectory.        %
% uref: m by nu array with samples from reference input trajectory.       %
%                                                                         %
% Output                                                                  %
%                                                                         %
% U: n by nu array with input applied by the finite time horizon LQR      %
% controller.                                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function U = TVLQRThreeLinkInputs(T,X,Ktvmat,tgrid,xref,uref)

% Number of states.
nx = size(xref,2);

% Number of inputs.
nu = size(uref,2);

% Time samples in simulation.
nt = length(T);

% Interpolation of reference state trajectories.
Xbar = zeros(nt,nx);
for k=1:nx
    Xbar(:,k) = csapi(tgrid,xref(:,k),T);
end

% Interpolation of reference input trajectories.
Ubar = zeros(nt,nu);
for k=1:nu
    Ubar(:,k) = csapi(tgrid,uref(:,k),T);
end

% Compute LQR input from simulation data.
dx = X(:,1:nx)-Xbar;
U = zeros(nt,nu);
for k = 1:nt
    Kt = eval2obj(Ktvmat,T(k));
    du = -Kt.MatArray*dx(k,:)';
    if nu == 1
        U(k,:) = Ubar(k,:) + [du,0,0,0];
    elseif nu == 4
        U(k,:) = Ubar(k,:) + du';
    end
end