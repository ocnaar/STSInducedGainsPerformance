%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez-Aroche                                                  %
% Berkeley Center for Control and Identification                          %
% Spring 2017                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Perform simulations, and show plots for the state, input, and output    %
% trajectories for two sit-to-stand (STS) movements, STS 1 and STS 2, as  %
% presented in:                                                           %
% O. Narvaez-Aroche, A. Packard, M. Arcak, “Finite time robust control of %
% the Sit-To-Stand movement for powered lower limb orthoses”, Proceedings %
% of The American Control Conference, June 2018.                          %
% http://dx.doi.org/10.23919/ACC.2018.8431465                             %
%                                                                         %
% Input                                                                   %
%                                                                         %
% Q1: 6 by 6 matrix weight to compute LQR gain for STS 1.                 %
% R1: nu by nu matrix weight to compute LQR gain for STS 1.               %
% S1: 6 by 6 matrix weight to compute LQR gain for STS 1.                 %
% Q2: 6 by 6 matrix weight to compute LQR gain for STS 2.                 %
% R2: nu by nu matrix weight to compute LQR gain for STS 2.               %
% S2: 6 by 6 matrix weight to compute LQR gain for STS 2.                 %
% Atv1: 6 by 6 LTV object with system matrices from the Jacobian          %
% 	linearization of the three-link robot model about the reference       %
% 	trajectory for STS 1.                                                 %
% B2tv1: 6 by nu LTV object with input matrices from the Jacobian         %
% 	linearization of the three-link robot model about the reference       %
% 	trajectory for STS 1.                                                 %
% Atv2: 6 by 6 LTV object with system matrices from the Jacobian          %
% 	linearization of the three-link robot model about the reference       %
% 	trajectory for STS 2.                                                 %
% B2tv2: 6 by nu LTV object with input matrices from the Jacobian         %
% 	linearization of the three-link robot model about the reference       %
% 	trajectory for STS 2.                                                 %
% tgrid: m by 1 time array for samples in reference trajectories.         %
% xbar1: m by 6 array with reference state trajectory for STS 1.          %
% ubar1: m by nu array with reference input trajectory for STS 1.         %
% xbar2: m by 6 array with reference state trajectory for STS 2.          %
% ubar2: m by nu array with reference input trajectory for STS 2.         %
% UQpar: 1 by n cell array with structures containing values for the      %
% 	parameters of the three-link robot. The dynamics of the system in     %
% 	closed-loop with the LQR controllers for the STS movements are solved %
% 	for each of these n values. The reference trajectories in the plots   %
% 	are taken from the simulation for UQpar{end}.                         %
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
% Output                                                                  %
%                                                                         %
% Show plots on screen.                                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plotACC2018(Q1,R1,S1,Q2,R2,S2,Atv1,B2tv1,Atv2,B2tv2,tgrid,xbar1,ubar1,xbar2,ubar2,UQpar)

% Number of states nx, and control inputs nu.
[nx,nu] = size(B2tv1);

% Weight matrices. 
Q = {Q1,Q2};
R = {R1,R2};
S = {S1,S2};

% Number of STS movements for simulation.
ne = numel(Q);

% System matrices.
A = {Atv1,Atv2};

% Input matrices.
B = {B2tv1,B2tv2};

% Reference state trajectories.
xbar = {xbar1,xbar2};

% Reference input trajectories.
ubar = {ubar1,ubar2};

% Number of parameter values for the three-link robot used for simulations.
nuq = numel(UQpar);

% Cell arrays for storing data from simulations.
TLQR = cell(ne,nuq); % Time. 
XLQR = cell(ne,nuq); % State. 
ULQR = cell(ne,nuq); % Input.
Z = cell(ne,nuq);    % Output. 

% Bounds for state during ode45 simulations. 
xmin = [70;-130;-10;-25;-20;-80]*pi/180;
xmax = [130;10;150;15;70;30]*pi/180;

% Set options for ode45 events function.
options = odeset('Events',@(t,x)TVLQRSTS3LinkEvents(t,x,xmin,xmax));

% Simulate closed-loop systems for each controller. 
tStart = tic;
for i=1:ne
    fprintf('\nSimulating STS %d for %d parameter values...\n',i,nuq);
    
    % Compute LTV LQR gain for the STS movement. 
    K = lqrfh(A{i},B{i},Q{i},R{i},S{i},[tgrid(1) tgrid(end)]);
    
    % Simulate three-link robot system under LQR control, using parallel computation tools.
    tic
    parfor j=1:nuq
		% State trajectory.
        [TLQR{i,j},XLQR{i,j}] = ode45(@(t,x) TVLQRSTS3Link(t,x,tgrid,xbar{i},ubar{i},K,UQpar{j}),[tgrid(1) tgrid(end)],xbar{i}(1,:)',options);
		% Input trajectory.
		ULQR{i,j} = TVLQRThreeLinkInputs(TLQR{i,j},XLQR{i,j},K,tgrid,xbar{i},ubar{i});
		% Output trajectory.
		Z{i,j} = xpar2CoMpv(XLQR{i,j}',[UQpar{j}.m1; UQpar{j}.m2; UQpar{j}.m3; UQpar{j}.I1; UQpar{j}.I2; UQpar{j}.I3; UQpar{j}.l1; UQpar{j}.l2; UQpar{j}.l3; UQpar{j}.lc1;UQpar{j}.lc2; UQpar{j}.lc3]);
    end
    toc
end

% Offset for figure numbering. 
noff = 0;

% State Trajectories.
 
% Initial plots for adding legends.
for i=1:nx
    figure(noff+i)
    hold on
    plot(TLQR{1,nuq},XLQR{1,nuq}(:,i)*180/pi,'b-',TLQR{2,nuq},XLQR{2,nuq}(:,i)*180/pi,'r-','LineWidth',2)
    legend('STS 1','STS 2','Location','Best')
end

% Plot STS 1 state trajectories for sampled parameters.
for j = 1:nuq-1
    for i=1:nx
        figure(noff+i)
        plot(TLQR{1,j},XLQR{1,j}(:,i)*180/pi,'b-','LineWidth',2)
    end
end

% Plot STS 2 state trajectories for sampled parameters.
for j = 1:nuq-1
    for i=1:nx
        figure(noff+i)
        plot(TLQR{2,j},XLQR{2,j}(:,i)*180/pi,'r-','LineWidth',2)
    end
end

% Plot reference state trajectories.
for i=1:nx
    figure(noff+i)
    plot(TLQR{1,nuq},XLQR{1,nuq}(:,i)*180/pi,'w--',TLQR{2,nuq},XLQR{2,nuq}(:,i)*180/pi,'w--','LineWidth',3);
    xlim([tgrid(1) tgrid(end)]);
    grid
	% Add labels. 
    xlabel('$t\:[s]$','Interpreter','Latex')
    if i <= (nx/2)
        lbl = ['$\theta_{' num2str(i) '}(t)\:[^\circ]$'];
        ylabel(lbl,'Interpreter','Latex')
    else
        lbl = ['$\dot{\theta}_{' num2str(i-nx/2) '}(t)\:[^\circ/s]$'];
        ylabel(lbl,'Interpreter','Latex')
    end
end

% Input trajectories

% Initial plots for adding legends.
for i=1:nu
    figure(noff+nx+i)
    hold on
    plot(TLQR{1,nuq},ULQR{1,nuq}(:,i),'b-',TLQR{2,nuq},ULQR{2,nuq}(:,i),'r-','LineWidth',2);
    legend('STS 1','STS 2','Location','Best')
end

% Plot STS 1 input trajectories for sampled parameters.
for j = 1:nuq-1
    for i = 1:nu
        figure(noff+nx+i)
        plot(TLQR{1,j},ULQR{1,j}(:,i),'b-','LineWidth',2)
    end
end

% Plot STS 2 input trajectories for sampled parameters.
for j = 1:nuq-1
    for i=1:nu
        figure(noff+nx+i)
        plot(TLQR{2,j},ULQR{2,j}(:,i),'r-','LineWidth',2)
    end
end

% Plot reference input trajectories.
for i=1:nu
    figure(noff+nx+i)
    plot(TLQR{1,nuq},ULQR{1,nuq}(:,i),'w--',TLQR{2,nuq},ULQR{2,nuq}(:,i),'w--','LineWidth',3)
    xlim([tgrid(1), tgrid(end)])
    grid
	% Add labels.
    xlabel('$t\:[s]$','Interpreter','Latex');
    if i==1
        ylabel('$\tau_{1}(t)\:[N.m]$','Interpreter','Latex');
    elseif i==2
        ylabel('$\tau_{2}(t)\:[N.m]$','Interpreter','Latex');
    elseif i==3
        ylabel('$F_x(t)\:[N]$','Interpreter','Latex');
    else
        ylabel('$F_y(t)\:[N]$','Interpreter','Latex');
    end
end

% Output trajectories.

% Initial plots for adding legends.
j = nuq;
for k = 1:nx
    figure(noff+nx+nu+k)
    hold on
    if k == 1
        % x coordinate for the position of the Center of Mass (CoM) of the robot.
        plot(TLQR{1,j},Z{1,j}(2,:),'b-',TLQR{2,j},Z{2,j}(2,:),'r-','LineWidth',2);
    elseif k == 2
        % y coordinate for the position of the CoM of the robot.
        plot(TLQR{1,j},Z{1,j}(3,:),'b-',TLQR{2,j},Z{2,j}(3,:),'r-','LineWidth',2);
    elseif k == 3
        % Position of the CoM of the robot.
        plot(Z{1,j}(2,:),Z{1,j}(3,:),'b-',Z{2,j}(2,:),Z{2,j}(3,:),'r-','LineWidth',2);
    elseif k == 4
        % x coordinate for the velocity of the CoM of the robot. 
        plot(TLQR{1,j},Z{1,j}(5,:),'b-',TLQR{2,j},Z{2,j}(5,:),'r-','LineWidth',2);
    elseif k == 5
        % y coordinate for the velocity of the CoM of the robot.
        plot(TLQR{1,j},Z{1,j}(6,:),'b-',TLQR{2,j},Z{2,j}(6,:),'r-','LineWidth',2);
    else
        % Velocity of the CoM of the robot.
        plot(Z{1,j}(5,:),Z{1,j}(6,:),'b-',Z{2,j}(5,:),Z{2,j}(6,:),'r-','LineWidth',2);
    end
    legend('STS 1','STS 2','Location','Best')
end

% Plot STS 1 output trajectories for sampled parameters.
for j = 1:nuq
    for k = 1:nx
        figure(noff+nx+nu+k)
        if k == 1
            % x coordinate for the position of the CoM of the robot.
            plot(TLQR{1,j},Z{1,j}(2,:),'b-','LineWidth',2);
        elseif k == 2
            % y coordinate for the position of the CoM of the robot.
            plot(TLQR{1,j},Z{1,j}(3,:),'b-','LineWidth',2);
        elseif k == 3
            % Position of the CoM of the robot.
            plot(Z{1,j}(2,:),Z{1,j}(3,:),'b-','LineWidth',2);
        elseif k == 4
            % x coordinate for the velocity of the CoM of the robot.
            plot(TLQR{1,j},Z{1,j}(5,:),'b-','LineWidth',2);
        elseif k == 5
            % y coordinate for the velocity of the CoM of the robot.
            plot(TLQR{1,j},Z{1,j}(6,:),'b-','LineWidth',2);
        else
            % Velocity of the CoM of the robot.
            plot(Z{1,j}(5,:),Z{1,j}(6,:),'b-','LineWidth',2);
        end
    end
end

% Plot STS 2 output trajectories for sampled parameters.
for j = 1:nuq
    for k = 1:nx
        figure(noff+nx+nu+k)
        if k == 1
            % x coordinate for the position of the CoM of the robot.
            plot(TLQR{2,j},Z{2,j}(2,:),'r-','LineWidth',2);
        elseif k == 2
            % y coordinate for the position of the CoM of the robot.
            plot(TLQR{2,j},Z{2,j}(3,:),'r-','LineWidth',2);
        elseif k == 3
            % Position of the CoM of the robot.
            plot(Z{2,j}(2,:),Z{2,j}(3,:),'r-','LineWidth',2);
        elseif k == 4
            % x coordinate for the velocity of the CoM of the robot.
            plot(TLQR{2,j},Z{2,j}(5,:),'r-','LineWidth',2);
        elseif k == 5
            % y coordinate for the velocity of the CoM of the robot.
            plot(TLQR{2,j},Z{2,j}(6,:),'r-','LineWidth',2);
        else
            % Velocity of the CoM of the robot.
            plot(Z{2,j}(5,:),Z{2,j}(6,:),'r-','LineWidth',2);
        end
    end
end

% Plot output reference trajectories.
for k = 1:nx
    figure(noff+nx+nu+k)
    if k == 1
        % x coordinate for the position of the CoM of the robot.
        plot(TLQR{1,nuq},Z{1,nuq}(2,:),'w--',TLQR{2,nuq},Z{2,nuq}(2,:),'w--','LineWidth',3);
        grid
        xlabel('$t\:[s]$','Interpreter','Latex')
        ylabel('$x_{CoM}\:[m]$','Interpreter','Latex')
    elseif k == 2
        % y coordinate for the position of the CoM of the robot.
        plot(TLQR{1,nuq},Z{1,nuq}(3,:),'w--',TLQR{2,nuq},Z{2,nuq}(3,:),'w--','LineWidth',3);
        grid
        xlabel('$t\:[s]$','Interpreter','Latex')
        ylabel('$y_{CoM}\:[m]$','Interpreter','Latex')
    elseif k == 3
        % Position of the CoM of the robot.
        plot(Z{1,nuq}(2,:),Z{1,nuq}(3,:),'w--',Z{2,nuq}(2,:),Z{2,nuq}(3,:),'w--','LineWidth',3);
        grid
        xlabel('$x_{CoM}\:[m]$','Interpreter','Latex')
        ylabel('$y_{CoM}\:[m]$','Interpreter','Latex')
    elseif k == 4
        % x coordinate for the velocity of the CoM of the robot.
        plot(TLQR{1,nuq},Z{1,nuq}(5,:),'w--',TLQR{2,nuq},Z{2,nuq}(5,:),'w--','LineWidth',3);
        grid
        xlabel('$t\:[s]$','Interpreter','Latex')
        ylabel('$\dot{x}_{CoM}\:[m/s]$','Interpreter','Latex')
    elseif k == 5
        % y coordinate for the velocity of the CoM of the robot.
        plot(TLQR{1,nuq},Z{1,nuq}(6,:),'w--',TLQR{2,nuq},Z{2,nuq}(6,:),'w--','LineWidth',3);
        grid
        xlabel('$t\:[s]$','Interpreter','Latex')
        ylabel('$\dot{y}_{CoM}\:[m/s]$','Interpreter','Latex')
    else
        % Velocity of the CoM of the robot.
        plot(Z{1,nuq}(5,:),Z{1,nuq}(6,:),'w--',Z{2,nuq}(5,:),Z{2,nuq}(6,:),'w--','LineWidth',3);
        grid
        xlabel('$\dot{x}_{CoM}\:[m/s]$','Interpreter','Latex')
        ylabel('$\dot{y}_{CoM}\:[m/s]$','Interpreter','Latex')
    end
end

% Show total running time on screen.
tTotal = toc(tStart);
fprintf('\nFinished plot generation.\nTotal running time was %d[s] (%d[h]).\n',tTotal,tTotal/3600);