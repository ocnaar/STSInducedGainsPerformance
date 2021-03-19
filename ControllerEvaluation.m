%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez-Aroche                                                  %
% Berkeley Center for Control and Identification                          %
% Summer 2017                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Technique for assessing the robustness against parameter uncertainty of %
% a batch of controllers designed for tracking the sit-to-stand (STS)     %
% movement of a powered lower limb orthosis (PLLO), as published in:      %
%                                                                         %
% O. Narváez-Aroche, A. Packard, M. Arcak, “Finite time robust control of % 
% the sit-to-stand movement for powered lower limb orthoses”, Proceedings %
% of The American Control Conference, June 2018.                          %
% https://doi.org/10.23919/ACC.2018.8431465                               %
%                                                                         %
% Please cite our work accordingly.                                       %
%                                                                         %
% Make sure to install the LTV toolbox at https://z.umn.edu/LTVTools      %
% and Matlab's Parallel Computing Toolbox before running this script.     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set a clean workspace.
clear
clc
close all

%% Motion planning for two sit-to-stand movements (STS1 & STS2).

% Nominal values for the parameters of the three-link robot.
par.m1 = 9.68;      % Mass of link 1 [kg].
par.m2 = 12.59;     % Mass of link 2 [kg].
par.m3 = 44.57;     % Mass of link 3 [kg].
par.I1 = 1.16456;   % Moment of inertia of link 1 about its CoM [kg.m^2].
par.I2 = 0.518821;  % Moment of inertia of link 2 about its CoM [kg.m^2].
par.I3 = 2.55731;   % Moment of inertia of link 3 about its CoM [kg.m^2].
par.l1 = 0.533;     % Length of link 1 [m].
par.l2 = 0.406;     % Length of link 2 [m].
par.l3 = 0.52;      % Length of link 3 [m].
par.lc1 = par.l1/2; % Distance from ankle joint to CoM of link 1 [m].
par.lc2 = par.l2/2; % Distance from knee joint to CoM of link 2 [m].
par.lc3 = par.l3/2; % Distance from hip joint to CoM of link 3 [m].

% Invariant values through simulations. 
par.g = 9.81;             % Acceleration of gravity [m/s^2].
par.uconfig = 3;          % Input configuration. 
par.T = 3.5;              % Simulation time.
par.W = diag([1,1,10,1]); % Input weights in control allocation problem.  

% Coefficients of cubic polynomials for defining reference trajectories in 
% the space of z.
par.P = [-2,3,0,0;-2,3,0,0;-2,3,0,0];

% Input weights in the control allocation problem.
par.wa = 1;   % Ankle torque.
par.wk = 1;   % Knee torque.
par.wh = 1;   % Hip torque.
par.ws = 1;   % Shoulder torque.
par.wFx = 10; % Horizontal force at the shoulder.
par.wFy = 1;  % Vertical force at the shoulder.

% Lower and upper bounds for the input, assuming a fully-actuated architecture.
% The order of the values is:
% [ankle torque; knee torque; hip torque; shoulder torque;... 
% horizontal force at the shoulder; vertical force at the shoulder].
par.lb = [-Inf; -Inf; -Inf; -Inf; -Inf; 0];     % Input lower bound.
par.ub = [Inf; Inf; Inf; Inf; Inf; Inf];        % Input upper bound.

% Design LQR controller for the synthetic input in the feedback 
% linearization control scheme.
A = [0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1; zeros(3,6)];
B = [zeros(3,3);eye(3)];
Q = eye(6);
R = 1/100*eye(3);
par.K = lqr(A,B,Q,R);

% Simulation time.
T = par.T;

% Configuration of the three-link robot at the end of the ascension 
% phase in the space of theta.
thf = [90; -5; 13.9353; zeros(6,1)]*pi/180;

% Configuration of the three-link robot at the end of the ascension 
% phase in the space of z.
zf = theta2z3link(thf,par);

% Configuration of the three-link robot at the beginning of the 
% ascension phase for STS 1 (dynamic strategy) in the space of theta.
th1 = [pi/2; -pi/2; pi/2; zeros(6,1)];

% Configuration of the three-link robot at the beginning of the 
% ascension phase for STS 1 (dynamic strategy) in the space of z.
z1 = theta2z3link(th1,par);

% Configuration of the three-link robot at the beginning of the 
% ascension phase for STS 2 (quasi-static strategy) in the space of theta.
th2 = [120; -120; 110.87; zeros(6,1)]*pi/180;

% Configuration of the three-link robot at the beginning of the 
% ascension phase for STS 2 (quasi-static strategy) in the space of z.
z2 = theta2z3link(th2,par);

%% Compute reference trajectories for state, and input.

% Solve ODEs for STS 1.
fprintf('\nUsing a computed torque approach for obtaining the reference trajectories of STS 1...\n')
tic
[T1nom,TH1nom] = ode45(@(t,th) FLThreeLinkUQ(t,th,z1,zf,par,par),[0 T],th1(1:6));
toc

% Solve ODEs for STS 2.
fprintf('\nUsing a computed torque approach for obtaining the reference trajectories of STS 2...\n')
tic
[T2nom,TH2nom] = ode45(@(t,th) FLThreeLinkUQ(t,th,z2,zf,par,par),[0 T],th2(1:6));
toc

% Reference angular trajectories for STS 1.
TH1ref = z2theta3link(desiredz(par.P,T1nom,T,z1(1:3),zf(1:3)),par)';

% Obtain reference input for STS 1 from computed torque approach.
U1nom = FLThreeLinkInputs(TH1nom,TH1ref,par);

% Reference angular trajectories for STS 2.
TH2ref = z2theta3link(desiredz(par.P,T2nom,T,z2(1:3),zf(1:3)),par)';

% Obtain reference input for STS 2 from computed torque approach.
U2nom = FLThreeLinkInputs(TH2nom,TH2ref,par);

%% Jacobian Linearization about reference trajectories for the state, and input, and nominal value of the parameter. 

% Sampling time for linearization.
ts = 0.01;

% Time array for linearization. 
tgrid = 0:ts:T;

% Interpolate STS 1 reference state trajectory at the values 
% in tgrid using cubic splines.
nx = size(TH1nom,2);
xbar1 = zeros(length(tgrid),nx);
for i=1:nx
    xbar1(:,i) = csapi(T1nom,TH1nom(:,i),tgrid);
end

% Interpolate STS 1 reference input trajectory at the values 
% in tgrid using cubic splines.
nu = size(U1nom,2);
ubar1 = zeros(length(tgrid),nu);
for i=1:nu
    ubar1(:,i) = csapi(T1nom,U1nom(:,i),tgrid);
end

% Interpolate STS 2 reference state trajectory at the values 
% in tgrid using cubic splines.
xbar2 = zeros(length(tgrid),nx);
for i=1:nx
    xbar2(:,i) = csapi(T2nom,TH2nom(:,i),tgrid);
end

% Interpolate STS 2 reference input trajectory at the values 
% in tgrid using cubic splines.
ubar2 = zeros(length(tgrid),nu);
for i=1:nu
    ubar2(:,i) = csapi(T2nom,U2nom(:,i),tgrid);
end

% Array with nominal values of the parameter used in linearization. 
pbar = [par.m1,par.m2,par.m3,par.I1,par.I2,par.I3,...
    par.l1,par.l2,par.l3,par.lc1,par.lc2,par.lc3];
pbar = repmat(pbar,[length(tgrid),1]);

% Initial spacing for the two points in the finite difference 
% formula used for linearization. 
epsilon = 1e-3;

% Error tolerance in succesive computations of matrix columns 
% during linearization.
err = 1e-20;

% Maximum number of iterations allowed for computing matrix 
% columns during linearization. 
maxiter = 1;

% Perform Jacobian linearization about STS 1. 
fprintf('\nJacobian linearization about STS 1 reference state trajectory...\n')
tic
[Alin1, B1lin1, B2lin1] = ParameterLinearization(@STSThreeLinkPar,xbar1,pbar,ubar1,par,epsilon,err,maxiter);
toc

% Perform Jacobian linearization about STS 2. 
fprintf('\nJacobian linearization about STS 2 reference state trajectory...\n')
tic
[Alin2, B1lin2, B2lin2] = ParameterLinearization(@STSThreeLinkPar,xbar2,pbar,ubar2,par,epsilon,err,maxiter);
toc

% Build LTV toolbox objects from Jacobian matrices for STS 1.
Atv1 = tvmat(Alin1,tgrid);
B1tv1 = tvmat(B1lin1,tgrid);
B2tv1 = tvmat(B2lin1,tgrid);

% Build LTV toolbox objects from Jacobian matrices for STS 2.
Atv2 = tvmat(Alin2,tgrid);
B1tv2 = tvmat(B1lin2,tgrid);
B2tv2 = tvmat(B2lin2,tgrid);

%% Jacobian linearization of the mapping from the space of theta to the space of z (function \zeta in ACC paper).

% For STS 1. 
[Clin1,D1lin1] = x2zLinearization(@xpar2CoMpv,xbar1,pbar,epsilon*100,err,maxiter);
Ctv1 = tvmat(Clin1,tgrid);
D1tv1 = tvmat(D1lin1,tgrid);

% For STS 2. 
[Clin2,D1lin2] = x2zLinearization(@xpar2CoMpv,xbar2,pbar,epsilon*100,err,maxiter);
Ctv2 = tvmat(Clin2,tgrid);
D1tv2 = tvmat(D1lin2,tgrid);

%% Compute performance metrics for finite time horizon LQR controllers.

% Number of control inputs.
[~,nc] = size(B2tv1);

% Number of random variables for sampling triplets of weight matrices (Q,R,S).
if nc == 1
    % For control of the hips.
    nw = 13;
elseif nc == 4
    % For control of the hips, and upper body.
    nw = 16;
end

% Triplets of weight matrices to be tested.
nwe = 1350;

% Control random number generation for consistency of experiments.
rng(2);

% Latin Hypercube Sampling (LHS).
LHw = lhsdesign(nwe,nw);

% Arrays for weight bounds.
wmin = zeros(nw,1);
wmax = zeros(nw,1);

% Lower bounds.
wmin(1:nx) = ones(nx,1);            % For diagonal Q.
wmin(nx+1:nx+nc) = 1e-2*ones(nc,1); % For diagonal R. 
wmin(nx+nc+1:nw) = ones(nx,1);      % For diagonal S. 

% Upper bounds.
wmax(1:nx) = 100*ones(nx,1);        % For diagonal Q.
wmax(nx+1:nx+nc) = ones(nc,1);      % For diagonal R. 
wmax(nx+nc+1:nw) = 100*ones(nx,1);  % For diagonal S.

% Intermediate time of the STS movements for computing L-2 to 2-Norm gain 
% (first term of performance metric J_{RP}) in [s].
tm = 2;

% Time constant for low pass filters with 50[Hz] bandwidth applied to input.
tau = 1/(100*pi);

% Weight for L-2 to 2-Norm gain computed at the end of the STS movements. 
% (second term of performance metric J_{RP}.  
wtf = 0.7;

% Weight matrix for output e.
We = diag([1, 1, 1, 10, 10, 10]);

% Lower and upper bounds for parameter uncertainty.
parmin = [par.m1-0.1; par.m2-0.1; par.m3-0.1; par.I1-0.1; par.I2-0.1; par.I3-0.1;...
          par.l1-0.01; par.l2-0.01; par.l3-0.01; par.lc1-0.01; par.lc2-0.01; par.lc3-0.01];
parmax = [par.m1+0.1; par.m2+0.1; par.m3+0.1; par.I1+0.1; par.I2+0.1; par.I3+0.1;...
          par.l1+0.01; par.l2+0.01; par.l3+0.01; par.lc1+0.01; par.lc2+0.01; par.lc3+0.01];

% Compute performance metrics for finite time horizon LQR controllers obtained from sampled triplets of weight matrices for STS 1.
fprintf('\nComputing performance metrics for STS 1...\n');
J1 = LTVLQRPerformanceACC(LHw,wmin,wmax,Atv1,B1tv1,B2tv1,Ctv1,D1tv1,tgrid,tm,wtf,tau,We,parmin,parmax);

% Compute performance metrics for finite time horizon LQR controllers obtained from sampled triplets of weight matrices for STS 2.
fprintf('\nComputing performance metrics for STS 2...\n');
J2 = LTVLQRPerformanceACC(LHw,wmin,wmax,Atv2,B1tv2,B2tv2,Ctv2,D1tv2,tgrid,tm,wtf,tau,We,parmin,parmax);

% Weights leading to the best performance metric for STS 1.
[~,norm2min1] = min(J1);
wp1 = wmin + (wmax-wmin).*LHw(norm2min1,:)';

% Weights leading to the best performance metric for STS 2.
[~,norm2min2] = min(J2);
wp2 = wmin + (wmax-wmin).*LHw(norm2min2,:)';

% Triplets of weight matrices leading to the best LQR controllers.
Q1 = diag(wp1(1:6));
Q2 = diag(wp2(1:6));
if nw == 13
    R1 = diag(wp1(7));
    S1 = diag(wp1(8:end));
	R2 = diag(wp2(7));
    S2 = diag(wp2(8:end));
elseif nw == 16
    R1 = diag(wp1(7:10));
    S1 = diag(wp1(11:end));
	R2 = diag(wp2(7:10));
    S2 = diag(wp2(11:end));
end

%% Assess Robustness to Parameter Uncertainty via Simulations

% Number of parameter values within the bounded uncertainty for performing simulations.
nuq = 200;

% Number of parameters for the three-link robot.
np = numel(parmin);

% Control random number generation for Latin Hypercube Sampling (LHS).
rng(1);

% Perform LHS.
LHuq = lhsdesign(nuq,np);

% Cell array with structures containing random parameter values within the bounded uncertainty.
UQpar = LH2par3link(LHuq,parmin,parmax);

% Add structure with nominal values for the parameter into the cell array.
UQpar{end+1} = par;

% Simulate STS movements with the LQR gains that lead to the best performance, plot results.
plotACC2018(Q1,R1,S1,Q2,R2,S2,Atv1,B2tv1,Atv2,B2tv2,tgrid,xbar1,ubar1,xbar2,ubar2,UQpar)

%% Weight matrices used for simulations in the ACC paper.
  
% For STS 1:
% Q1 = 1.0e+03*diag([3.2372, 5.5340, 6.5464, 7.9180, 4.0030, 8.5159]);
% R1 = diag([0.3659, 0.0155, 0.1433, 0.1553]);
% S1 = 1.0e+03*diag([1.0681, 5.3963, 1.3236, 9.4668, 3.9753, 5.8190]);

% For STS 2:
% Q2 = 1.0e+03*diag([3.7656, 9.5499, 2.9325, 8.3784, 9.5515, 9.2424]);
% R2 = diag([0.1119, 0.0252, 0.3600, 0.3045]);
% S2 = 1.0e+03*diag([9.5652, 0.8197, 5.3155, 5.7790, 6.0829, 8.8773]);