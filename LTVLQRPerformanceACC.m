%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez-Aroche                                                  %
% Berkeley Center for Control and Identification                          %
% Summer 2017                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Compute a performance metric for finite time horizon LQR controllers    %
% used to track a sit-to-stand (STS) movement. The performance metric is  %
% based on the calculation of L2-to-Euclidean gains, it is defined in:    %
% O. Narvaez-Aroche, A. Packard, M. Arcak, “Finite time robust control of %
% the Sit-To-Stand movement for powered lower limb orthoses”, Proceedings %
% of The American Control Conference, June 2018.                          %
% http://dx.doi.org/10.23919/ACC.2018.8431465                             %
%                                                                         %
% Input                                                                   %
%                                                                         %
% LHW: ne by nw array with normalized values for the diagonal entries of  %
% 	the weight matrices (Q, R, S) used in the design of the LQR gains.    % 
% wmin: nw by 1 vector with lower bounds for the entries of the weight    % 
% 	matrices.                                                             %  
% wmax: nw by 1 vector with upper bounds for the entries of the weight    % 
% 	matrices.                                                             % 
% A: nx by nx by nt LTV object with state matrices from the Jacobian      %
% 	linearization of the system about the reference trajectory for the    %
% 	STS movement.                                                         %
% B1: nx by nz by nt LTV object with parameter matrices.                  %
% B2: nx by nu by nt LTV object with input matrices.                      %
% C: nz by nx by nt LTV object with output matrices.                      %
% D1: nz by np by nt LTV object with parameter feedthrough matrices.      %
% tgrid: time array with nt samples in [s].                               %
% tm: intermediate time of the STS movement in [s] where the first        %
% 	finite-horizon L2-to-Euclidean gain calculation takes place.          %
% wtf: value between 0 and 1 for weighting the intermediate, and final    %
% 	induced gains in the performance metric.                              %    
% tau: time constant for the low-pass filter applied to the time-varying  %
% 	signals that model the parameter uncertainty of the system.           %
% We: weight matrix for the output of the extended LTV system.            %
% parmin: np by 1 vector with lower bound for parameter uncertainty.      %
% parmax: np by 1 vector with upper bound for parameter uncertainty.      %
%                                                                         %
% Output                                                                  %
%                                                                         %
% J: ne by 1 array with the values for the performance metric achieved by %
% 	the LQR gains obtained from the weights in LHW array.                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = LTVLQRPerformanceACC(LHW,wmin,wmax,A,B1,B2,C,D1,tgrid,tm,wtf,tau,We,parmin,parmax)

% Number of weight matrices triplets (Q, R, S) ne, and number of matrix entries used for their definition nw.
[ne,nw] = size(LHW);

% Number of states nx, and control inputs nc. 
[nx,nc] = size(B2);

fprintf('\nComputing performance metrics for %d finite time horizon LQR controllers for tracking a sit-to-stand movement...\n',ne);

% Array for storing the diagonal elements of the weight matrices triplets.
weights = zeros(nw,ne);

% Array for storing robust performance metric. 
J = zeros(ne,1);

% State space realization of the open-loop dynamics for the three-link robot model.
G = sstv(A,[B1 B2],eye(nx),0);

% State space realization of the static mapping from the space of theta into the space of z.
C = sstv(0,zeros(1,6),zeros(6,1),C);
D1 = sstv(0,zeros(1,12),zeros(6,1),D1);

% Low-pass filters for time-varying signals modeling the parameter uncertainty of the system.
filtsys = tf(1,[tau 1]);
Wd = filtsys*diag((parmax-parmin)/2);

% Cell for storing the extended LTV systems used for computing the induced gains in the performance metric. 
H = cell(1,ne);

% Compute diagonal entries for weight matrices triplets (Q, R, S).
for i=1:ne
    weights(:,i) = wmin + (wmax-wmin).*LHW(i,:)';
end

% Loop for computing robust performance metrics.
tStart = tic;
for i=1:ne
    % The weight matrices are given by:
    % Q = diag(weights(1:nx));
    % R = diag(weights(nx+1:nx+nc));
    % S = diag(weights(nx+nc+1:end));
    
    % The finite time horizon LQR gain is obtained from:
    % K = lqrfh(A,B2,diag(weights(1:nx,i)),diag(weights(nx+1:nx+nc,i)),diag(weights(nx+nc+1:end,i)),[tgrid(1) tgrid(end)]);

    % The state space realization of the extended LTV system is:
    % T = feedback(G,K,13:16,1:6);
    % R = parallel(D1,series(T,C),1:12,1:12,1:6,1:6);
    % H = eval2obj(series(Wd,R,1:12,1:12),tgrid);
	
	% The extraction of the weight matrices, calculation of the LQR gain, and the state space realization of the extended LTV system are all conducted in the following line to reduce running time, and memory usage.  
    H{i} = eval2obj(We*series(Wd,parallel(D1,series(feedback(G,lqrfh(A,B2,diag(weights(1:nx,i)),diag(weights(nx+1:nx+nc,i)),diag(weights(nx+nc+1:end,i)),[tgrid(1) tgrid(end)]),13:16,1:6),C),1:12,1:12,1:6,1:6),1:12,1:12),tgrid);
    
    % Compute robust performance metric for LQR controller from induced gains. 
    J(i) = (1-wtf)*L2to2normRicc(H{i},'horizon',[tgrid(1),tm]) + wtf*L2to2normRicc(H{i},'horizon',[tgrid(1),tgrid(end)]);
end

% Show total running time on screen.
tTotal = toc(tStart);
fprintf('\nFinished computation of the robust performance metrics.\nTotal running time was %d[s] (%d[h]).\n',tTotal,tTotal/3600);