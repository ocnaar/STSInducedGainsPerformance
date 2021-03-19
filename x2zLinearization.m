%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez-Aroche                                                  %
% Berkeley Center for Control and Identification                          %
% Summer 2017                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Jacobian linearization of function f about state, and parameter         %
% trajectories xbar, pbar using two-point finite difference formulas.     %
%                                                                         %
% Input                                                                   %
%                                                                         %
% f: function handle.                                                     %
% xbar: m by nx array with state trajectory.                              %
% pbar: m by np array with paramater trajectory.                          %
% eps: initial spacing for the two points in finite difference formula.   %
% err: error tolerance in succesive computations of matrix columns.       %
% maxiter: maximum number of iterations allowed for computing matrix      %
%   columns.                                                              %
%                                                                         %
% Output                                                                  %
%                                                                         %
% A: nx by nx by m array of state matrices.                               %
% B: nx by np by m array of parameter matrices.                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,B] = x2zLinearization(f,xbar,pbar,eps,err,maxiter)

% Number of sampled points from state, and parameter trajectories.
% Number of states.
[T,nx] = size(xbar);

% Number of parameters.
[~,np] = size(pbar);

% Array for storing system matrices.
A = zeros(nx,nx,T);

% Array for storing parameter matrices.
B = zeros(nx,np,T);

% Perform linearization about sampled points.
for k = 1:T
    % Compute system matrix at time k.
    f0 = f(xbar(k,:)',pbar(k,:)');
    I = eye(nx);
    for i=1:nx
        j=1;
        a = ones(nx,1);
        aprev = zeros(nx,1);
        while norm(a-aprev)>err && j<=maxiter
            delta = eps/j;
            aprev = a;
            f1 = f(xbar(k,:)'+delta*I(:,i),pbar(k,:)');
            a=(f1-f0)/delta;
            j = j+1;
        end
        A(:,i,k) = a;
    end
    % Compute parameter matrix at time k.
    I = eye(np);
    for i=1:np
        j = 1;
        b = ones(nx,1);
        bprev = zeros(nx,1);
        while norm(b-bprev)>err && j<=maxiter
            delta = eps/j;
            bprev = b;
            f1 = f(xbar(k,:)',pbar(k,:)'+delta*I(:,i));
            b = (f1-f0)/delta;
            j = j+1;
        end
        B(:,i,k) = b;
    end
end