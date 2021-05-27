function [postParticles] = Estimator(prevPostParticles, sens, act, estConst, km)
% The estimator function. The function will be called in two different
% modes: If km==0, the estimator is initialized. If km > 0, the
% estimator does an iteration for a single sample time interval using the 
% previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurement and control inputs.
%
% Inputs:
%   prevPostParticles   previous posterior particles at time step k-1
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%                           
%   sens                Sensor measurement z(k), scalar
%
%   act                 Control inputs u(k-1), [1x2]-vector
%                       act(1): u_f, forward control input
%                       act(2): u_phi, angular control input
%
%   estConst            estimator constants (as in EstimatorConst.m)
%
%   km                  time index k, scalar
%                       corresponds to continous time t = k*Ts
%                       If km==0 initialization, otherwise estimator
%                       iteration step.
%
% Outputs:
%   postParticles       Posterior particles at time step k
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch

% Set number of particles:
N_particles = 10; % obviously, you will need more particles than 10.

%% Mode 1: Initialization
if (km == 0)
    % Do the initialization of your estimator here!
    mu_0 = [[estConst.pA(1); estConst.pB(1)]; [estConst.pA(2); estConst.pB(2)]; 0; 0];
    var_0 = [estConst.d; estConst.d; estConst.phi_0; estConst.l];
    x_0 = sampleInitial(mu_0, var_0, N_particles);
    postParticles.x_r =  x_0(:, 1)';% 1xN_particles matrix
    postParticles.y_r = x_0(:, 2)';% 1xN_particles matrix
    postParticles.phi = x_0(:, 3)';% 1xN_particles matrix
    postParticles.kappa = x_0(:, 4)'% 1xN_particles matrix
    
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If km > 0, we perform a regular update of the estimator.

% Implement your estimator here!


% Prior Update:


% Posterior Update:

postParticles.x_r = ...
postParticles.y_r = ...
postParticles.phi = ...
postParticles.kappa = ...

end % end estimator


function x_bar = sampleInitial(mu, var, N) %Sample uniform around (p0-x_bar) <= d 
    %Stupid implementation for now, fix later
    stateSize = numel(var);
    u_bar = rand(N, stateSize); %Mutual independence of the initial conditions
    x_bar = zeros(N, stateSize);
    for i = 1:N
        decision = rand;
        for j = 1:stateSize
            if j==1
                if decision < 0.5
                    x_bar(i, j) = var(j)*(2*u_bar(i, j) - 1) + mu(j);
                else
                    x_bar(i, j) = var(j)*(2*u_bar(i, j) - 1) + mu(j+1);
                end
            elseif j==2
                bound = 2*u_bar(i, j-1) - 1;
                if decision < 0.5
                    x_bar(i, j) = sqrt(var(j)*var(j)*(1 - bound*bound))*(2*u_bar(i, j) - 1) + mu(j+1);
                else
                    x_bar(i, j) = sqrt(var(j)*var(j)*(1 - bound*bound))*(2*u_bar(i, j) - 1) + mu(j+2);
                end
            else
                x_bar(i, j) = var(j)*(2*u_bar(i, j) - 1) + mu(j+2);
            end
        end
    end
end