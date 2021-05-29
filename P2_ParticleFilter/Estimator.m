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
N_particles = 1000; % obviously, you will need more particles than 10.
%% Mode 1: Initialization
if (km == 0)
    % Do the initialization of your estimator here!
    mu_0 = [[estConst.pA(1); estConst.pB(1)]; [estConst.pA(2); estConst.pB(2)]; 0; 0];
    var_0 = [estConst.d; estConst.d; estConst.phi_0; estConst.l];
    x_0 = sampleInitial(mu_0, var_0, N_particles);
    postParticles.x_r =  x_0(:, 1)';% 1xN_particles matrix
    postParticles.y_r = x_0(:, 2)';% 1xN_particlses matrix
    postParticles.phi = x_0(:, 3)';% 1xN_particles matrix
    postParticles.kappa = x_0(:, 4)';% 1xN_parunrticles matrix    
    
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If km > 0, we perform a regular update of the estimator.
% Implement your estimator here!

procNoise = zeros(2, N_particles); % Process noise, [2xN_particles]
priors = zeros(4, N_particles);
points = zeros(size(estConst.contour, 1), size(estConst.contour, 2), N_particles); % contour map for each kappa, [10x2xN_particles]
p_z = zeros(1, N_particles); %Probability distribution of the observations based on the prior

%Propagate particles through dynamics
for i = 1:N_particles
    procNoise(:, i) = getProcessNoise([estConst.sigma_f;
                         estConst.sigma_phi]); 

    % Prior Update:
    priors(:, i) = q([prevPostParticles.x_r(i);
                prevPostParticles.y_r(i);
                prevPostParticles.phi(i);
                prevPostParticles.kappa(i)], act, procNoise(:, i));
    
    %Kappa update on contour
    points(:, :, i) = estConst.contour;
    points(8, 1, i) =  priors(4, i);
    points(9, 1, i)  = priors(4, i); 
    
    %Find intersection 
    interPoint = findIntersectionPoints(priors(3, i), points(:, :, i), priors(1:2, i)); %[x_c; y_c] - coordinates of the laser intersection of the contour
    df = priors(1:2, i) - interPoint; 
    w_k = sens - sqrt(sum(df .* df));
    p_z(i) = pdf_wk(w_k, estConst.epsilon);
end

alpha = 1/sum(p_z);
beta = alpha * p_z;

postParticles.x_r = zeros(1, N_particles);
postParticles.y_r = zeros(1, N_particles);
postParticles.phi = zeros(1, N_particles);
postParticles.kappa = zeros(1, N_particles);

%Measurement update + resample
idx = getResample(beta);
postParticles.x_r = priors(1, idx);
postParticles.y_r = priors(2, idx);
postParticles.phi = priors(3, idx);
postParticles.kappa = priors(4, idx);

end % end estimator

%Disclaimer - probably can be cleaned up quite a bit, e.g. using 
% getSample more generally to sample a pdf f with a cdf F

function xkp1 = q(xp, u, v)
    xkp1 = [xp(1) + (u(1) + v(1)) .* cos(xp(3));
            xp(2) + (u(1) + v(1)) .* sin(xp(3));
            xp(3) + u(2) + v(2);
            xp(4)];
end

function zk = w(position, interPoint, wk)
    diff = (position-interPoint);
    tmp = diff .* diff;
    zk = sqrt(tmp(1) + tmp(2)) + wk;
end

function interPoint = findIntersectionPoints(phi, points, position)
    interPoint = zeros(2, numel(phi));
    for j = 1:size(points, 1)
        if j == size(points, 1)
            p2 = 1;
        else
            p2 = j+1;
        end
        p1 = j;
        phi1 = findAngle(position, points(p1, :));
        phi2 = findAngle(position, points(p2, :));
        if ((phi1 <= phi) && (phi < phi2)) || ((phi1 < phi) && (phi <= phi2))
            line = getLine(points(p1, :), points(p2, :));
            line_c = [tan(phi), position(2) - tan(phi)*position(1)];
            if ((line(1) == 0) && (line(2) == 0)) % Singular --> vertical contour
                y = line_c(1)*points(p1, 1) + line_c(2);
                out = [points(p1, 1); y];
            else
                out = findIntersection(line, line_c);
            end
            interPoint = [out(1); out(2)];
            break;
        end
    end
end

function phi = findAngle(p1, p2)
    phi = atan2(p2(2) - p1(2), p2(1) - p1(1));
end

function res = getLine(p1, p2)
    if det([p1(1) 1; p2(1) 1]) ~= 0
        res = [p1(1) 1; p2(1) 1] \ ([p1(2); p2(2)]); 
    else
        res = [0; 0];
    end
end

function res = findIntersection(line1, line2)
    x = (line2(2) - line1(2))/(line1(1) - line2(1));
    y = line1(1) * x + line1(2);
    res = [x; y];
end

function procNoise = getProcessNoise(sigma)
    procNoise = rand(2, 1);
    procNoise(1) = (sigma(1)/2)*(2*procNoise(1) - 1);
    procNoise(2) = (sigma(2)/2)*(2*procNoise(2) - 1);
end

function idx = getResample(beta) 
    idx = zeros(1, numel(beta));
    u = rand(1, numel(beta));
    for j = 1:numel(beta)
        k = 0;
        sum = 0;
        while sum < u(j)
            k = k+1;
            sum = sum + beta(k);
        end
        idx(j) = k;
    end
end

function pwk = pdf_wk(x, epsilon)
    if x<0
        x = -x;
    end
    if x<=2*epsilon
        pwk = (2/(5*epsilon)) * (1 - (1/(2*epsilon)) * x);
    elseif (x>2*epsilon) && (x<= 2.5*epsilon)
        pwk = (2/(5*epsilon*epsilon)) * (x - 2*epsilon);
    elseif (x>2.5*epsilon) && (x<=3*epsilon)
        pwk = 1/(5*epsilon) - (2/(5*epsilon*epsilon)) * (x - 2.5*epsilon);
    else
        pwk = 0;
    end
end

function x_bar = sampleInitial(mu, var, N) %Sample uniform around (p0-x_bar) <= d 
    %Stupid implementation for now, fix later
    stateSize = numel(var);
    u_bar = rand(N, stateSize); %Mutual independence of the initial conditions
    x_bar = zeros(N, stateSize);
    for k = 1:N
        decision = rand;
        for j = 1:stateSize
            if j==1
                if decision < 0.5
                    x_bar(k, j) = var(j)*(2*u_bar(k, j) - 1) + mu(j);
                else
                    x_bar(k, j) = var(j)*(2*u_bar(k, j) - 1) + mu(j+1);
                end
            elseif j==2
                bound = 2*u_bar(k, j-1) - 1;
                if decision < 0.5
                    x_bar(k, j) = sqrt(var(j)*var(j)*(1 - bound*bound))*(2*u_bar(k, j) - 1) + mu(j+1);
                else
                    x_bar(k, j) = sqrt(var(j)*var(j)*(1 - bound*bound))*(2*u_bar(k, j) - 1) + mu(j+2);
                end
            else
                x_bar(k, j) = var(j)*(2*u_bar(k, j) - 1) + mu(j+2);
            end
        end
    end
end