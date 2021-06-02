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
N_particles = 3000; 
% Roughening constant
K = 0.001;
% Resampling when probabilities collapse
reSampSize = 10*N_particles;

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
    postParticles.beta = 0;
    postParticles.prevPriors = zeros(4, N_particles);
    
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
w_k = zeros(1, N_particles);

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
    d = get_distance(priors(:,i), points(:,:,i)) ;
    w_k(i)=sens - d; 
    p_z(i) = pdf_wk(w_k(i), estConst.epsilon);
end

if sum(p_z) == 0
    p_z_rr = zeros(1, reSampSize);
    while (sum(p_z_rr) == 0) 
        mu = [mean(prevPostParticles.x_r); mean(prevPostParticles.y_r); mean(prevPostParticles.phi); mean(prevPostParticles.kappa)];
        bound = [0.2; 0.2; 2*pi; 0.2];
        priors_rr = uniformResample(mu, bound, reSampSize, estConst.contour);
        w_rr = zeros(1, reSampSize);
        points_rr = zeros(size(estConst.contour, 1), size(estConst.contour, 2), reSampSize);
        for i = 1:reSampSize
            points_rr(:, :, i) = estConst.contour;
            points_rr(8, 1, i) =  priors_rr(4, i);
            points_rr(9, 1, i)  = priors_rr(4, i); 
            d_rr = get_distance(priors_rr(:,i), points_rr(:,:,i)) ;
            w_rr(i)=sens - d_rr; 
            p_z_rr(i) = pdf_wk(w_rr(i), estConst.epsilon);
        end
        [~, idx_rr] = maxk(p_z_rr, N_particles);
    end
    priors = priors_rr(:, idx_rr);
    p_z = p_z_rr(idx_rr);
end
alpha = 1/sum(p_z);
beta = alpha * p_z;

postParticles.x_r = zeros(1, N_particles);
postParticles.y_r = zeros(1, N_particles);
postParticles.phi = zeros(1, N_particles);
postParticles.kappa = zeros(1, N_particles);

rough_sigmas = getRougheningSigma(K, priors, N_particles);

%Measurement update + resample
idx = getResample(beta);

postParticles.x_r = priors(1, idx) + getNormalSample(0, rough_sigmas(1), [1, N_particles]);
postParticles.y_r = priors(2, idx) + getNormalSample(0, rough_sigmas(2), [1, N_particles]);
postParticles.phi = priors(3, idx) + getNormalSample(0, rough_sigmas(3), [1, N_particles]);
postParticles.kappa = priors(4, idx) + getNormalSample(0, rough_sigmas(4), [1, N_particles]);
end % end estimator

%Disclaimer - probably can be cleaned up quite a bit, e.g. using 
% getSample more generally to sample a pdf f with a cdf F

function samps = uniformResample(mu, bounds, N_particles, contour)
    samps = zeros(4, N_particles);
    for i = 1:N_particles
        if i == 1
            samps(:, i) = mu + (rand(4, 1)*2 - 1) .* bounds;
        else
            while ~inpolygon(samps(1, i), samps(2, i), contour(:, 1), contour(:, 2))
                samps(:, i) = mu + (rand(4, 1)*2 - 1) .* bounds;
            end
        end
    end
end


function xkp1 = q(xp, u, v)
    xkp1 = [xp(1) + (u(1) + v(1)) .* cos(xp(3));
            xp(2) + (u(1) + v(1)) .* sin(xp(3));
            xp(3) + u(2) + v(2);
            xp(4)];
end

function sigmas = getRougheningSigma(K, priors, N)
    sigmas = K * (N^(-1/4)) * (max(priors, [], 2) - min(priors, [], 2));
end

function samp = getNormalSample(mu, sigma, size)
    samp = mu + sqrt(2)*sigma*erfinv(2*rand(size(1), size(2))-1);
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

function pwk = pdf_wk(x_in, epsilon)
    if x_in<0
        x_in = -x_in;
    end
    if x_in<=2*epsilon
        pdf = @(x) ((2/(5*epsilon)) * (1 - ((1/(2*epsilon))* x)));
    elseif (x_in>2*epsilon) && (x_in<= 2.5*epsilon)
        pdf = @(x) (((1/(5*epsilon))/(0.5*epsilon))*(x - 2*epsilon));
    elseif (x_in>2.5*epsilon) && (x_in<=3*epsilon)
        pdf = @(x) (-((1/(5*epsilon))/(0.5*epsilon))*(x - 3*epsilon));
    else
        pwk = 0;
        return;
    end
    delta = epsilon/1000; %Completely random 
    pwk = 0.5*delta*(pdf(x_in+0.5*delta) + pdf(x_in-0.5*delta)); % Trapezoid integration lol long time no see
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

function [d]= get_distance(priors,contour) % NOTE: to vectorize!
  % solve system of equations with two unknows: fraction of length (point 1
  % to intersection point) and distance from robot to intersection point. 
  % A[ fraction;distance] = b with A and b known 
  
    fraction = zeros(1,10);
    distance  = zeros(1,10);
    % read corner coordinate
    for i = 1:10
        x1 = contour( i ,1);
        y1 = contour( i ,2);
        
        if i==10
            x2 = contour(1 ,1);
            y2 = contour(1,2);
        else
            x2 = contour( i+1 ,1);
            y2 = contour( i+1,2);
        end
        A = [ x1-x2, -cos(priors(3));
              y1-y2 -sin(priors(3))];
        b = [priors(1)-x2; priors(2)-y2];
        
%         solution = inv(A)\b;
%         fraction(i) = solution(1);
%         distance(i) = solution (2);
        % >> Too slow, try to make it faster by expanding manually
%         
        determinant = sin(priors(3) )*(x2-x1)+cos(priors(3) )*(y1-y2);
        fraction(i) = (sin(priors(3) )*(x2-priors(1) )+cos(priors(3) )*(priors(2) -y2))/determinant;
        distance(i)  = ( (y2-y1)*(priors(1) -x2)+(x1-x2)*(priors(2) -y2) )/determinant;
    end
    % determine the wall intersected
    % fraction of the segment
    mask = fraction >= 0 & fraction <= 1 & distance >= 0;
    B = distance.*mask;
    d = min(B(B~=0));
    if isempty(d)
        d = inf;
    end
end