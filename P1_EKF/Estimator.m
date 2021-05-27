function [posEst,linVelEst,oriEst,windEst,driftEst,...
          posVar,linVelVar,oriVar,windVar,driftVar,estState] = ...
    Estimator(estState,actuate,sense,tm,estConst)
% [posEst,linVelEst,oriEst,windEst,driftEst,...
%    posVar,linVelVar,oriVar,windVar,driftVar,estState] = 
% Estimator(estState,actuate,sense,tm,estConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_t, thrust command
%                   actuate(2): u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time t_k, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConst.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   windEst         wind direction estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   windVar         variance of wind direction estimate(time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch

%% Initialization
if (tm == 0)
    % Do the initialization of your estimator here!
    
    % initial state mean
    posEst =  [0 0]; % 1x2 matrix
    linVelEst = [ 0 0 ]; % 1x2 matrix
    oriEst = 0; % 1x1 matrix
    driftEst = 0; % 1x1 matrix
    
    R0 = estConst.StartRadiusBound;
    % initial state variance
    posVar = [0.25*R0^2  0.25*R0^2]; % 1x2 matrix
    linVelVar = [0 0]; % 1x2 matrix
    
    phi_bar = estConst.RotationStartBound;
    rho_bar = estConst.WindAngleStartBound;

    oriVar = (1/3)*phi_bar^2; % 1x1 matrix
    windVar = (1/3)*rho_bar^2; % 1x1 matrix
    
    % There is no drift for initial measurement so we know drift est=0
    % without uncertainty. 
    driftVar = 0; % 1x1 matrix
    
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar linVelVar oriVar windVar driftVar]); % diagonal matrix 7x7
    % estimator state
    estState.xm = [posEst'; linVelEst'; oriEst; windEst; driftEst]; % 7x1 matrix
    % time of last update
    estState.tm = tm;
end

%% Estimator iteration.
% get time since last estimator update
dt = tm - estState.tm;
estState.tm = tm; % update measurement update time

% Easier to write variables
% constants
C_dh = estConst.dragCoefficientHydr;
C_da = estConst.dragCoefficientAir;
C_r = estConst.rudderCoefficient; 
C_w = estConst.windVel;

% Process noise variance
Q_d = estConst.DragNoise ;
Q_r = estConst.RudderNoise ;
Q_rho= estConst.WindAngleNoise ;
Q_b = estConst.GyroDriftNoise ;

% Measurement noise variance
s_a2 = estConst.DistNoiseA;
s_b2= estConst.DistNoiseB;
s_c2= estConst.DistNoiseC;
s_g2= const.GyroNoise;
s_n2 = const.CompassNoise;




%% prior update

% Prepare useful matrices
Qc=diag([estConst.DragNoise estconst.RudderNoise estConst.GyroDriftNoise])




%% measurement update

% Get resulting estimates and variances
% Output quantities
posEst = estState.xm(1:2);
linVelEst = ...
oriEst = ...
windEst = ...
driftEst = ...

posVar = ...
linVelVar = ...
oriVar = ...
windVar = ...
driftVar = ...

end


% Function for system dynamics
function [q]= get_q(x_hat,estConst,u)
    % input: xhat vector: px, py, sx, sx, phi, rho, b
    
    C_dh = estConst.dragCoefficientHydr;
    C_da = estConst.dragCoefficientAir;
    C_w = estConst.windVel;
    C_r = estConst.rudderCoefficient; 
    
    % make it easier to read
    sx = x_hat(3);
    sy = x_hat(4);
    phi = x_hat(5);
    rho = x_hat(6);
    
    % returns a 7x1 vector
    % x_dot = q(x_hat, 0, t)
    % all noises = 0 beacuse they are zero mean 
    q = zeros(7,1);
    q(1)= sx; %sx
    q(2)= sy; %sy
    q(3)= cos(phi)*(tanh(u(1))- C_dh(sx^2+sy^2))-C_da(sx-Cw*cos(rho))*sqrt((sx-C_w*cos(rho))^2+(sy-C_w*sin(rho))^2);
    q(4)= sin(phi)*(tanh(u(1))- C_dh(sx^2+sy^2))-C_da(sy-Cw*sin(rho))*sqrt((sx-C_w*cos(rho))^2+(sy-C_w*sin(rho))^2);
    q(5)= C_r*u(2);
    %q(6)= 0  and q(7) =0
    
end 

% Functions for A nad L matrices

function [A] = get_A(x_hat, estConst, u)
    % A is 7x7 matrix: A = dq/dx
    C_dh = estConst.dragCoefficientHydr;
    C_da = estConst.dragCoefficientAir;
    C_w = estConst.windVel;
    C_r = estConst.rudderCoefficient; 
    
    % make it easier to read
    sx = x_hat(3);
    sy = x_hat(4);
    phi = x_hat(5);
    rho = x_hat(6);
    
    A = zeros(7); %7x7 zeros matrix
    A(1,3)= 1; % dq(1)/d(sx)
    A(2,4) = 1; % dq(2)/d(sx)
    
    % 13 = c_da; 7 = C_dh; A = tanh(ut)
    % holy shit this is going to be one very long partial derivative to get
    
    A(3,3) = -(C_da*(sx - C_w*cos(rho))^2)/sqrt((sx - C_w*cos(rho))^2 + (sy - C_w*sin(rho))^2) - C_da*sqrt((sx - C_w*cos(rho))^2 + (sy - C_w*sin(rho))^2) - 2*C_dh*sx*cos(phi); %dq(3)/d(sx)
    A(3,4) = -(C_da*(x - C_w*cos(rho))*(sy - C_w*sin(rho)))/sqrt((sx - C_w*cos(rho))^2 + (sy - C_w*sin(rho))^2) - 2*C_dh*sy*cos(phi); %dq(3)/d(sy)
    A(3,5) = -sin(phi)*(tanh(u(1)) - C_hd*(sx^2 + sy^2)); %dq(3)/d(phi)
    A(3,6) = -C_da*C_w*sin(rho)*sqrt((sx - C_w*cos(rho))^2 + (sy - C_w*sin(rho))^2) - (C_da*(sx - C_w*cos(rho))*(2*C_w*sin(rho)*(sx - C_w*cos(rho)) - 2*C_w*cos(rho)*(sy - C_w*sin(rho))))/(2*sqrt((sx - C_w*cos(rho))^2 + (sy - C_w*sin(rho))^2)); %dq/d(rho)
    
    A(4,3) = -(C_da*(sx - C_w*cos(rho))*(sy - C_w*sin(rho)))/sqrt((sx - C_w*cos(rho))^2 + (sy - C_w*sin(rho))^2) - 2*C_dh*sx*sin(phi);%dq4/d(sx)
    A(4,4) = -(C_da*(sy - C_w*sin(rho))^2)/sqrt((sx - C_w*cos(rho))^2 + (sy - C_w*sin(rho))^2) - C_da*sqrt((sx - C_w*cos(rho))^2 + (sy - C_w*sin(rho))^2) - 2*C_dh*sy*sin(phi);%dq4/d(sy)
    A(4,5) = cos(phi)*(tanh(u(1)) - C_da*(sx^2 + sy^2));
    A(4,6) = C_da*C_w*cos(rho)*sqrt((sx - C_w*cos(rho))^2 + (sy - C_w*sin(rho))^2) - (C_da*(sy - C_w*sin(rho))*(2*C_w*sin(rho)*(sx - C_w*cos(rho)) - 2*C*cos(rho)*(sy - C_w*sin(rho))))/(2*sqrt((sx - C_w*cos(rho))^2 + (sy - C_w*sin(rho))^2)) ; %dq/d(rho)
end 

function [L]= get_L((x_hat, estConst, u)
    % process noise vector v= [vd; vr; vrho; vb]
    
    C_dh = estConst.dragCoefficientHydr;
    C_da = estConst.dragCoefficientAir;
    C_w = estConst.windVel;
    C_r = estConst.rudderCoefficient; 
    
    % make it easier to read
    sx = x_hat(3);
    sy = x_hat(4);
    phi = x_hat(5);
    rho = x_hat(6);
    
    % So L = dq/dv is 7x4 matrix
    L = zeros(7,4);
    %L(1,:), L(2,:) = 0 position not related to noise
    L(3,1) = -cos(phi)*C_dh*(sx^2+sy^2);
    L(4,1) = -sin(phi)*C_dh*(sx^2+sy^2);
    L(5,2)= C_r*u(2);
    L(6,3)= 1;
    L(7,4)= 1;
end