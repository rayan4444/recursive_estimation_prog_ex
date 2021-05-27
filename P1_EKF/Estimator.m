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
    windEst = 0;
    
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
    return;
end

%% Estimator iteration.
% get time since last estimator update
dt = tm - estState.tm;
estState.tm = tm; % update measurement update time

%% prior update

% Initial state
xP_0 =[estState.xm;reshape(estState.Pm,[49,1])];

% Solve ODE
[~ , sol ] = ode45( @(t,estState_simple) prior_update_ODE(t, estState_simple,estConst, actuate),[tm-dt tm], xP_0);

% Predicted intermediary values
xp = sol(end, 1:7)';
Pp = reshape(xP(end,8:56), [7,7]);

%% measurement update


px =estState.xm(1);
py = estState.xm(2);
phi = estState.xm(5);
b = estState.xm(7);

%get R matrix
R  = diag([estConst.DistNoiseA estConst.DistNoiseB estConst.DistNoiseC estConst.GyroNoise estConst.CompassNoise]);

% get M matrix
M = eye(5); %5x5 matrix

% get h matrix ( measurement model)
h = [sqrt((px- estConst.pos_radioA(1))^2+(py-estConst.pos_radioA(2))^2);
    sqrt((px- estConst.pos_radioB(1))^2+(py-estConst.pos_radioB(2))^2);
    sqrt((px- estConst.pos_radioC(1))^2+(py-estConst.pos_radioC(2))^2);
    phi+b;
    phi]; %5x1 matrix

% get H matrix
H = zeros(5,7); %5x7 matrix
H(1,1) = (px- estConst.pos_radioA(1))/h(1);
H(1,2) = (py-estConst.pos_radioA(2))/h(1);
H(2,1) = (px- estConst.pos_radioB(1))/h(2);
H(2,2) = (py-estConst.pos_radioB(2))/h(2);
H(3,1) = (px- estConst.pos_radioC(1))/h(3);
H(3,2) = (py-estConst.pos_radioC(2))/h(3);
H(4,5) = 1; 
H(4,7)= 1;
H(5,5)= 1;

% Measurement at station C is not always available
% measurement = inf when not available
if isinf(sense(3))
    % replace affected values in matrices by a blank
    R(3,:) = [];
    R(:,3) = [];
    
    M(3,:) = [];
    M(:,3) = [];
    
    H(3,:) = [];
    
    sense(3) = [];
    h(3) = [];
    
end

% get Kalman gain matrix
K= Pp*H'/(H*Pp*H'+M*R*M'); %7x5 matrix

% Update estimate with measurement
xm = xp + K*(sense'-h); %7x1 matrix
Pm = (eye(7)-K*H)*Pp ;%7x7 matrix

% Get resulting estimates and variancese
estState.xm = xm;
estState.Pm = Pm;

% Output quantities
posEst = estState.xm(1:2);
linVelEst = estState.xm(3:4);
oriEst = estState.xm(5);
windEst = estState.xm(6);
driftEst = estState.xm(7);

posVar = [estState.Pm(1,1) estState.Pm(2,2)];
linVelVar = [estState.Pm(3,3) estState.Pm(4,4)];
oriVar = estState.Pm(5,5);
windVar = estState.Pm(6,6);
driftVar = estState.Pm(7,7);

end


%% Functions used in the prior update step 
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

function [L]= get_L(x_hat, estConst, u)
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

function [xP_dot]= prior_update_ODE(t,estState_simple,estConst, actuate)
    
    xm = estState_simple(1:7);
    x_dot = get_q(xm, estConst,actuate); %7x1 matrix
    
    % find A and L matrices
    A = get_A(xm, estConst, actuate); %7x7 matrix
    L = get_L(xm, estConst, actuate); %7x4 matrix 
    Qc =diag([estConst.DragNoise estconst.RudderNoise estConst.GyroDriftNoise]); % 4x4 diagonal matrix
    
    P = reshape(estState_simple(8:56), [7,7]);
    P_dot = reshape(A*P + P*A'+L*Qc*L', [49, 1]); % flatten 7x7 matrix in to 49x1
    
    % Combine and reshape matrices for solving ODE
    
    xP_dot = [x_dot; P_dot]; %(7+49)x1 matrix
end 

