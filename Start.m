%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Taak stangen (8 stangen pantograaf)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data initialization (all data is converted to SI units)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% program data
fig_kin_8bar = 1;        % draw figures of kinematic analysis if 1
fig_dyn_8bar = 1;        % draw figures of dynamic analysis if 1

% kinematic parameters (link lengths)
r1 = 0.493;
r2 = 0.151;
r3 = 0.081;
r6 = 0.494;
r7 = 0.407;
AE = 0.409;
FC = 0.494;
EB = 0.494;
BF = 0.405;
phi1 = 0;
r4 = AE+EB;
r5 = BF+FC;

% r0 = 0.65;
% r1 = 0.3;
% r2 = 0.15;
% r5 = 0.40;
% r6 = 0.58;
% AE = 0.5;
% FC = 0.5;
% EB = 0.5;
% BF = 0.5;
% phi0 = 0;
% r3 = 1;
%r4 = 1;

%DYNAMISCHE PARAMETERS NOG TE DEFINIËREN!!!
% dynamic parameters, defined in a local frame on each of the bars.
X2 = r2/2;               % X coordinates of cog (centre of gravity)
X3 = r3/2;
X4 = r4/2;

Y2 = 0;                  % Y coordinates of cog
Y3 = 0.0102362;
Y4 = 0;

m2 = r3*1.76;
m3 = r3*1.76;
m4 = r4*0.54;

J2 = m2*r3^2/12;
J3 = m3*r3^2/12;
J4 = m4*r4^2/12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Determination of Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t_begin = 0;                   % start time of simulation
t_end = 20;                    % end time of simulation
Ts = 0.05;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector

% initialization of driver
omega = pi/10;
A = pi; % amplitude van hoek
% phi1=A*sin(omega*t+theta);
% dphi1=omega*A*cos(omega*t);
% ddphi1=-omega^2*A*sin(omega*t); % omega ct
phi2 = -omega*t+pi/2 -0.25;
dphi2 = -omega*ones(size(t));
ddphi2 =0*t;

% position analysis
y8_init = 0;
phi3_init = pi-phi2(1);
phi4_init = pi/4;    % initial condition for first step of position analysis with fsolve (phi3 and phi4)
phi5_init = 2*pi/3;  % VERY IMPORTANT because it determines which branch of the mechanism you're in
phi6_init = pi/4;
phi7_init = 2*pi/3;
x8_init = r3*cos(phi3_init)+r6*cos(phi6_init)-FC*cos(phi5_init);
% calculation of the kinematics (see kin_8bar.m)
[y8,phi4,phi5,phi6,phi7,x8,dy8,dphi4,dphi5,dphi6,dphi7,dx8,ddy8,ddphi4,ddphi5,ddphi6,ddphi7,ddx8] = kinematics_8bar(r1,r2,r3,r4,r5,r6,r7,AE,FC,phi2,dphi2,ddphi2,y8_init,phi4_init,phi5_init,phi6_init,phi7_init,x8_init,t,fig_kin_8bar);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % STEP 2. Dynamics Calculation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % calculation of the dynamics (see dyn_4bar.m)
% [F_P_x,F_Q_x,F_R_x,F_S_x,F_P_y,F_Q_y,F_R_y,F_S_y,M_P] = dynamics_4bar(phi2,phi3,phi4,dphi2,dphi3,dphi4,ddphi2,ddphi3,ddphi4,r2,r3,r4, ...
%   m2,m3,m4,X2,X3,X4,Y2,Y3,Y4,J2,J3,J4,t,fig_dyn_4bar);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % STEP 3. Movie
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% load fourbar_movie Movie
% movie(Movie)
