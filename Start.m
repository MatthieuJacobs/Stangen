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
fig_kin_8bar = 0;        % draw figures of kinematic analysis if 1
fig_dyn_8bar = 0;        % draw figures of dynamic analysis if 1
fig_forward = 1;

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
X5 = r5/2;
X6 = r6/2;
X7 = r7/2;

m2 = r2*1/2;
m3 = r3*1/2;
m4 = r4*1/2;
m5 = r5*1/2;
m6 = r6*1/2;
m7 = r7*1/2;
mt2 = pi*r1*1/2;
mt3 = mt2;
m8 = 1;

J2 = m2*r2^2/12;
J3 = m3*r3^2/12;
J4 = m4*r4^2/12;
J5 = m5*r5^2/12;
J6 = m6*r6^2/12;
J7 = m7*r7^2/12;
Jt2 = mt2*r1^2/8;
Jt3 = mt3*r1^2/8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Determination of Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t_begin = 0;                   % start time of simulation
t_end = 20;                    % end time of simulation
Ts = 0.05;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector

% initialization of driver
omega = pi/5;
A = pi; % amplitude van hoek
% phi1=A*sin(omega*t+theta);
% dphi1=omega*A*cos(omega*t);
% ddphi1=-omega^2*A*sin(omega*t); % omega ct
phi2 = -omega*t+pi/2 -0.25;
dphi2 = -omega*ones(size(t));
ddphi2 =0*t;
phi3 = pi-phi2;
dphi3 = -dphi2;
ddphi3 = -ddphi2;
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

%calculation of the dynamics (see dyn_4bar.m)
[F2x,F2y,F3x,F3y,F42x,F42y,F63x,F63y,F73x,F73y,F54x,F54y,F74x,F74y,F65x,F65y,M2,M3] = dynamics_8bar(phi2,phi3,phi4,phi5,phi6,phi7,x8,y8,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dx8,dy8,...
    ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddx8,ddy8,...
    r2,r3,r4,r5,r6,r7,AE,FC,r1, ...
    m2,m3,m4,m5,m6,m7,m8,...
    J2,J3,J4,J5,J6,J7,Jt2,Jt3,t,fig_dyn_8bar);

% calculation of the dynamics with gravity implemented
% [F2x,F2y,F3x,F3y,F42x,F42y,F63x,F63y,F73x,F73y,F54x,F54y,F74x,F74y,F65x,F65y,M2,M3] = dynamics_8bar_gravity(phi2,phi3,phi4,phi5,phi6,phi7,x8,y8,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dx8,dy8,...
%     ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddx8,ddy8,...
%     r2,r3,r4,r5,r6,r7,AE,FC,r1, ...
%     m2,m3,m4,m5,m6,m7,m8,...
%     J2,J3,J4,J5,J6,J7,Jt2,Jt3,t,fig_dyn_8bar);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % STEP 3. Movie
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
% load fourbar_movie Movie
% movie(Movie)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Forward Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute velocities and internal forces using driving moments and initial
% position
[phi2,phi3,phi4,phi5,phi6,phi7,x8,y8,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dx8,dy8,...
    ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddx8,ddy8,...
    F42x,F42y,F63x,F63y,F73x,F73y,F65x,F65y,F58x,F58y,F74x,F74y,F54x,F54y] = ...
    forward_dynamics(M2,M3,...
    r2,r3,r4,r5,r6,r7,AE,FC,r1,m2,m3,m4,m5,m6,m7,m8,...
    J2,J3,J4,J5,J6,J7,Jt2,Jt3,t,fig_forward);