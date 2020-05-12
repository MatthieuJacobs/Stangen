% CAM ASSIGNMENT subtask 4.3
% Computation of contact forces
clear variables
close all 


cam = load('campower_ecc32');
camfol = load('stiffness_values');

%define parameters
h=30;
m=25;
zeta = 0.063;
kf = 4.55*10^7;
t1 = 1/6;
t2 = 10/9;
t3 = 2*(360-180)/360;
ks = 4;
lambda = (t1/(2*pi))*sqrt((camfol.kf+ks)/m);


%define transfer function
dividend = (2*pi*lambda)^2;
divisor = [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];
sys = tf(dividend, divisor);

tau = linspace(0,1,36000);
theta = cam.S/h; 
gamma = lsim(sys,theta,tau);
gamma = transpose(gamma);

%dimensionless deviation from input
dev = gamma-theta;
devabs = dev*0.030;
extranorm = (devabs)*kf./cos(cam.pressure_angle);

%add this to existing normal force 
normal = cam.normalforce_tot;
normtot = normal-extranorm;

%plots
x= cam.thetadegree;

figure
tiledlayout(2,1)

nexttile
hold on
plot(x,devabs);
legend('absolute deviation [m]')
hold off

nexttile
hold on
plot(x,extranorm*(-1));
plot(x,normal)
plot(x,normtot);
legend('addition contact force','initial contact force','total contact force')
hold off