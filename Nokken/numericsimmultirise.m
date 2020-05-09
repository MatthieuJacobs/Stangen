% CAM ASSIGNMENT subtask 4.3
% approximate single rise analysis of cam-follower
clear variables
close all

cam = load('campower_ecc32');
camfol = load('stiffnessapprox.mat');

%define parameters
h = 30;
thfall = 150;
thdwell = 180;
thend = 20;
zeta = 0.063;
m = 25;
ks = cam.springconstant;
kf = camfol.kf;
lambda = camfol.lam;
tcyc = 2;
omega = pi;
t1 = camfol.t1;
t2 = camfol.t2;

%define transfer function
dividend =  (2*pi*lambda)^2;
divisor =  [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];

%define start conditions
theta0 = 1;
theta_dot0 = 0; 
[A,B,C,D] = tf2ss(dividend,divisor);
X0 = [1/C(2)*theta_dot0;1/C(2)*theta0];

%define exact input and time 
tau1 = linspace(0,t2/t1,23000); %230 degrees we want to analyse
theta1 = cam.S(15000:36000)/h;
theta1 = [theta1,zeros(1,1999)]; %data points of cam profile from 150° to 380°
gamma1 = lsim(A,B,C,D,theta1,tau1,X0);
gamma1 = transpose(gamma1);
