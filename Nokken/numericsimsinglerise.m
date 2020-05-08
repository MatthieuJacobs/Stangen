% CAM ASSIGNMENT subtask 4.3
% numerical simulation of cam-follower
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
sys = tf(dividend,divisor);

%define 
tau = [0,t2/t1,230]; %230 degrees we want to analyse
theta = 0;%degrees 



