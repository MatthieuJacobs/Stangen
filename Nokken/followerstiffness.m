%CAMS ASSIGNMENT subtask 4.2
% Determine kf and check condition for approximate single-rise analysis
clear variables
close all 

cam = load('campower_ecc32');

%define parameters
h = 30*0.001;
thfall = 150;
thdwell = 180;
thend = 20;
zeta = 0.063;
m = 25;
ks = cam.springconstant;
tcyc = 2;
omega = pi;


%calculate t1 and t2
t1 = tcyc*(thdwell-thfall)/360;
t2 = tcyc*(thend+360-thdwell)/360;

%calculate kf bordervalue and lambda
%kf = ((0.75*2*pi)/(zeta*t1))^2*m-ks;
kf = 5.1*10^6;
lam = (t1/(2*pi))*sqrt((kf+ks)/m);

%check condition
t2_min = t1*(1+log(0.02)/(-zeta*2*pi*lam));
cond = zeta*lam; 


