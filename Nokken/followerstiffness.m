%CAMS ASSIGNMENT subtask 4.2
% Determine kf and check condition for approximate single-rise analysis
clear variables
close all 

cam = load('campower_ecc32');

%%%%%%%%%%%%%%%%% Segment 1 of critical fall (150-160)%%%%%%%%%%%%%
%define parameters
h = 10*0.001;
thfall1 = 150;
thdwell1 = 160;
thend1 = 0;
zeta = 0.063;
m = 25;
ks = cam.springconstant*10^3;
tcyc = 2;
omega = pi;

%calculate t1 and t2
t11 = tcyc*(thdwell1-thfall1)/360;
t21 = tcyc*(thend1+360-thdwell1)/360;

%calculate kf bordervalue and lambda
kf_b1 = (m*((0.75*2*pi)/(zeta*t11))^2-ks);
kf1 = 4.55*10^7;
lam_b1 = (t11/(2*pi))*sqrt((kf_b1+ks)/m);
lam1 = (t11/(2*pi))*sqrt((kf1+ks)/m);

%check two conditions
t2_min_b1 = t11*(1+log(0.02)/(-zeta*2*pi*lam_b1));
t2_min1 = t11*(1+log(0.02)/(-zeta*2*pi*lam1));
cond_b1 = zeta*lam_b1; 
cond1 = zeta*lam1;

%%%%%%%%%%%%%%%%% Segment 2 of critical fall (160-180)%%%%%%%%%%%%%
h = 20*0.001;
thfall2 = 160;
thdwell2 = 180;
thend2 = 20;
zeta = 0.063;
m = 25;
ks = cam.springconstant;
tcyc = 2;
omega = pi;

%calculate t1 and t2
t12 = tcyc*(thdwell2-thfall2)/360;
t22 = tcyc*(thend2+360-thdwell2)/360;

%calculate kf bordervalue and lambda
kf_b2 = (m*((0.75*2*pi)/(zeta*t12))^2-ks);
kf2 = 1.14*10^7;
lam_b2 = (t12/(2*pi))*sqrt((kf_b2+ks)/m);
lam2 = (t12/(2*pi))*sqrt((kf2+ks)/m);

%check two conditions
t2_min_b2 = t12*(1+log(0.02)/(-zeta*2*pi*lam_b2));
t2_min2 = t12*(1+log(0.02)/(-zeta*2*pi*lam2));
cond_b2 = zeta*lam_b2; 
cond2 = zeta*lam2;

%%%%%%%%%%%%%%%%% Total critical fall (150-180) with chosen kf %%%%%%%%%%%%%
t1 = 1/6;
kf = 4.55*10^7;
lambda = (t1/(2*pi))*sqrt((kf+ks)/m);
cond = zeta*lambda;
t2_min = t1*(1+log(0.02)/(-zeta*2*pi*lambda));





