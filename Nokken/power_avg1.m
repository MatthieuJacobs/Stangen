% CAMS ASSIGNMENT Subtask 3
% compute average power
%--------------------------------------------------------%
clear variables
close all 

%load profile with eccentricity
cam_ecc=load('campower_ecc32');

%load normal force, pressure angle and eccentricity
N=cam_ecc.normalforce_tot;
alpha=cam_ecc.pressure_angle;
e=cam_ecc.exc*0.001;
omega=cam_ecc.w;

%compute pitch radius
R0=(cam_ecc.xpitch.^2+cam_ecc.ypitch.^2).^(1/2)*0.001;

%compute power with eccentricity
P1=(N.*cos(alpha).*e + N.*sin(alpha).*(R0.^2-e.^2).^(1/2)).*omega;

%compute energy consumed in one rotation and average power 
P_avg = mean(P1); 
t = 60/cam_ecc.rpm;
P_avg_list = repmat(P_avg/t,1,36000);
x= cam_ecc.thetadegree;

hold on 
plot(x,P1)
plot(x,P_avg_list)
legend('instanteneous power','average power')
