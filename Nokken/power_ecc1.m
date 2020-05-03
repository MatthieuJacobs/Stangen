% CAMS ASSIGNMENT Subtask 3
% compute instanteneous power with eccentric profile
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

%--------------------------------------------------------%
%load profile without eccentricity
cam_norm=load('campower_ecc0');

%load normal force, pressure angle and pitch radius
N=cam_norm.normalforce_tot;
alpha=cam_norm.pressure_angle;
omega=cam_norm.w;

%compute pitch radius
R0=(cam_norm.xpitch.^2+cam_norm.ypitch.^2).^(1/2)*0.001;

%compute power without eccentricity 
P2=N.*sin(alpha).*R0*omega;

%--------------------------------------------------------%
%plot both cases 
figure
hold on 
plot(cam_ecc.thetadegree,P1)
plot(cam_norm.thetadegree,P2,'-')
title('Combination power with/without excentricity')
legend('with excentricity','without excentricity')
hold off



