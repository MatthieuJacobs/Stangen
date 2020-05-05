% CAMS ASSIGNMENT Subtask 3
% Determine flywheel dimensions
%--------------------------------------------------------%
clear variables
close all 

%determine torque 
cam_pow=load('campower_ecc32');

%load normal force, pressure angle, eccentricity and pitch radius
N=cam_pow.normalforce_tot;
alpha=cam_pow.pressure_angle;
e=cam_pow.exc*0.001;
omega=cam_pow.w;
R0=(cam_pow.xpitch.^2+cam_pow.ypitch.^2).^(1/2)*0.001;


%compute torque:
T1 = (N.*cos(alpha).*e + N.*sin(alpha).*(R0.^2-e.^2).^(1/2)); %instanteneous torque 
T_avg = mean(T1); %torque the engine will produce in regime 
x= cam_pow.thetadegree;

%compute speed variation and work surplus
delta_T = T1-T_avg;

A = zeros(size(x));
for i = 2:36000
    A(i) = (trapz(x(1:i),delta_T(1:i))); 
end

[~,cM] = max(A);
[~,cm] = min(A);
tM = cM/100;
tm = cm/100;

%compute Amax and I 
A_max = trapz(x(min(tm,tM):max(tm,tM)),T1(min(tm,tM):max(tm,tM))-T_avg);
K=0.1;
I = abs((A_max)/((omega^2)*K));

%compute weight of flywheel
R_max = max(R0);
m = I*2/((R_max)^2);


%compute kinetic energy in regime 
E_kin = I*(omega^2)/2;


%plot torque demand and average power 
figure 
tiledlayout(2,1)

nexttile
hold on
plot(x,T1)
plot(x,delta_T)
legend('Instanteneous torque','Torque variation')
hold off

nexttile
hold on
plot(x,A)
xline(tM)
xline(tm)
legend('work function')
hold off






