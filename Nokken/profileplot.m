% CAMS ASSIGNMENT subtask 1
% plot motion law 
% -------------------------------
clear variables
close all 
% load correct motion law
load('motionlaw_pg1','S','Vdegree','Adegree','thetadegree');

% acceleration = acc
% velocity = Vdegree
% profile = S

x = thetadegree;

tiledlayout(3,1)

%profile
nexttile
plot1 = plot(x,S);
legend('S')
title('cam profile')

%velocity
nexttile
plot2 = plot(x,Vdegree);
legend('velocity')
title('velocity')

%acceleration
nexttile
plot3 = plot(x,Adegree);
legend('acceleration')
title('acceleration')
box on;




















 