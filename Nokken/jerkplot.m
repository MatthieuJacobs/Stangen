% CAMS ASSIGNMENT Subtask 1
% compute jerk and Djerk
% -------------------------------
clear variables


load('motionlaw_pg1')
x = thetadegree; 

%compute jerk and Djerk
jerk = diff(A);
jerk = [jerk,0];

Djerk = diff(jerk);
Djerk = [Djerk,0];

%jerk plot
tiledlayout(2,1)
nexttile;
plot(x,jerk); 
title('jerk plot')
legend('jerk');
xlabel('theta [degrees]');
box on; 

%Djerk plot
nexttile;
plot(x,Djerk);
title('Djerk plot')
legend('Djerk');
xlabel('theta [degrees]');
box on; 

