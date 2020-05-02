% CAMS ASSIGNMENT Subtask 3
% Determination optimal spring constant and preload

clear variables
close all 

load('springconstant1.mat')
bound = 

% all forces independent of spring or preload
F_fixed = normalforce_acc + normalforce_load; 
camprofile = S;
spring_length = camprofile;

% contact equation

CE(F_pre,k) = F_fixed + k*spring_length + F_pre == 0 ; 
plot(CE(20,10),k);



