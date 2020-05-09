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

%define start conditions
theta0 = 1;
theta_dot0 = 0; 
[A,B,C,D] = tf2ss(dividend,divisor);
X0 = [1/C(2)*theta_dot0;1/C(2)*theta0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% numerical solution with transfer function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define exact input and time 
tau1 = linspace(0,t2/t1,23000); %230 degrees we want to analyse
theta1 = cam.S(15000:36000)/h;
theta1 = [theta1,zeros(1,1999)]; %data points of cam profile from 150° to 380°
gamma1 = lsim(A,B,C,D,theta1,tau1,X0);
gamma1 = transpose(gamma1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% approximate analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define parameters
Q = (2*pi)^2;
N = 3;
tau2 = linspace(0,t2/t1,23000);
theta2 = ((Q.*(tau2-1).^N)/factorial(N)+1)/h;
gamma2 = lsim(A,B,C,D,theta2,tau2,X0);
gamma2 = transpose(gamma2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure numerical solution
figure;
tiledlayout(1,2)

nexttile
hold on
title('numerical solution')
plot(tau1,gamma1);
plot(tau1,theta1);
legend('dimensionless output: \gamma_n(\tau)','dimensionless input: \theta_n(\tau)')
xlabel('\tau(-)')
hold off

nexttile
plot(tau1,gamma1-theta1);
title('difference input/output numerical solution');
ylabel('\gamma_n(\tau) -  \theta_n(\tau)')
xlabel('\tau(-)')

%figure approximate solution
figure;
tiledlayout(1,2);

nexttile
hold on
title('approximate solution')
plot(tau2,gamma2);
plot(tau2,theta2);
legend('dimensionless output: \gamma_a(\tau)','dimensionless input: \theta_a(\tau)')
xlabel('\tau(-)')
hold off

nexttile
plot(tau2,gamma2-theta2);
ylabel('\gamma_a(\tau) -  \theta_a(\tau)')
xlabel('\tau(-)')
title('difference input/output approximate solution');

%figure difference between methods plot
figure;

title('Difference approximate and numerical output')
plot(tau1, gamma1-gamma2);
legend(': \gamma_a  -  \gamma_n')
xlabel('\tau(-)')



