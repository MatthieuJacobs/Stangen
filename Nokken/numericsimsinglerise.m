% CAM ASSIGNMENT subtask 4.3
% numerical simulation of cam-follower
clear variables
close all

cam = load('campower_ecc32');
camfol = load('stiffness_values.mat');

%define parameters
h=30;
zeta = 0.063;
kf = 4.55*10^8;
lambda = 35.78;
t1 = 1/6;
t2 = 10/9;
t3 = 2*(360-180)/360;

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
tau1 = linspace(0,t3/t1,21000); %230 degrees we want to analyse
theta1 = [cam.S(15000:35999)/h];
gamma1 = lsim(A,B,C,D,theta1,tau1,X0);
gamma1 = transpose(gamma1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% approximate analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define parameters
Q = (2*pi)^2;
N = 3;
tau2 = linspace(0,t3/t1,21000);
theta2 = ((Q.*(tau2-1).^N)/factorial(N)+1)/h;
gamma2 = lsim(A,B,C,D,theta2,tau2,X0);
gamma2 = transpose(gamma2);
A1 = Q/((2*pi*lambda)^N);


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
yline(0);
legend(': \gamma_a  -  \gamma_n')
xlabel('\tau(-)')



