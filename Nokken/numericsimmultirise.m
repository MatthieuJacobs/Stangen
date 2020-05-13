% CAM ASSIGNMENT subtask 4.3
% multirise numerical simulation of complete cam-follower
clear variables
close all 


cam = load('campower_ecc32');
camstiff = load('stiffness_values');

%define parameters
h=30; %mm
m=25; %kg
zeta = 0.063; 
kf = camstiff.kf; 
t1 = 1/6; %time to execute fall (150° to 180°)
t2 = 10/9;%time to execute dwell after fall (180° to 360°+20°)
wn = sqrt(kf/m);
tn = 2*pi/wn;
t3 = 2*(360-180)/360;
ks = 4000;
lambda = t1/tn;


%define transfer function
dividend =  (2*pi*lambda)^2;
divisor =  [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% numerical solution single rise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plots from 150° to 360°
theta0 = 1;
theta_dot0 = 0; 
[A,B,C,D] = tf2ss(dividend,divisor);
X0 = [1/C(2)*theta_dot0;1/C(2)*theta0];

tau_s = linspace(0,t2/t1,21000); 
theta_s = cam.S(15000:35999)/h;
gamma_s = lsim(A,B,C,D,theta_s,tau_s,X0);
gamma_s = transpose(gamma_s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% numerical solution multirise one period %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sys = tf(dividend, divisor); 
tau = linspace(0,1,36000); 
theta = cam.S/h;

gamma = lsim(sys,theta,tau);
gamma = transpose(gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% numerical solution multirise four periods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau2 = linspace(0,4,36000*4); 
theta2 = [cam.S/h,cam.S/h,cam.S/h,cam.S/h]; 

gamma2 = lsim(sys,theta2,tau2);
gamma2 = transpose(gamma2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot difference multi/single %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
taud = linspace(0,t3/t1,21000);
gamma1 = gamma2(1,36000*3+15000:36000*3+35999);
gammad = gamma_s - gamma1;


figure;
hold on
plot(tau,theta)
plot(tau,gamma)
legend('input','output')
hold off

%four periods 
figure;
hold on 
plot(tau2,gamma2-theta2)
xlabel('\tau(-)')
ylabel('\gamma(\tau)')
hold off

figure;
tiledlayout('flow')
nexttile
hold on
plot(taud,gamma_s)
plot(taud,gamma(1,15000:35999))
plot(taud,theta_s)
legend('single rise','multi rise','input')
hold off

nexttile
hold on 
plot(taud,gammad)
legend('difference multi-rise/single-rise')
hold off

%one period
figure
hold on 
plot(tau,gamma-theta)
legend('deviation from input')
hold off

