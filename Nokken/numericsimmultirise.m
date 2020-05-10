% CAM ASSIGNMENT subtask 4.3
% multirise numerical simulation of complete cam-follower
clear variables


cam = load('campower_ecc32');
camfol = load('stiffnessapprox.mat');

%define parameters
h = 30;
thfall = 150;
thdwell = 180;
thend = 20;
zeta = 0.063;
m = 25;
t1 = camfol.t1;
t2 = camfol.t2;
ks = cam.springconstant;
kf = camfol.kf;
lambda = (t1/(2*pi))*sqrt((kf*10+ks)/m);
tcyc = 2;
omega = pi;


%define transfer function
dividend =  (2*pi*lambda)^2;
divisor =  [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% numerical solution single rise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lots from 150° to 360°
theta0 = 1;
theta_dot0 = 0; 
[A,B,C,D] = tf2ss(dividend,divisor);
X0 = [1/C(2)*theta_dot0;1/C(2)*theta0];

taus = linspace(0,t2/t1,21000); 
thetas = cam.S(15000:35999)/h;
gammas = lsim(A,B,C,D,thetas,taus,X0);
gammas = transpose(gammas);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% numerical solution multirise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dividend = (2*pi*lambda)^2;
divisor = [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];
sys = tf(dividend, divisor);

tau = linspace(0,1,36000); % vul zelf in, dimensieloze tijd
theta = cam.S/h; % vul zelf in, dimensieloze heffing theta als functie van tau
gamma = lsim(sys,theta,tau);
gamma = transpose(gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot difference multi/single %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
taud = linspace(0,t2/t1,21000);
gamma1 = gamma(1,15000:35999);
gammad = gammas - gamma1;

figure;
hold on
plot(tau,theta)
plot(tau,gamma)
legend('input','output')
hold off

figure;
tiledlayout('flow')
nexttile
hold on
plot(taud,gammas)
plot(taud,gamma(1,15000:35999))
plot(taud,thetas)
legend('single rise','multi rise','input')
hold off

nexttile
hold on 
plot(taud,gammad)
hold off
legend('difference')

