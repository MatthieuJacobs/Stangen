clear variables
close all

out = load('campower_ecc32.mat');

%%% CRITICAL RISE

beta = 0.5236; %[rad] 30 graden
w = out.w;
t1 = beta/w;

%%% NUMERICAL ANALYSIS

m = out.mass; %kg
h = 30; %mm
kf = 4.55*10^7; %N/m
wn = sqrt(kf/m);
tn = 2*pi/wn; 
lambda = t1/tn;
zeta = 0.063;

% critical_rise = out.S(4501:11376);
% we look at the critical rise  + following dwell: 
% delta_theta_rise = 68,75° , delta_theta_tot = 125,00°

tau = (0 :1/3000 : 3/3000);

% theta = heffing/h

theta = out.S(4500:17000)/h;

teller = (2*pi*lambda)^2;
noemer = [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];
theta0 = 1;
theta_dot0 = 0; 
[A,B,C,D] = tf2ss(teller,noemer);
X0 = [1/C(2)*theta_dot0;1/C(2)*theta0];

gamma = lsim(A,B,C,D,theta,tau,X0);

figure
plot(tau,theta); hold on; plot(tau,gamma);
xlabel('\tau [-]'),
ylabel('\theta [-] en \gamma [-]')
hold off

error = gamma.' - theta;
figure
plot(tau,error)
axis tight
xlabel('\tau [-]'),
ylabel('verschil \gamma en \theta')

%------------------------------------------------------------------------------------------

%%% A ANALYTISCH BEREKENEN


lambda_d = lambda*sqrt(1-zeta^2);
x0 = gamma(6875)-1; 
dx0 = (gamma(6876)-gamma(6874))/(2/6875);

A1_analytisch = sqrt(((x0*2*pi*lambda_d)^2+(dx0+zeta*2*pi*lambda*x0)^2)/(2*pi*lambda_d)^2);

phi_analytisch = atan(- (dx0 + zeta*(2*pi*lambda)*x0)/(x0*(2*pi*lambda_d)));

gamma_analytisch =  1 + A1_analytisch.*exp(-zeta*(2*pi*lambda).*(tau-1)).*cos((2*pi*lambda_d).*(tau-1)+ phi_analytisch);


figure
plot(tau,gamma_analytisch);
axis tight
xlabel('\tau [-]'),
ylabel('\gamma_{analytisch}')

%%% APPROXIMATE ANALYSIS

Q = (2*pi)^2;
N = 3;   

A1_approx = Q/(2*pi*lambda)^N * sqrt(1/(1-zeta^2));

epsilon = abs((A1_analytisch-A1_approx)/A1_analytisch);

% --------------------------------------------------------------------------------------------

%%% MULTI RISE ANALYSIS

T = 2*pi/w;
lambda2 = T/tn;

teller2 = (2*pi*lambda2)^2;
noemer2 = [1, 2*zeta*(2*pi*lambda2), (2*pi*lambda2)^2];
sys2 = tf(teller2, noemer2);

tau2 = (0:1/36000:25-1/36000); % from 0 to 25 instead of 0 to 1 because we look at the steady state respons
theta2 = repmat(out.S/h,1,25); %25 times S to get the steady state respons at the end
gamma2 =lsim(sys2,theta2,tau2);

int_mr = 24*36000:25*36000-1; %interval to look at steady state respons

figure 
plot(tau2(int_mr),gamma2(int_mr),'b',tau2(int_mr),theta2(int_mr),'r')
title('Multi-rise analysis')
xlabel('tau [-]')
ylabel('\theta [-] en \gamma [-]')
legend('\gamma','\theta')

%%% Comparision single rise and multi rise

gamma_mr = gamma2(24*36000+4500:24*36000+17000);
figure
plot(tau,gamma-gamma_mr)
title('Difference \gamma_{single rise} and \gamma_{multi rise}')
xlabel('tau [-]')
ylabel('\gamma [-]')

% graph of the transient respons 

gamma0 = gamma2(1:1*36000);
gamma24 = gamma2(24*36000+1:25*36000);
figure
plot(tau2(1:1*36000),gamma0-gamma24) % = output first cycle - output 24th cycle
title('Transient response')
xlabel('tau [-]')
ylabel('\gamma [-]')


%%% Vibrations as a contact force

tau3 = 0:1/36000:1-1/36000;
theta3 = out.S;
gamma3 =lsim(sys,theta3,tau3);
size(theta3)
size(gamma3)
force_follower = kf .* (-gamma3+transpose(theta3))*0.025;
size(force_follower) 

figure
plot(tau3, force_follower)
xlabel('\tau [-]'),
ylabel('F_{volger} [N]')

force_normal =  out.normalforce_tot;
size(force_normal)
c = cos(out.pressure_angle);
size(c)
x = force_follower .'/ cos(out.pressure_angle);
size(x)
F_tot = force_normal + force_follower .'/ cos(out.pressure_angle);
figure
plot(tau3, F_tot)
xlabel('\tau [-]'),
ylabel('F_{tot} [N]')

