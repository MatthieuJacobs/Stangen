%%%%%%%%%%%%%%%%% Forward dynamic analysis %%%%%%%%%%%%%%%%%%%%
function [phi2,phi3,phi4,phi5,phi6,phi7,x8,y8,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dx8,dy8,...
    ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddx8,ddy8,...
    F42x,F42y,F63x,F63y,F73x,F73y,F65x,F65y,F58x,F58y,F74x,F74y,F54x,F54y] = ...
    forward_dynamics(M2,M3,...
    r2,r3,r4,r5,r6,r7,AE,FC,r1,m2,m3,m4,m5,m6,m7,m8,...
    J2,J3,J4,J5,J6,J7,Jt2,Jt3,t,fig_forward)

% fsolve options (help fsolve, help optimset)
optim_options = optimset('Display','off','TolFun',10^-12);
% *** loop over positions ***
Ts = t(2) - t(1);      % timestep
t_size = size(t,1);    % number of simulation steps

% All unknowns
phi2 = zeros(size(M2));
phi3 = zeros(size(phi2));
phi4 = zeros(size(phi2));
phi5 = zeros(size(phi2));
phi6 = zeros(size(phi2));
phi7 = zeros(size(phi2));
x8 = zeros(size(phi2));
y8 = zeros(size(phi2));
dphi2 = zeros(size(phi2));
dphi3 = zeros(size(phi2));
dphi4 = zeros(size(phi2));
dphi5 = zeros(size(phi2));
dphi6 = zeros(size(phi2));
dphi7 = zeros(size(phi2));
dx8 = zeros(size(phi2));
dy8 = zeros(size(phi2));
ddphi2 = zeros(size(phi2));
ddphi3 = zeros(size(phi2));
ddphi4 = zeros(size(phi2));
ddphi5 = zeros(size(phi2));
ddphi6 = zeros(size(phi2));
ddphi7 = zeros(size(phi2));
ddx8 = zeros(size(phi2));
ddy8 = zeros(size(phi2));
F2x = zeros(size(phi2));
F2y = zeros(size(phi2));
F3x = zeros(size(phi2));
F3y = zeros(size(phi2));
F42x = zeros(size(phi2));
F42y = zeros(size(phi2));
F63x = zeros(size(phi2));
F63y = zeros(size(phi2));
F73x  = zeros(size(phi2));
F73y = zeros(size(phi2));
F54x = zeros(size(phi2));
F54y = zeros(size(phi2));
F74x = zeros(size(phi2));
F74y = zeros(size(phi2));
F65x = zeros(size(phi2));
F65y = zeros(size(phi2));
F58x = zeros(size(phi2));
F58y = zeros(size(phi2));
% Initialize positions
phi2_init = pi/2-0.25;
dphi2_init = -pi/5;
phi3_init = pi-phi2_init;
phi4_init = pi/4;    
phi5_init = 2*pi/3; 
phi6_init = pi/4;
phi7_init = 2*pi/3;
x8_init = r3*cos(phi3_init)+r6*cos(phi6_init)-FC*cos(phi5_init);
y8_init = 0;
for k=1:t_size
    
%% Step1: Compute next position and speed of the mechanism
% Find next positions
    [x, fval, exitflag]=fsolve('loop_closure_eqs_8bar',[y8_init,phi4_init phi5_init,phi6_init,phi7_init,x8_init]',optim_options,phi2_init,r1,r2,r3,r4,r5,r6,r7,AE,FC);
    if (exitflag ~= 1)
        display 'The fsolve exit flag was not 1, probably no convergence!'
        exitflag
        k
    else
    end
    
    % save results of fsolve
    y8(k) = x(1);
    phi4(k)=x(2);
    phi5(k)=x(3);
    phi6(k) = x(4);
    phi7(k) = x(5);
    x8(k) = x(6); 
    phi2(k) = phi2_init;
    phi3(k) = pi-phi2(k);


% Find next velocities
dphi2(k) = dphi2_init;
dphi3(k) = -dphi2(k);
    A = [0,-AE*sin(phi4(k)),0,0,r7*sin(phi7(k)),0;
         0,AE*cos(phi4(k)),0,0,-r7*cos(phi7(k)),0;
         0,-r4*sin(phi4(k)),r5*sin(phi5(k)),0,0,-1;
         -1,r4*cos(phi4(k)),-r5*cos(phi5(k)),0,0,0;
         0,0,FC*sin(phi5(k)),-r6*sin(phi6(k)),0,-1;
         -1,0,-FC*cos(phi5(k)),r6*cos(phi6(k)),0,0];

    B = [ r2*sin(phi2(k))*dphi2(k)-r3*sin(phi3(k))*dphi3(k);
         -r2*cos(phi2(k))*dphi2(k)+r3*cos(phi3(k))*dphi3(k);
          r2*sin(phi2(k))*dphi2(k);
         -r2*cos(phi2(k))*dphi2(k);
         r3*sin(phi3(k))*dphi3(k);
         -r3*cos(phi3(k))*dphi3(k)];
     
    x = A\B;
    
    % save results
    dy8(k) = x(1);
    dphi4(k) = x(2);
    dphi5(k) = x(3);
    dphi6(k) = x(4);
    dphi7(k) = x(5);
    dx8(k) = x(6);
    
    % *** calculate initial values for next iteration step ***
    y8_init = y8(k)+Ts*dy8(k);
    phi4_init = phi4(k)+Ts*dphi4(k);
    phi5_init = phi5(k)+Ts*dphi5(k);
    phi6_init = phi6(k)+Ts*dphi6(k);
    phi7_init = phi7(k)+Ts*dphi7(k);
    x8_init = x8(k)+Ts*dx8(k);
    phi2_init = phi2_init+Ts*dphi2_init;


%% Compute unkown internal forces and all accelerations
x2 = r2*cos(phi2(k));
x3 = r3*cos(phi3(k));
x4 = r4*cos(phi4(k));
x5 = r5*cos(phi5(k));
x6 = r6*cos(phi6(k));
x7 = r7*cos(phi7(k));
y2 = r2*sin(phi2(k));
y3 = r3*sin(phi3(k));
y4 = r4*sin(phi4(k));
y5 = r5*sin(phi5(k));
y6 = r6*sin(phi6(k));
y7 = r7*sin(phi7(k));
AEx = x4*(r4/2-AE)/r4;
AEy = y4*(r4/2-AE)/r4;
FCx = x5*(FC-r5/2)/r5;
FCy = y5*(FC-r5/2)/r5;

% Change moments of inertia 2 and 3 to be around fixed point and not cog
J2f = J2+m2*r2^2/4;
J3f = J3+m3*r3^2/4;

A = zeros(26,26);
A(1,1) = m2*y2/2; A(1,9) = 1; A(1,13) = -1;
A(2,1) = -m2*x2/2; A(2,10) = 1; A(2,14) = -1;
A(3,1) = -J2f; A(3,14) = -x2; A(3,13) = y2;
A(4,2) = m3*y3/2; A(4,11) = 1; A(4,15) = 1; A(4,17) = 1;
A(5,2) = -m3*x3/2; A(5,12) = 1; A(5,16) = -1; A(5,18) = -1;
A(6,2) = -J3f; A(6,15) = -y3; A(6,17) = -y3; A(6,16) = -x3; A(6,18) = -x3;
A(7,1) = m4*y2; A(7,3) = y4/2*m4; A(7,13) = 1; A(7,19) = -1; A(7,21) = -1;
A(8,1) = -m4*x2; A(8,3) = -x4/2*m4; A(8,14) = 1; A(8,20) = 1; A(8,22) = 1;
A(9,3) = -J4; A(9,13) = y4/2; A(9,14) = -x4/2; A(9,19) = y4/2; A(9,20) = x4/2; A(9,21) = -AEy; A(9,22) = -AEx;
A(10,7) = -m8; A(10,25)= 1;
A(11,8) = -m8; A(11,26) = 1;
A(12,4) = m5*y5/2; A(12,7) = -m5; A(12,19) = 1; A(12,23) = 1; A(12,25) = -1;
A(13,4) = -m5*x5/2; A(13,8) = -m5; A(13,20) = -1; A(13,24) = -1; A(13,26) = -1;
A(14,4) = -J5; A(14,19) = -y5/2; A(14,20) =-x5/2; A(14,23) = -FCy; A(14,24) = -FCx; A(14,25) = -y5/2; A(14,26) = x5/2;
A(15,2) = m6*y3; A(15,5) = m6*y6/2; A(15,15) = -1; A(15,23) = -1;
A(16,2) = -m6*x3; A(16,5) = -m6*x6/2; A(16,16) = 1; A(16,24) = 1;
A(17,5) = -J6; A(17,15) = -y6/2; A(17,16) = -x6/2; A(17,23) = y6/2; A(17,24) = x6/2;
A(18,2) = m7*y3; A(18,6) = m7*y7/2; A(18,17) = -1; A(18,21) = 1;
A(19,2) = -m7*x3; A(19,6) = -m7*x7/2; A(19,18) = 1; A(19,22) = -1;
A(20,6) = -J7; A(20,17) = -y7/2; A(20,18) = -x7/2; A(20,21) = -y7/2; A(20,22) = -x7/2;
A(21,1) = y2; A(21,2) = -y3; A(21,3) = AE/r4*y4; A(21,6) = -y7;
A(22,1) = -x2; A(22,2) = x3; A(22,3) = -AE/r4*x4; A(22,6) = x7;
A(23,1) = y2; A(23,3) = y4; A(23,4) = -y5; A(23,7) = 1;
A(24,1) = -x2; A(24,3) = -x4; A(24,4) = x5; A(24,8) = 1;
A(25,2) = y3; A(25,4) = -FC/r5*y5; A(25,5) = y6; A(25,7) = 1; 
A(26,2) = -x3; A(26,4) = FC/r5*x5; A(26,5) = -x6; A(26,8) = 1;

B = [-m2*x2/2*dphi2(k)^2;
    -m2*y2/2*dphi2(k)^2;
    M2(k);
    -m3*x3/2*dphi3(k)^2;
    -m3*y3/2*dphi3(k)^2;
    -M3(k);
    -m4*x4/2*dphi4(k)^2-m4*x2*dphi2(k)^2;
    -m4*y4/2*dphi4(k)^2-m4*y2*dphi2(k)^2;
    0;
    0;
    0;
    -m5*x5/2*dphi5(k)^2;
    -m5*y5/2*dphi5(k)^2;
    0;
    -m6*x3*dphi3(k)^2-m6*x6/2*dphi6(k)^2;
    -m6*y3*dphi3(k)^2-m6*y6/2*dphi6(k)^2;
    0;
    -m7*x3*dphi3(k)^2-m7*x7/2*dphi7(k)^2;
    -m7*y3*dphi3(k)^2-m7*y7/2*dphi7(k)^2;
    0;
    -x2*dphi2(k)^2-AE/r4*x4*dphi4(k)^2+x3*dphi3(k)^2+x7*dphi7(k)^2;
    -y2*dphi2(k)^2-AE/r4*y4*dphi4(k)^2+y3*dphi3(k)^2+y7*dphi7(k)^2;
    -x2*dphi2(k)^2-x4*dphi4(k)^2+x5*dphi5(k)^2;
    -y2*dphi2(k)^2-y4*dphi4(k)^2+y5*dphi5(k)^2;
    -x3*dphi3(k)^2-x6*dphi6(k)^2+FC/r5*x5*dphi5(k)^2;
    -y3*dphi3(k)^2-y6*dphi6(k)^2+FC/r5*y5*dphi5(k)^2];

x = A\B;

ddphi2(k) = x(1);
ddphi3(k) = x(2);
ddphi4(k) = x(3);
ddphi5(k) = x(4);
ddphi6(k) = x(5);
ddphi7(k) = x(6);
ddx8(k) = x(7);
ddy8(k) = x(8);
F2x(k) = x(9);
F2y(k) = x(10);
F3x(k) = x(11);
F3y(k) = x(12);
F42x(k) = x(13);
F42y(k) = x(14);
F63x(k) = x(15);
F63y(k) = x(16);
F73x(k) = x(17);
F73y(k) = x(18);
F54x(k) = x(19);
F54y(k) = x(20);
F74x(k) = x(21);
F74y(k) = x(22);
F65x(k) = x(23);
F65y(k) = x(24);
F58x(k) = x(25);
F58y(k) = x(26);

dphi2_init = dphi2_init+Ts*ddphi2(k);
% Controls: Shaking Forces and Moments
omega2 = [0,0, dphi2(k)];
omega3 = [0,0, dphi3(k)];
omega4 = [0,0, dphi4(k)];
omega5 = [0,0, dphi5(k)];
omega6 = [0,0, dphi6(k)];
omega7 = [0,0, dphi7(k)];
alpha2 = [0, 0, ddphi2(k)];
alpha3 = [0, 0, ddphi3(k)];
alpha4 = [0, 0, ddphi4(k)];
alpha5 = [0, 0, ddphi5(k)];
alpha6 = [0, 0, ddphi6(k)];
alpha7 = [0, 0, ddphi7(k)];
% 3D model vectors
O1_cog2_vec = 1/2*[x2    y2  0];
O2_cog3_vec = 1/2*[x3  y3   0];
A_cog4_vec = 1/2*[x4   y4   0];
D_cog7_vec = 1/2*[x7 y7 0];
D_cog6_vec = 1/2*[x6 y6 0];
C_cog5_vec = 1/2*[x5 y5 0];
B_cog5_vec = -1/2*[x5,y5,0];
O1A_vec =[x2 y2 0];
AB_vec =[x4 y4 0];
AE_vec =AE/r4*AB_vec;
O2C_vec =[x8(k) y8(k) 0];
CB_vec =[x5 y5 0];
CF_vec = FC/r5*CB_vec;
DF_vec =[x6 y6 0]; 
DE_vec =[x7 y7 0]; 
O1O2_vec = [r1 0 0];
% Define arms from point 1 for shaking moment check
arm2 = O1_cog2_vec;
arm3 = O1O2_vec+O2_cog3_vec;
arm4 = 2*arm2+A_cog4_vec;
arm5 = O1O2_vec+O2C_vec+C_cog5_vec;
arm6 = O1O2_vec+2*O2_cog3_vec+D_cog6_vec;
arm7 = O1O2_vec+2*O2_cog3_vec+D_cog7_vec;
arm8 = O1O2_vec+O2C_vec;
% speed vectors
vel_2 = cross(omega2,O1_cog2_vec);
vel_3 = cross(omega3,O2_cog3_vec);
vel_C = [dx8(k),dy8(k),0];
vel_A = 2*vel_2;
vel_D = 2*vel_3;
vel_4 = vel_A+cross(omega4,A_cog4_vec);
vel_5 = vel_C+cross(omega5,C_cog5_vec);
vel_6 = vel_D+cross(omega6,D_cog6_vec);
vel_7 = vel_D + cross(omega7,D_cog7_vec);
vel_8 = vel_C;

% acceleration vectors
acc_2 =       cross(omega2,cross(omega2,O1_cog2_vec))+cross(alpha2,O1_cog2_vec);
acc_A =       cross(omega2,cross(omega2,O1A_vec))+cross(alpha2,O1A_vec);
acc_3 =       cross(omega3,cross(omega3,O2_cog3_vec))+cross(alpha3,O2_cog3_vec);
acc_4 = acc_A+cross(omega4,cross(omega4,A_cog4_vec))+cross(alpha4,A_cog4_vec);
acc_D =     2*(cross(omega3,cross(omega3,O2_cog3_vec))+cross(alpha3,O2_cog3_vec));
acc_C = [ddx8(k),ddy8(k),0];
acc_5 = acc_C+cross(omega5,cross(omega5,C_cog5_vec))+cross(alpha5,C_cog5_vec);
acc_6 = acc_D+cross(omega6,cross(omega6,D_cog6_vec))+cross(alpha6,D_cog6_vec);
acc_7 = acc_D+cross(omega7,cross(omega7,D_cog7_vec))+cross(alpha7,D_cog7_vec);
acc_8 = acc_C;
acc_2x = acc_2(:,1);
acc_2y = acc_2(:,2);
acc_3x = acc_3(:,1);
acc_3y = acc_3(:,2);
acc_4x = acc_4(:,1);
acc_4y = acc_4(:,2);
acc_5x = acc_5(:,1);
acc_5y = acc_5(:,2);
acc_6x = acc_6(:,1);
acc_6y = acc_6(:,2);
acc_7x = acc_7(:,1);
acc_7y = acc_7(:,2);
acc_8x = acc_8(:,1);
acc_8y = acc_8(:,2);

mass = [m2,m3,m4,m5,m6,m7,m8];
inertia = [J2,J3,J4,J5,J6,J7,0];
arm = [arm2;arm3;arm4;arm5;arm6;arm7;arm8];
vel = [vel_2;vel_3;vel_4;vel_5;vel_6;vel_7;vel_8];
acc = [acc_2;acc_3;acc_4;acc_5;acc_6;acc_7;acc_8];
omega = [omega2(3),omega3(3),omega4(3),omega5(3),omega6(3),omega7(3),0];
alpha = [alpha2(3),alpha3(3),alpha4(3),alpha5(3),alpha6(3),alpha7(3),0];

Fshakx(k) = -(F2x(k)+F3x(k));
Fshaky(k) = -(F2y(k)+F3y(k));
Mshak(k) = -(-M2(k)+M3(k)+r1*F3y(k));


Xtest(k) = Fshakx(k);
Ytest(k) = Fshaky(k);
Mtest(k) = Mshak(k);
for i=1:7
    Xtest(k) = Xtest(k)+mass(i)*acc(i,1);
    Ytest(k) = Ytest(k)+mass(i)*acc(i,2);
    Mr = mass(i)*cross(transpose(arm(i,:)),transpose(acc(i,:)));
    Mtest(k) = Mtest(k)+inertia(i)*alpha(i)+Mr(3);
end
end

%% Plots

if fig_forward
    figure
    subplot(4,1,1)
    plot(t,phi2);
    ylabel('phi2 [rad]')
    xlabel("time [s]")
    subplot(4,1,2)
    plot(t,dphi2+pi/5)
    ylabel('omega2+pi/5 [rad/s]')
    xlabel("time [s]")
    subplot(4,1,3)
    plot(t,ddphi2)
    ylabel('alpha2 [rad/s^2]')
    xlabel("time [s]")
    subplot(4,1,4)
    plot(t,dphi2+pi/5)
    ylabel('Absolute error')
    xlabel("time [s]")
    
    figure
    subplot(2,2,1)
    plot(t,Fshakx);
    ylabel("Shaking Force X-direction [N]")
    xlabel("time [s]")
    subplot(2,2,2)
    plot(t,Xtest)
    ylabel("Absolute error [N]")
    xlabel("time [s]")
    subplot(2,2,3)
    plot(t,Fshaky);
    ylabel("Shaking Force Y-direction [N]")
    xlabel("time [s]")
    subplot(2,2,4)
    plot(t,Ytest)
    ylabel("Absolute error [N]")
    xlabel("time [s]")
    
    figure
    subplot(2,1,1)
    plot(t,Mshak);
    ylabel("Shaking Moment [Nm]")
    xlabel("time [s]")
    subplot(2,1,2)
    plot(t,Mtest)
    ylabel("Absolute error [Nm]")
    xlabel("time [s]")
end