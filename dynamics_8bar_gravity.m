%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Voorbeeldanalyse van een vierstangenmechanisme.
%
% Bram Demeulenaere <bram.demeulenaere@mech.kuleuven.be>
% Maarten De Munck <maarten.demunck@mech.kuleuven.be>
% Johan Rutgeerts <johan.rutgeerts@mech.kuleuven.be>
% Wim Meeussen <wim.meeussen@mech.kuleuven.be>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [F2x,F2y,F3x,F3y,F42x,F42y,F63x,F63y,F73x,F73y,F54x,F54y,F74x,F74y,F65x,F65y,M2,M3] = ...
dynamics_8bar_gravity(phi2,phi3,phi4,phi5,phi6,phi7,x8,y8,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dx8,dy8,...
    ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddx8,ddy8,...
    r2,r3,r4,r5,r6,r7,AE,FC,r1, ...
    m2,m3,m4,m5,m6,m7,m8,...
    J2,J3,J4,J5,J6,J7,Jt2,Jt3,t,fig_dyn_8bar)


% a lot of definitions to make the matrix A and B a bit clear.
% skip the definitions for now (move down to "force analysis")
% and check them when you need them.


% cogi_P_x, cogn_P_y = vector from the centre of gravity of bar i to point P
x2 = r2*cos(phi2);
x3 = r3*cos(phi3);
x4 = r4*cos(phi4);
x5 = r5*cos(phi5);
x6 = r6*cos(phi6);
x7 = r7*cos(phi7);
y2 = r2*sin(phi2);
y3 = r3*sin(phi3);
y4 = r4*sin(phi4);
y5 = r5*sin(phi5);
y6 = r6*sin(phi6);
y7 = r7*sin(phi7);
AEx = x4*(r4/2-AE)/r4;
AEy = y4*(r4/2-AE)/r4;
FCx = x5*(FC-r5/2)/r5;
FCy = y5*(FC-r5/2)/r5;
xt2 = r1/2;
xt3 = -xt2;
% Change moments of inertia 2 and 3 to be around fixed point and not cog
J2f = J2+m2*r2^2/4;
J3f = J3+m3*r3^2/4;
% 3D omega (dphi) and alpha (ddphi) vectors)
omega2 = [zeros(size(phi2)) zeros(size(phi2)) dphi2];
omega3 = [zeros(size(phi2)) zeros(size(phi2)) dphi3];
omega4 = [zeros(size(phi2)) zeros(size(phi2)) dphi4];
omega5 = [zeros(size(phi2)) zeros(size(phi2)) dphi5];
omega6 = [zeros(size(phi2)) zeros(size(phi2)) dphi6];
omega7 = [zeros(size(phi2)) zeros(size(phi2)) dphi7];
alpha2 = [zeros(size(phi2)) zeros(size(phi2)) ddphi2];
alpha3 = [zeros(size(phi2)) zeros(size(phi2)) ddphi3];
alpha4 = [zeros(size(phi2)) zeros(size(phi2)) ddphi4];
alpha5 = [zeros(size(phi2)) zeros(size(phi2)) ddphi5];
alpha6 = [zeros(size(phi2)) zeros(size(phi2)) ddphi6];
alpha7 = [zeros(size(phi2)) zeros(size(phi2)) ddphi7];
% 3D model vectors
O1_cog2_vec = 1/2*[x2    y2    zeros(size(phi2))];
O2_cog3_vec = 1/2*[x3  y3    zeros(size(phi2))];
A_cog4_vec = 1/2*[x4   y4    zeros(size(phi2))];
D_cog7_vec = 1/2*[x7 y7 zeros(size(phi2))];
D_cog6_vec = 1/2*[x6 y6 zeros(size(phi2))];
C_cog5_vec = 1/2*[x5 y5 zeros(size(phi2))];
B_cog5_vec = -1/2*[x5,y5,zeros(size(phi2))];
O1A_vec =[x2 y2 zeros(size(phi2))];
AB_vec =[x4 y4 zeros(size(phi2))];
AE_vec =AE/r4*AB_vec;
O2C_vec =[x8 y8 zeros(size(phi2))];
CB_vec =[x5 y5 zeros(size(phi2))];
CF_vec = FC/r5*CB_vec;
DF_vec =[x6 y6 zeros(size(phi2))]; 
DE_vec =[x7 y7 zeros(size(phi2))]; 
O1O2_vec = [r1*ones(size(phi2)) zeros(size(phi2)) zeros(size(phi2))];
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
vel_C = [dx8,dy8,zeros(size(phi2))];
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
acc_C = [ddx8,ddy8,zeros(size(phi2))];
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
% Add gravity
g = 9.81;
% Testing
Xtest = zeros(size(phi2));
Ytest = zeros(size(phi2));
Mtest = zeros(size(phi2));
Fshakx = zeros(size(phi2));
Fshaky = zeros(size(phi2));
Mshak = zeros(size(phi2));
% **********************
% *** force analysis ***
% **********************

% allocate matrices for force (F) and moment (M)

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
M2 = zeros(size(phi2));
M3 = zeros(size(phi2));
F58x = zeros(size(phi2));
F58y =  zeros(size(phi2));
M = -M2+M3;

% calculate dynamics for each time step
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
 A = zeros(20,20);
 A(1,1) = 1; A(1,5) = -1;
 A(2,2) = 1; A(2,6) = -1;
 A(3,3) = 1; A(3,7) = 1; A(3,9) = 1;
 A(4,4) = 1; A(4,8) = -1; A(4,10) = -1;
 A(5,5) = 1; A(5,13) = -1; A(5,11) = -1;
 A(6,6) = 1; A(6,14) = 1; A(6,12) = 1;
 A(7,15) = 1; A(7,11) = 1; A(7,19) = -1;
 A(8,16) = -1; A(8,12) = -1; A(8,20) = -1;
 A(9,7) = -1; A(9,15) = -1;
 A(10,8) = 1; A(10,16) = 1;
 A(11,9) = -1; A(11,13) = 1;
 A(12,10) = 1; A(12,14) = -1;
 A(13,17) = -1; A(13,6) = -x2(k); A(13,5) = y2(k);
 A(14,18) = 1; A(14,10) = -x3(k); A(14,8) = -x3(k); A(14,9) = -y3(k); A(14,7) = -y3(k);
 A(15,14) = -AEx(k); A(15,13) = -AEy(k); A(15,12) = x4(k)/2; A(15,11) = y4(k)/2; A(15,5) = y4(k)/2; A(15,6) = -x4(k)/2;
 A(16,12) = -x5(k)/2; A(16,11) = -y5(k)/2; A(16,16) = -FCx(k); A(16,15) = -FCy(k); A(16,20) = x5(k)/2; A(16,19) = -y5(k)/2;
 A(17,15) = y6(k)/2; A(17,16) = x6(k)/2; A(17,7) = -y6(k)/2; A(17,8) = -x6(k)/2;
 A(18,14) = -x7(k)/2; A(18,10) = -x7(k)/2; A(18,13) = -y7(k)/2; A(18,9) = -y7(k)/2;
 A(19,19) = 1;
 A(20,20) = 1;
  %[R, bj] = rref(A')     
  B = [ m2*acc_2x(k);
        m2*acc_2y(k)+m2*g;
        m3*acc_3x(k);
        m3*acc_3y(k)+m3*g;
        m4*acc_4x(k);
        m4*acc_4y(k)+m4*g;
        m5*acc_5x(k);
        m5*acc_5y(k)+m5*g;
        m6*acc_6x(k);
        m6*acc_6y(k)+m6*g;
        m7*acc_7x(k);
        m7*acc_7y(k)+m7*g;
        J2f*ddphi2(k)+m2*g*x2(k)/2;
        J3f*ddphi3(k)+m3*g*x3(k)/2;
        J4*ddphi4(k);
        J5*ddphi5(k);
        J6*ddphi6(k);
        J7*ddphi7(k);
        m8*acc_8x(k);
        m8*acc_8y(k)+m8*g];
        
    
    x = A\B;
    
    % save results
F2x(k) =x(1);
F2y(k) =x(2);
F3x(k)= x(3);
F3y(k) = x(4);
F42x(k) =x(5);
F42y(k) = x(6);
F63x(k) =x(7);
F63y(k) = x(8);
F73x(k) = x(9);
F73y(k) = x(10);
F54x(k) = x(11);
F54y(k) = x(12);
F74x(k) = x(13);
F74y(k) = x(14);
F65x(k) = x(15);
F65y(k) =x(16);
M2(k) = x(17);
M3(k) = x(18);
F58x(k) = x(19);
F58y(k) = x(20);
M(k) = -M2(k)-M3(k);

%%%%% Dynamic checks %%%%%%%%%%%%

mass = [m2,m3,m4,m5,m6,m7,m8];
inertia = [J2,J3,J4,J5,J6,J7,0];
arm = [arm2(k,:);arm3(k,:);arm4(k,:);arm5(k,:);arm6(k,:);arm7(k,:);arm8(k,:)];
acc = [acc_2(k,:);acc_3(k,:);acc_4(k,:);acc_5(k,:);acc_6(k,:);acc_7(k,:);acc_8(k,:)];
alpha = [alpha2(k,3),alpha3(k,3),alpha4(k,3),alpha5(k,3),alpha6(k,3),alpha7(k,3),0];


% Method : Shaking Forces and Moments
Fshakx(k) = -(F2x(k)+F3x(k));
Fshaky(k) = -(F2y(k)+F3y(k));
Mshak(k) = -(-M2(k)+M3(k)+r1*F3y(k));


Xtest(k) = Fshakx(k);
Ytest(k) = Fshaky(k);
Mtest(k) = Mshak(k);
for i=1:7
    acc(i,2) = acc(i,2)+g;
    Xtest(k) = Xtest(k)+mass(i)*acc(i,1);
    Ytest(k) = Ytest(k)+mass(i)*acc(i,2);
    Mr = mass(i)*cross(transpose(arm(i,:)),transpose(acc(i,:)));
    Mtest(k) = Mtest(k)+inertia(i)*alpha(i)+Mr(3);
end
end


% **********************
% *** plot figures ***
% **********************

if fig_dyn_8bar
    
    figure
    plot(t,M),grid
    xlabel("time [s]")
    ylabel('Driving Moment [Nm]')
    axis tight
    figure
    subplot(2,1,1)
    plot(F73x,F73y),grid
    xlabel('F73x [N]')
    ylabel('F73y [N]')
    axis tight
    subplot(2,1,2)
    plot(F65x,F65y),grid
    ylabel('F65x')
    xlabel("time [s]")
    axis tight
    
    
%     figure
%     subplot(2,1,1)
%     plot(t,M3)
%     ylabel('M3 [Nm]')
%     xlabel("time [s]")
%     subplot(2,1,2)
%     plot(t,M2);
%     ylabel('M2 [Nm]')
%     xlabel("time [s]")
% 
%     figure
%     subplot(2,2,1)
%     plot(t,Fshakx);
%     ylabel("Shaking Force X-direction [N]")
%     xlabel("time [s]")
%     subplot(2,2,2)
%     plot(t,Xtest)
%     ylabel("Absolute error [N]")
%     xlabel("time [s]")
%     subplot(2,2,3)
%     plot(t,Fshaky);
%     ylabel("Shaking Force Y-direction [N]")
%     xlabel("time [s]")
%     subplot(2,2,4)
%     plot(t,Ytest)
%     ylabel("Absolute error [N]")
%     xlabel("time [s]")
%     
%     figure
%     subplot(2,1,1)
%     plot(t,Mshak);
%     ylabel('Shaking Moment [Nm]')
%     xlabel("time [s]")
%     subplot(2,1,2)
%     plot(t,Mtest)
%     ylabel('Absolute error [Nm]')
%     xlabel("time [s]")
    
    
end


