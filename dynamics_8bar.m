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


function [M2t2,F2t2x,F2t2y,Ft2x,Ft2y,F3t3x,F3t3y,M3t3,Ft3y,Ft3x,F42x,F42y,F73x,F73y,F63y,F63x,F74x,F74y,F54x,F54y,F65x,F65y,M2,M3] = ...
dynamics_8bar(phi2,phi3,phi4,phi5,phi6,phi7,x8,y8,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dx8,dy8,...
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
x7 = r6*cos(phi7);
y2 = r2*sin(phi2);
y3 = r3*sin(phi3);
y4 = r4*sin(phi4);
y5 = r5*sin(phi5);
y6 = r6*sin(phi6);
y7 = r7*sin(phi7);
AEx = x4*AE/r4;
AEy = y4*AE/r4;
FCx = x5*FC/r5;
FCy = y5*FC/r5;
xt2 = r1/2;
xt3 = -xt2;
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
O1A_vec =[x2 y2 zeros(size(phi2))];
AB_vec =[x4 y4 zeros(size(phi2))];
AE_vec =AE/r4*AB_vec;
O2C_vec =[x8 y8 zeros(size(phi2))];
CB_vec =[x5 y5 zeros(size(phi2))];
CF_vec = FC/r5*CB_vec;
DF_vec =[x6 y6 zeros(size(phi2))]; 
DE_vec =[x7 y7 zeros(size(phi2))]; 
O1O2_vec = [r1*ones(size(phi2)) zeros(size(phi2)) zeros(size(phi2))];

% acceleration vectors
acc_2 =       cross(omega2,cross(omega2,O1_cog2_vec))+cross(alpha2,O1_cog2_vec);
acc_A =       cross(omega2,cross(omega2,O1A_vec))+cross(alpha2,O1A_vec);
acc_3 =       cross(omega3,cross(omega3,O2_cog3_vec))+cross(alpha3,O2_cog3_vec);
acc_4 = acc_A+cross(omega4,cross(omega4,A_cog4_vec))+cross(alpha4,A_cog4_vec);
acc_D =     2*cross(omega3,cross(omega3,O2_cog3_vec))+cross(alpha3,O2_cog3_vec);
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


% **********************
% *** force analysis ***
% **********************

% allocate matrices for force (F) and moment (M)
M2t2 = zeros(size(phi2));
F2t2x = zeros(size(phi2));
F2t2y = zeros(size(phi2));
Ft2y = zeros(size(phi2));
Ft2x = zeros(size(phi2));
F3t3x = zeros(size(phi2));
F3t3y = zeros(size(phi2));
M3t3 = zeros(size(phi2));
Ft3y = zeros(size(phi2));
Ft3x = zeros(size(phi2));
F42y = zeros(size(phi2));
F42x = zeros(size(phi2));
F73x  = zeros(size(phi2));
F73y = zeros(size(phi2));
F63y = zeros(size(phi2));
F63x = zeros(size(phi2));
F74x = zeros(size(phi2));
F74y = zeros(size(phi2));
F54x = zeros(size(phi2));
F54y = zeros(size(phi2));
F65x = zeros(size(phi2));
F65y = zeros(size(phi2));
% F85x = zeros(size(phi2));
% F85y = zeros(size(phi2));
M2 = zeros(size(phi2));
M3 = zeros(size(phi2));
M = M2+M3;

% calculate dynamics for each time step
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
  A = [ 0           -1           0            0            0            0            0           0           0      0      0      1    0    0    0    0    0    0    0    0    0    0    0    0  ;
        0           0            -1           0            0            0            0           0           0     0       1      0    0    0    0    0    0    0    0    0    0    0    0    0  ;
        0           0            0            0            0            -1           0           0           0      0      0      0    1    0    0    1    0    0    0    0    0    0    0    0  ;
        0           0            0            0            0            0           -1           0           0      0      0      0    0    1    1    0    0    0    0    0    0    0    0    0   ;
        0           0            0            0            0            0            0           0           0      0      0     -1    0    0    0    0    1    0    1    0    0    0    0    0  ;
        0           0            0            0            0            0            0           0           0      0     -1      0    0    0    0    0    0    1    0    1    0    0    0    0  ;
        0           0            0            0            0            0            0           0           0      0      0      0    0    0    0    0    0    0    -1   0   1    0    0   0   ;  
        0           0            0            0            0            0            0           0           0      0      0      0    0    0    0    0    0     0    0   -1   0    1    0   0  ;
        0           0            0            0            0            0            0           0           0      0      0      0    0    0    0    -1    0    0    0    0   -1   0    0    0  ;
        0           0            0            0            0            0            0           0           0      0      0      0    0    0    -1    0    0    0    0    0    0   -1    0    0 ;   
        0           0            0            0            0            0            0           0           0      0      0      0    -1    0    0    0    -1    0    0    0   0   0    0    0  ;
        0           0            0            0            0            0            0           0           0      0      0      0    0    -1    0    0    0    -1    0    0   0   0    0    0  ;
        0           1            0            0            1            0            0           0           0      0      0      0    0    0    0     0    0    0    0    0    0   0    0    0  ;
        0           0            1            1            0            0            0           0           0      0      0      0    0    0    0     0    0    0    0    0    0   0    0    0  ;
        0           0            0            0            0            1            0           0           0      1      0      0    0    0    0     0    0    0    0    0    0   0    0    0  ;  
        0           0            0            0            0            0            1           0           1      0      0      0    0    0    0      0    0   0    0    0    0    0    0   0  ;    
        1           0            0            0            0            0            0           0           0      0      0      0    0    0    0     0    0    0    0    0    0   0    -1    0 ;
        0           0            0            0            0            0            0           -1           0      0      0      0    0    0    0     0    0    0    0    0    0   0    0    1 ;
        1           0            0            0            0            0            0           0           0      0     -x2(k) y2(k)    0    0    0     0    0    0    0    0    0   0    0    0   ;
        0           0            0            0            0            0            0           -1           0      0      0      0   -y3(k)  -x3(k) -x3(k)  -y3(k)  0    0    0    0   0   0    0    0   ;
        0           0            0            0            0            0            0           0           0      0   x4(k)/2    y4(k)/2    0    0    0    0   -(AEy(k)-y4(k)/2)     -(AEx(k)-x4(k)/2)    -y4(k)/2    x4(k)/2   0   0    0    0  ;
        0           0            0            0            0            0            0           0           0      0      0      0    0    0    0    0    0    0    y5(k)/2  x5(k)/2  -(FCy(k)-y5(k)/2) (FCx(k)-x5(k)/2)   0   0 ;
        0           0            0            0            0            0            0           0           0      0      0      0    0    0   x6(k)/2  -y6(k)/2    0    0    0    0  -y6(k)/2   x6(k)/2    0    0   ;
        0           0            0            0            0            0            0           0           0      0      0      0   -y7(k)/2  x7(k)/2    0    0    y7(k)/2    x7(k)/2    0    0   0   0    0    0   ];
  [R, bj] = rref(A')     
  B = [ m2*acc_2x(k);
        m2*acc_2y(k);
        m3*acc_3x(k);
        m3*acc_3y(k);
        m4*acc_4x(k);
        m4*acc_4y(k);
        m5*acc_5x(k);
        m5*acc_5y(k);
        m6*acc_6x(k);
        m6*acc_6y(k);
        m7*acc_7x(k);
        m7*acc_7y(k);
         0;
         0;
         0;
         0;
        Jt2*ddphi2(k);
        Jt3*ddphi3(k);
        J2*ddphi2(k);
        J3*ddphi3(k);
        J4*ddphi4(k);
        J5*ddphi5(k);
        J6*ddphi6(k);
        J7*ddphi7(k)];
        
    
    x = A\B;
    
    % save results
M2t2(k) = x(1);
F2t2x(k) = x(2);
F2t2y(k) = x(3);
Ft2y(k) =x(4);
Ft2x(k) =x(5);
F3t3x(k)= x(6);
F3t3y(k) = x(7);
M3t3(k) =x(8);
Ft3y(k) =x(9);
Ft3x(k) = x(10);
F42y(k) =x(11);
F42x(k) = x(12);
F73x(k) =x(13);
F73y(k) = x(14);
F63y(k) = x(15);
F63x(k) = x(16);
F74x(k) = x(17);
F74y(k) = x(18);
F54x(k) = x(19);
F54y(k) = x(20);
F65x(k) = x(21);
F65y(k) =x(22);
M2(k) = x(23);
M3(k) = x(24);
M(k) = M2(k)+M3(k);
end



% **********************
% *** plot figures ***
% **********************

if fig_dyn_8bar
    
    figure
    subplot(221)
    plot(t,M),grid
    xlabel('time')
    ylabel('Maandrijf')
    axis tight
    subplot(222)
    plot(F73x,F73y),grid
    xlabel('F73x [N]')
    ylabel('F73y [N]')
    axis tight
%     subplot(223)
%     plot(F_R_x,F_R_y),grid
%     xlabel('F_R_x [N]')
%     ylabel('F_R_y [N]')
%     axis tight
%     subplot(224)
%     plot(F_S_x,F_S_y),grid
%     xlabel('F_S_x [N]')
%     ylabel('F_S_y [N]')
%     axis tight
    
    figure
    plot(t,M)
    ylabel('M [N-m]')
    xlabel('t [s]')
    
end


