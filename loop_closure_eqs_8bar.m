%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
%LOOP CLOSURE EQUATIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function F=loop_closure_eqs(phi_init,phi1,r0,r1,r2,r5,r6,AE,FC,BF,EB)

% first argument: the initial values of the unknown angles and unknown
% distance x7
% argument phi1: input angle phi1 for which we want to calculate the
% unknowns
% arguments: constants

% copy initial values of unknown angles and x7
phi2 = phi_init(1);
phi3=phi_init(2);
phi4=phi_init(3);
phi5 = phi_init(4);
phi6 = phi_init(5);
x7 = phi_init(6);

% loop closure equations:
F(1)=AE*cos(phi3)-r6*cos(phi6)-r2*cos(phi2)-r0+r1*cos(phi1);
F(2)=AE*sin(phi3)+r1*sin(phi1)-r6*sin(phi6)-r2*sin(phi2);
F(3)=r2*cos(phi2)+r5*cos(phi5)+FC*cos(phi4)-x7;
F(4)=r2*sin(phi2)+r5*sin(phi5)+FC*sin(phi4);
F(5)=r6*cos(phi6)+EB*cos(phi3)+BF*cos(phi4)-r5*cos(phi5);
F(6)=r6*sin(phi6)+EB*sin(phi3)+BF*sin(phi4)-r5*sin(phi5);

