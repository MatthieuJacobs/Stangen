%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test git
% Kinematica en werkuigendynamica.
%
%LOOP CLOSURE EQUATIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function F=loop_closure_eqs_8bar(phi_init,phi2,r1,r2,r3,r4,r5,r6,r7,AE,FC)

% first argument: the initial values of the unknown angles and unknown
% distance x8
% argument phi1: input angle phi1 for which we want to calculate the
% unknowns
% arguments: constants

% copy initial values of unknown angles and x8
phi3 = phi_init(1);
phi4=phi_init(2);
phi5=phi_init(3);
phi6 = phi_init(4);
phi7 = phi_init(5);
x8 = phi_init(6);

% loop closure equations:
F(1)=r2*cos(phi2)+AE*cos(phi4)-r7*cos(phi7)-r3*cos(phi3)-r1;
F(2)=r2*sin(phi2)+AE*sin(phi4)-r7*sin(phi7)-r3*sin(phi3);
F(3)=r2*cos(phi2)+r4*cos(phi4)-r1-x8-r5*cos(phi5);
F(4)=r2*sin(phi2)+r4*sin(phi4)-r5*sin(phi5);
F(5)=r3*cos(phi3)+r6*cos(phi6)-x8-FC*cos(phi5);
F(6)=r3*sin(phi3)+r6*sin(phi6)-FC*sin(phi5);
