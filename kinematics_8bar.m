%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phi2,phi3,phi4,phi5,phi6,x7,dphi2,dphi3,dphi4,dphi5,dphi6,dx7,ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddx7] = kinematics_8bar(r0,r1,r2,r5,r6,AE,FC,BF,EB,phi0,phi1,dphi1,ddphi1,phi2_init,phi3_init,phi4_init,phi5_init,phi6_init,x7_init,t,fig_kin_8bar)

% allocation of the result vectors (this results in better performance because we don't have to reallocate and
% copy the vector each time we add an element.
phi2 = zeros(size(t));
phi3 = zeros(size(t));
phi4 = zeros(size(t));
phi5 = zeros(size(t));
phi6 = zeros(size(t));
x7 = zeros(size(t));
dphi2 = zeros(size(t));
dphi3 = zeros(size(t));
dphi4 = zeros(size(t));
dphi5 = zeros(size(t));
dphi6 = zeros(size(t));
dx7 = zeros(size(t));
ddphi2 = zeros(size(t));
ddphi3 = zeros(size(t));
ddphi4 = zeros(size(t));
ddphi5 = zeros(size(t));
ddphi6 = zeros(size(t));
ddx7 = zeros(size(t));
Control1 = zeros(size(t));
Control2 = zeros(size(t));
Control3 = zeros(size(t));
Control4 = zeros(size(t));
Control5 = zeros(size(t));
Control6 = zeros(size(t));

numdphi2 = zeros(size(t));
control_numeric_dphi2 =  zeros(size(t));
numdphi3 = zeros(size(t));
control_numeric_dphi3 =  zeros(size(t));
numddphi2 = zeros(size(t));
control_numeric_ddphi2 =  zeros(size(t));
numddphi3 = zeros(size(t));
control_numeric_ddphi3 =  zeros(size(t));


% fsolve options (help fsolve, help optimset)
optim_options = optimset('Display','off');

% *** loop over positions ***
Ts = t(2) - t(1);      % timestep
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
    
    % *** position analysis ***
    
    % fsolve solves the non-linear set of equations
    % loop closure equations: see loop_closure_eqs.m
    % argument loop_closure_eqs: file containing closure equations
    % argument [..]': initial values of unknown angles
    % argument optim options: parameters for fsolve
    % argument phi1(k): input angle phi1 for which we want to calculate the unknown angles
    % argument a1 ... phi0: constants
    % return value x: solution for the unknown angles
    % return exitflag: indicates convergence of algorithm
    [x, fval, exitflag]=fsolve('loop_closure_eqs_8bar',[phi2_init,phi3_init phi4_init,phi5_init,phi6_init,x7_init]',optim_options,phi1(k),r0,r1,r2,r5,r6,AE,FC,BF,EB,phi0);
    if (exitflag ~= 1)
        display 'The fsolve exit flag was not 1, probably no convergence!'
        exitflag
        k
    else
    end
    
    % save results of fsolve
    phi2(k) = x(1);
    phi3(k)=x(2);
    phi4(k)=x(3);
    phi5(k) = x(4);
    phi6(k) = x(5);
    x7(k) = x(6);
    
    
    % *** velocity analysis ***
    
    A = [r2*sin(phi2(k)),-AE*sin(phi3(k)),0,0,r6*sin(phi6(k)),0;
         -r2*cos(phi2(k)),AE*cos(phi3(k)),0,0,-r6*cos(phi6(k)),0;
         -r2*sin(phi2(k)),0,-FC*sin(phi4(k)),-r5*sin(phi5(k)),0,-1;
         r2*cos(phi2(k)),0,FC*cos(phi4(k)),r5*cos(phi5(k)),0,0;
         0,-EB*sin(phi3(k)),-BF*sin(phi4(k)),r5*sin(phi5(k)),-r6*sin(phi6(k)),0;
         0,EB*cos(phi3(k)),BF*cos(phi4(k)),-r5*cos(phi5(k)),r6*cos(phi6(k)),0];
     
    B = [ r1*sin(phi1(k))*dphi1(k);
         -r1*cos(phi1(k))*dphi1(k);
         0;
         0;
         0;
         0];
     
    x = A\B;
    
    % save results
    dphi2(k) = x(1);
    dphi3(k) = x(2);
    dphi4(k) = x(3);
    dphi5(k) = x(4);
    dphi6(k) = x(5);
    dx7(k) = x(6);
    
    
    % *** acceleration analysis ***
    
    A = [r2*sin(phi2(k)),-AE*sin(phi3(k)),0,0,r6*sin(phi6(k)),0;
        -r2*cos(phi2(k)),AE*cos(phi3(k)),0,0,-r6*cos(phi6(k)),0;
        -r2*sin(phi2(k)),0,-FC*sin(phi4(k)),-r5*sin(phi5(k)),0,-1;
        r2*cos(phi2(k)),0,FC*cos(phi4(k)),r5*cos(phi5(k)),0,0;
        0,-EB*sin(phi3(k)),-BF*sin(phi4(k)),r5*sin(phi5(k)),-r6*sin(phi6(k)),0;
        0,EB*cos(phi3(k)),BF*cos(phi4(k)),-r5*cos(phi5(k)),r6*cos(phi6(k)),0];
     
    B = [AE*dphi3(k)^2*cos(phi3(k))+r1*sin(phi1(k))*ddphi1(k)-r6*cos(phi6(k))*dphi6(k)^2-r2*cos(phi2(k))*dphi2(k)^2+r1*cos(phi1(k))*dphi1(k)^2;
        AE*dphi3(k)^2*sin(phi3(k))-r1*cos(phi1(k))*ddphi1(k)-r6*sin(phi6(k))*dphi6(k)^2-r2*sin(phi2(k))*dphi2(k)^2+r1*sin(phi1(k))*dphi1(k)^2;
        r2*cos(phi2(k))*dphi2(k)^2+r5*cos(phi5(k))*dphi5(k)^2+FC*cos(phi4(k))*dphi4(k)^2;
        r2*sin(phi2(k))*dphi2(k)^2+r5*sin(phi5(k))*dphi5(k)^2+FC*sin(phi4(k))*dphi4(k)^2;
        r6*cos(phi6(k))*dphi6(k)^2+EB*cos(phi3(k))*dphi3(k)^2+BF*cos(phi4(k))*dphi4(k)^2-r5*cos(phi5(k))*dphi5(k)^2;
        r6*sin(phi6(k))*dphi6(k)^2+EB*sin(phi3(k))*dphi3(k)^2+BF*sin(phi4(k))*dphi4(k)^2-r5*sin(phi5(k))*dphi5(k)^2];
    
    x = A\B;
    % save results
    ddphi2(k) = x(1);
    ddphi3(k) = x(2);
    ddphi4(k) = x(3);
    ddphi5(k) = x(4);
    ddphi6(k) = x(5);
    ddx7(k) = x(6);
    
    
    % *** calculate initial values for next iteration step ***
    phi2_init = phi2(k);%+Ts*dphi2(k);
    phi3_init = phi3(k);%+Ts*dphi3(k);
    phi4_init = phi4(k);%+Ts*dphi4(k);
    phi5_init = phi5(k);%+Ts*dphi5(k);
    phi6_init = phi6(k);%+Ts*dphi6(k);
    x7_init = x7(k);%+Ts*dx7(k);
    
    
    % Controleren door kinematische lus
    Control1(k) = (r1*exp(j*phi1(k))+AE*exp(j*phi3(k)))-(r2*exp(j*phi2(k))+r6*exp(j*phi6(k))+r0);
    Control1(k) = abs(Control1(k));
    Control2(k) = (r2*exp(j*phi2(k))+r5*exp(j*phi5(k))+FC*exp(j*phi4(k)))-x7(k);
    Control2(k) = abs(Control2(k));
    
    % Controleren door numerieke afgeleide
    if k>2 && k<t_size
        numdphi2(k) = (phi2(k+1)-phi2(k+1))/(2*Ts);
    else
        numdphi2(k) = 0;
    end
    control_numeric_dphi2(k) = numdphi2(k)-dphi2(k);
end % loop over positions



% *** create movie ***

% point  O1 = fixed
O1 = 0;
% point O2 = fixed
O2 = r0;
% define which positions we want as frames in our movie
frames = 40;    % number of frames in movie
delta = floor(t_size/frames); % time between frames
index_vec = [1:delta:t_size]';

% Create a window large enough for the whole mechanisme in all positions, to prevent scrolling.
% This is done by plotting a diagonal from (x_left, y_bottom) to (x_right, y_top), setting the
% axes equal and saving the axes into "movie_axes", so that "movie_axes" can be used for further
% plots.
x_left = -1.5*r1;
y_bottom = -0.2;
x_right = r0+5*AE;
y_top = r1+2*AE;

figure(10)
hold on
plot([x_left, x_right], [y_bottom, y_top]);
axis equal;
movie_axes = axis;   %save current axes into movie_axes

% draw and save movie frame
for m=1:length(index_vec)
    index = index_vec(m);
    A = O1 + r1 * exp(j*phi1(index));
    E = A + AE * exp(j*phi3(index));
    B = E + EB * exp(j*phi3(index));
    D = O2 + r2*exp(j*phi2(index));
    %F = D + r5*exp(j*phi5(index));
    F = B + BF*exp(j*phi4(index));
    C = F + FC*exp(j*phi4(index));
    
    loop1 = [O1 A E B F C];
    loop2 = [O2 D E B F D O2];
    
    figure(10)
    clf
    hold on
    plot(real(loop1),imag(loop1),'-ro',real(loop2),imag(loop2),'-ro')
    
    axis(movie_axes);     % set axes as in movie_axes
    Movie(m) = getframe;  % save frame to a variable Film
end

% save movie
save fourbar_movie Movie
close(10)


% *** plot figures ***

if fig_kin_8bar
    
    %plot assembly at a certain timestep 
    index = 1; %select 1st timestep
    O1 = 0;
    O2 = r0;
    A = O1 + r1 * exp(j*phi1(index));
    E = A + AE * exp(j*phi3(index));
    B = E + EB * exp(j*phi3(index));
    D = O2 + r2*exp(j*phi2(index));
   %F = D + r5*exp(j*phi5(index));
    F = B + BF*exp(j*phi4(index));
    C = F + FC*exp(j*phi4(index));
    
    figure
    assembly1=[O1 A E B F C];
    assembly2= [O2 D E B F D O2];
    plot(real(assembly1),imag(assembly1),'ro-',real(assembly2),imag(assembly2),'ro-')
    xlabel('[m]')
    ylabel('[m]')
    title('assembly')
    axis equal
    
    figure
    subplot(611)
    plot(t,phi1)
    ylabel('\phi_1 [rad]')
    subplot(612)
    plot(t,phi2)
    ylabel('\phi_2 [rad]')
    subplot(613)
    plot(t,x7)
    ylabel('x_7 [m]')
    xlabel('t [s]')
    subplot(614)
    plot(t,phi4)
    ylabel('\phi_4 [rad]')
    subplot(615)
    plot(t,phi5)
    ylabel('\phi_5 [rad]')
    subplot(616)
    plot(t,phi6)
    ylabel('\phi_6 [rad]')
    
    figure
    subplot(311)
    plot(t,dphi1)
    ylabel('d\phi_1 [rad/s]')
    subplot(312)
    plot(t,dphi3)
    ylabel('d\phi_3 [rad/s]')
    subplot(313)
    plot(t,dx7)
    ylabel('dx_7 [m/s]')
    xlabel('t [s]')
    
    figure
    subplot(311)
    plot(t,ddphi1)
    ylabel('dd\phi_1 [rad/s^2]')
    subplot(312)
    plot(t,ddphi3)
    ylabel('dd\phi_3 [rad/s^2]')
    subplot(313)
    plot(t,ddx7)
    ylabel('ddx_7 [m/s^2]')
    xlabel('t [s]')
    
    % Plot differnce of position using 2 different kinematic chains
    figure
    subplot(211)
    plot(t,Control1)
    ylabel('xerror_E')
    xlabel('t [s]')
    subplot(212)
    plot(t,Control2)
    ylabel('xerror_C')
    xlabel('t [s]')
    
    % Plot difference of angle velocity computed numerically and exactly
    figure
    plot(t,control_numeric_dphi2)
    ylabel('phi2error')
    xlabel('t[s]')
    
end


