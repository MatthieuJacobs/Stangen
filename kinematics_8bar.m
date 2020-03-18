%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y8,phi4,phi5,phi6,phi7,x8,dy8,dphi4,dphi5,dphi6,dphi7,dx8,ddy8,ddphi4,ddphi5,ddphi6,ddphi7,ddx8] = kinematics_8bar(r1,r2,r3,r4,r5,r6,r7,AE,FC,phi2,dphi2,ddphi2,y8_init,phi4_init,phi5_init,phi6_init,phi7_init,x8_init,t,fig_kin_8bar)

% allocation of the result vectors (this results in better performance because we don't have to reallocate and
% copy the vector each time we add an element.
y8 = zeros(size(t));
phi4 = zeros(size(t));
phi5 = zeros(size(t));
phi6 = zeros(size(t));
phi7 = zeros(size(t));
x8 = zeros(size(t));
dy8 = zeros(size(t));
dphi4 = zeros(size(t));
dphi5 = zeros(size(t));
dphi6 = zeros(size(t));
dphi7 = zeros(size(t));
dx8 = zeros(size(t));
ddy8 = zeros(size(t));
ddphi4 = zeros(size(t));
ddphi5 = zeros(size(t));
ddphi6 = zeros(size(t));
ddphi7 = zeros(size(t));
ddx8 = zeros(size(t));
Control1 = zeros(size(t));
Control2 = zeros(size(t));
Control3 = zeros(size(t));
Control4 = zeros(size(t));
Control5 = zeros(size(t));
Control6 = zeros(size(t));

numdphi3 = zeros(size(t));
control_numeric_dphi3 =  zeros(size(t));
numdphi4 = zeros(size(t));
control_numeric_dphi4 =  zeros(size(t));
numddphi3 = zeros(size(t));
control_numeric_ddphi3 =  zeros(size(t));
numddphi4 = zeros(size(t));
control_numeric_ddphi4 =  zeros(size(t));


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
    [x, fval, exitflag]=fsolve('loop_closure_eqs_8bar',[y8_init,phi4_init phi5_init,phi6_init,phi7_init,x8_init]',optim_options,phi2(k),r1,r2,r3,r4,r5,r6,r7,AE,FC);
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
    % Set phi3 equal to opposite of phi2
    phi3(k) = pi-phi2(k);
    dphi3(k) = -dphi2(k);
    ddphi3(k) = -ddphi2(k);
    % *** velocity analysis ***
    
    A = [0,-AE*sin(phi4(k)),0,0,r7*sin(phi7(k)),0;
         0,AE*cos(phi4(k)),0,0,-r7*cos(phi7(k)),0;
         0,-r4*sin(phi4(k)),r5*sin(phi5(k)),0,0,-1;
         -1,r4*cos(phi5(k)),-r5*cos(phi5(k)),0,0,0;
         0,0,FC*sin(phi5(k)),-r6*sin(phi6(k)),0,-1;
         -1,0,-FC*cos(phi5(k)),r6*cos(phi6(k)),0,0];

    B = [ r2*sin(phi2(k))*dphi2(k)-r3*sin(phi3(k))*dphi3(k);
         -r2*cos(phi2(k))*dphi2(k)+r3*cos(phi3(k))*dphi3(k);
          r2*sin(phi2(k))*dphi2(k);
         -r2*cos(phi2(k))*dphi2(k);
         r3*sin(phi3(k));
         -r3*cos(phi3(k))];
     
    x = A\B;
    
    % save results
    dy8(k) = x(1);
    dphi4(k) = x(2);
    dphi5(k) = x(3);
    dphi6(k) = x(4);
    dphi7(k) = x(5);
    dx8(k) = x(6);
    
    
    % *** acceleration analysis ***
    
    A = [0,-AE*sin(phi4(k)),0,0,r7*sin(phi7(k)),0;
        0,AE*cos(phi4(k)),0,0,-r7*cos(phi7(k)),0;
        0,-r4*sin(phi4(k)),r5*sin(phi5(k)),0,0,-1;
        -1,r4*cos(phi4(k)),-r5*cos(phi5(k)),0,0,0;
        0,0,FC*sin(phi5(k)),-r6*cos(phi6(k)),0,-1;
        -1,0,-FC*cos(phi5(k)),-r6*sin(phi6(k)),0,0];
     
    B = [r2*sin(phi2(k))*ddphi2(k)+r2*cos(phi2(k))*dphi2(k)^2+-r3*cos(phi3(k))*dphi3(k)^2+AE*cos(phi4(k))*dphi4(k)^2-r7*cos(phi7(k))*dphi7(k)^2-r3*sin(phi3(k))*ddphi3(k);
        -r2*cos(phi2(k))*ddphi2(k)+r2*sin(phi2(k))*dphi2(k)^2-r3*sin(phi3(k))*dphi3(k)^2+AE*dphi4(k)^2*sin(phi4(k))-r7*sin(phi7(k))*dphi7(k)^2+r3*cos(phi3(k))*ddphi3(k);
        r2*sin(phi2(k))*ddphi2(k)+r2*cos(phi2(k))*dphi2(k)^2+r4*sin(phi4(k))*dphi4(k)^2-r5*cos(phi5(k))*dphi5(k)^2;
        -r2*cos(phi2(k))*ddphi2(k)+r2*sin(phi2(k))*dphi2(k)^2+r4*dphi4(k)^2*sin(phi4(k))-r5*sin(phi5(k))*dphi5(k)^2;
        r3*cos(phi3(k))*dphi3(k)^2-FC*cos(phi5(k))*dphi5(k)^2+r6*cos(phi6(k))*dphi6(k)^2+r3*sin(phi3(k))*ddphi3(k);
        r3*sin(phi3(k))*dphi3(k)^2-FC*sin(phi5(k))*dphi5(k)^2+r6*sin(phi6(k))*dphi6(k)^2-r3*cos(phi3(k))*ddphi3(k)];
    
    x = A\B;
    % save results
    ddy8(k) = x(1);
    ddphi4(k) = x(2);
    ddphi5(k) = x(3);
    ddphi6(k) = x(4);
    ddphi7(k) = x(5);
    ddx8(k) = x(6);
    
    
    % *** calculate initial values for next iteration step ***
    y8_init = y8(k)+Ts*dy8(k);
    phi4_init = phi4(k)+Ts*dphi4(k);
    phi5_init = phi5(k)+Ts*dphi5(k);
    phi6_init = phi6(k)+Ts*dphi6(k);
    phi7_init = phi7(k)+Ts*dphi7(k);
    x8_init = x8(k)+Ts*dx8(k);
    
    
    % Controleren door kinematische lus
    Control1(k) = (r2*exp(j*phi2(k))+AE*exp(j*phi4(k)))-(r3*exp(j*phi3(k))+r7*exp(j*phi7(k))+r1);
    Control1(k) = abs(Control1(k));
    Control2(k) = (r3*exp(j*phi3(k))+r6*exp(j*phi6(k))-FC*exp(j*phi5(k)))-x8(k);
    Control2(k) = abs(Control2(k));
    
%     % Controleren door numerieke afgeleide
%     if k>2 && k<t_size
%         numdphi3(k) = (phi3(k+1)-phi3(k+1))/(2*Ts);
%     else
%         numdphi3(k) = 0;
%     end
%     control_numeric_dphi3(k) = numdphi3(k)-dphi3(k);
end % loop over positions



% *** create movie ***

% point  O1 = fixed
O1 = 0;
% point O2 = fixed
O2 = r1;
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
x_right = r1+5*AE;
y_top = r1+2*AE;

figure(10)
hold on
plot([x_left, x_right], [y_bottom, y_top]);
axis equal;
movie_axes = axis;   %save current axes into movie_axes

% draw and save movie frame
for m=1:length(index_vec)
    index = index_vec(m);
    A = O1 + r2 * exp(j*phi2(index));
    E = A + AE * exp(j*phi4(index));
    B = A + r4 * exp(j*phi4(index));
    D = O2 + r3*exp(j*phi3(index));
    %F = D + r5*exp(j*phi5(index));
    C = r1+x8(index)+j*y8(index);
    F = C + FC*exp(j*phi5(index));

    %loop1 = [O1 A E D O2]
    loop1 = [O1 A E B F C];
    loop2 = [O2 D E B F D];
    
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
    O2 = r1;
    A = O1 + r2 * exp(j*phi2(index));
    E = A + AE * exp(j*phi4(index));
    B = A + r4 * exp(j*phi4(index));
    D = O2 + r3*exp(j*phi3(index));
    %F = D + r5*exp(j*phi5(index));
    y8(index)
    C = O2 + x8(index) + j*y8(index);
    %C = B - r5 *exp(j*phi5(index));
    F = C + FC*exp(j*phi5(index));
    
    
    figure
    %assembly1=[O1 A E B C];
    assembly1=[O1 A E B F C];
    assembly2= [O2 D E B F D O2];
    plot(real(assembly1),imag(assembly1),'ro-',real(assembly2),imag(assembly2),'ro-')
    xlabel('[m]')
    ylabel('[m]')
    title('assembly')
    axis equal
    
    figure
    subplot(611)
    plot(t,phi2)
    ylabel('\phi_2 [rad]')
    subplot(612)
    plot(t,phi3)
    ylabel('\phi_3 [rad]')
    subplot(613)
    plot(t,x8)
    ylabel('x_8 [m]')
    xlabel('t [s]')
    subplot(614)
    plot(t,phi5)
    ylabel('\phi_5 [rad]')
    subplot(615)
    plot(t,phi6)
    ylabel('\phi_6 [rad]')
    subplot(616)
    plot(t,phi7)
    ylabel('\phi_7 [rad]')
    
    figure
    subplot(311)
    plot(t,dphi2)
    ylabel('d\phi_2 [rad/s]')
    subplot(312)
    plot(t,dphi4)
    ylabel('d\phi_4 [rad/s]')
    subplot(313)
    plot(t,dx8)
    ylabel('dx_8 [m/s]')
    xlabel('t [s]')
    
    figure
    subplot(311)
    plot(t,ddphi2)
    ylabel('dd\phi_2 [rad/s^2]')
    subplot(312)
    plot(t,ddphi4)
    ylabel('dd\phi_3 [rad/s^2]')
    subplot(313)
    plot(t,ddx8)
    ylabel('ddx_8 [m/s^2]')
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
    plot(t,control_numeric_dphi3)
    ylabel('phi3error')
    xlabel('t[s]')
    
end


