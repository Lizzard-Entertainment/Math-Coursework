
clearvars
h =10;      % [m];
v0=5;     % [m*sec^-1];
theta=60;  % [deg];
ro=1.29;   % [kg*m^-3];
g=9.81;    % [kg*m]
d=.053;    % [m];
m=.5;     % [kg];
w=20;      % [m*sec^-1]; angular velocity of a spinning ball
eta=1;     % describes direction(+-) of rotation; 
           % eta=1 is for topspin.
alfa=ro*pi*d^2/(8*m);
time=0:.01:1.5;
ICs=[0, 1, v0*cos(theta*pi/180), v0*sin(theta*pi/180)];

u  = @(t,x)sqrt(x(3).^2+x(4).^2);
CD = @(t,x)(0.508+(1./(22.053+4.196*(u(t,x)./w).^(5/2))).^(2/5));
CM = @(t,x)(1/(2.022+.981*(u(t,x)./w)));
vacuum = @(t,x)([x(3); x(4); 0; -g]);
nospin = @(t,x)([x(3); x(4);
    (-1)*CD(t,x)*alfa*u(t,x)*x(3);
    (-1)*g-CD(t,x)*alfa*u(t,x)*x(4)]);
topspin= @(t,x)([x(3); x(4);
    (-1)*CD(t,x)*alfa*u(t,x)*x(3)+eta*CM(t,x)*alfa*u(t,x)*x(4);
    (-1)*g-CD(t,x)*alfa*u(t,x)*x(4)-eta*CM(t,x)*alfa*u(t,x)*x(3)]);


[tvac, XZvac]= ode23(vacuum, time, ICs, []);
[tns, XZns]  = ode45(nospin, time, ICs, []);
[tts, XZts]  = ode113(topspin, time, ICs, []);

% It is also important to find when the ball hits the ground  
% (x-axis) and after how many seconds.
% For three cases: Vacuum, nospin and top-spin

% Case #1. Vacuum 
z_i=find(abs(XZvac(:,2))<=min(abs(XZvac(:,2))));
t_gr=XZvac(z_i,1);
t_t = time(z_i);

% Case #2. No-spin 
z1_i =find(abs(XZns(:,2))<=min(abs(XZns(:,2))));
t_gr1=XZns(z1_i,1);
t_t1 = time(z1_i);

% Case #3. Top-spin 
z2_i =find(abs(XZts(:,2))<=min(abs(XZts(:,2))));
t_gr2=XZts(z2_i,1);
t_t2 = time(z2_i);
figure 
plot(XZvac(:,1), XZvac(:,2), 'bo', 'linewidth', 1), grid
hold on
plot(XZns(:,1), XZns(:,2), 'm', 'linewidth', 2)
plot(XZts(:,1), XZts(:,2), 'ko-', 'linewidth', 1.0)
legend('In Vacuum','No-Spin in Air','Top-Spin in Air',0)
ylim([0, 3.35])
title 'Trajectory of a Tennis Ball hit under \eta=15^0, h=1 m'
xlabel 'Distance, [m]', ylabel 'Height, [m]'

tt1=['VACUUM: Ball hits the ground (x-axis):' num2str(t_gr) '[m]'];
gtext(tt1);
tt1a=['VACUUM: Ball hits the ground after: ' num2str(t_t) '[sec]'];
gtext(tt1a);
tt2=['Nospin: Ball hits the ground: ' num2str(t_gr1) '[m]'];
gtext(tt2);
tt2a=['Nospin: Ball hits the ground: ' num2str(t_t1)  ' [sec]'];
gtext(tt2a);
tt3=['Topspin: Ball hits the ground: ' num2str(t_gr2)  '[m]'];
gtext(tt3);
tt3a=['Topspin: Ball hits the ground: ' num2str(t_t2) '[sec]'];
gtext(tt3a);
hold off
%
%% Simulink model
eta=1;
open('TENNIS_Ball.mdl')
set_param('TENNIS_Ball','Solver','ode3','StopTime','1.5');
[tout, SIMout]=sim('TENNIS_Ball.mdl');
% plot(tout, SIMout(:,1), 'gd-', tout, SIMout(:,3),'ko', tout, SIMout(:,5),'rh'); 
% title('Simulation of the Simulink model')
% xlabel('time, [sec]'), ylabel('Height, [m]')
% ylim([0, 3.35])
% legend('Vacuum','NO-spin','TOP-spin',0)
figure
plot(SIMout(:,2), SIMout(:,1),'gd-',SIMout(:,4),SIMout(:,3),'ko',SIMout(:,6), SIMout(:,5),'rh')
ylim([0, 3.35])
title('Simulation of the Simulink model')
xlabel('ground, [m]'), ylabel('Height, [m]')
legend('Vacuum','TOP-spin','NO-spin',0)
shg
