close all 
clear all

t0=0;      % start time
tf=100;    % end time
Y0=[0.5; 0; 0.7; 0]; %initial state (x,y,vx,vy)

tt=linspace(t0,tf,(tf-t0)*100+1);

mu = 1/82.45; % m2/m1
mus = 1-mu;

dY=@(t,y) F(t,y,mu);

%% Solving the system
% standard Runge-Kutta
opt1=odeset(RelTol=1e-7,AbsTol=1e-7,Stats="on");
[t,Y]=ode45(dY,tt,Y0,opt1);
% higher order method, can use lower tolerence
opt2=odeset(RelTol=1e-5,AbsTol=1e-7,Stats="on");
[t2,Y2]=ode89(dY,tt,Y0,opt2);

%% Plotting
figure
hold on
plot(Y(:,1),Y(:,2))
scatter([-mu,mus],[0,0],50,'filled')
hold off
axis equal
box on
ax=gca;
ax.YLim=[-0.8,0.8];
ax.XLim=[-0.8,1.2];

figure
hold on
plot(Y2(:,1),Y2(:,2))
scatter([-mu,mus],[0,0],50,'filled')
hold off
axis equal
box on
ax=gca;
ax.YLim=[-0.8,0.8];
ax.XLim=[-0.8,1.2];

save ("CR3BP_trajectory.mat","t","Y","mu");

%%
% The governing equations
function yp = F(t,y,mu)
mus = 1-mu;
r1 = norm([y(1)+mu, y(2)]);   % Distance to the m1
r2 = norm([y(1)-mus, y(2)]);  % Distance to the m2
r1_3=r1^3;
r2_3=r2^3;

yp = [y(3);
      y(4);
      2*y(4) + y(1) - mus*(y(1)+mu)/r1_3 - mu*(y(1)-mus)/r2_3;
     -2*y(3) + y(2) - mus*y(2)/r1_3 - mu*y(2)/r2_3;];

end