close all 
clear all

t0=0;
tf=6.54;
CI=[0.5; 0; 0.7; 0]; %initial condition [x;y;vx;vy]


tt=linspace(t0,tf,(tf-t0)*100+1);

mu = 1/82.45; % m2/m1
mus = 1-mu;

dY=@(t,Y) F(t,Y,mu);
J=@(t,Y) JF(t,Y,mu);

alpha = 0.999; %fractional order
tic
[t, Y] = flmm2(alpha,dY,J, t0,tf,CI,0.00001,[],3,1e-7,500);
toc
Y=Y.'; %so that the array is in the same form as ode45

%%
figure
hold on
plot(Y(:,1),Y(:,2))
scatter([-mu,mus],[0,0],50,'Filled')
hold off
axis equal
box on
ax=gca;
ax.YLim=[-0.8,0.8];
ax.XLim=[-0.8,1.2];

%save the result for post-processing
save ("FracCR3BP_trajectory.mat","t","Y","mu");

function yp = F(t,y,mu)
mus = 1-mu;
r1 = norm([y(1)+mu, y(2)]);   % Distance to the earth
r2 = norm([y(1)-mus, y(2)]);  % Distance to the moon
r1_3=r1^3;
r2_3=r2^3;

yp = [y(3);
      y(4);
      2*y(4) + y(1) - mus*(y(1)+mu)/r1_3 - mu*(y(1)-mus)/r2_3;
     -2*y(3) + y(2) - mus*y(2)/r1_3 - mu*y(2)/r2_3;];

end

function J=JF(t,Y,mu)
% mu = 1/82.45;
mus = 1-mu;

u1=Y(1);
u2=Y(2);
u3=Y(3);
u4=Y(4);

r1=[ u1+mu; u2];
r2=[u1-mus; u2];
d1 = sqrt(r1(1)^2+r1(2)^2);  % Distance to the earth
d2 = sqrt(r2(1)^2+r2(2)^2);  % Distance to the moon

J=[                                                                    0,                                                           0,  1, 0;
                                                                       0,                                                           0,  0, 1;
3*mu*(mus - u1)^2/d2^5 - mu/d2^3 - mus/d1^3 + 3*mus*(mu + u1)^2/d1^5 + 1,           3*mus*u2*(mu + u1)/d1^5 - 3*mu*u2*(mus - u1)/d2^5,  0, 2;
                       3*mus*u2*(mu + u1)/d1^5 - 3*mu*u2*(mus - u1)/d2^5, 3*mu*u2^2/d2^5 - mu/d2^3 - mus/d1^3 + (3*mus*u2^2)/d1^5 + 1, -2, 0;];
end