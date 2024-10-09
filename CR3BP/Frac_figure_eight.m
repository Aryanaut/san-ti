close all
clear all

v1x=0.3471168881189269;
v1y=0.5327249453880302;

CI = [-1; 0; ... %initial position of body 1
       1; 0; ... %initial position of body 2
       0; 0; ... %initial position of body 3
       v1x;v1y; ... %initial velocity of body 1
       v1x;v1y; ... %initial velocity of body 2
      -v1x*2;-v1y*2; ]; %initial velocity of body 3

m1 = 1 ; % mass obj 1 (assume kg)
m2 = 1 ; % ratio of mass obj 2 to mass 1
m3 = 1 ; % ratio of mass obj 3 to mass 1
G = 1 ; % for convinience

r1=[CI(1);CI(2)]; % x,y of obj 1
r2=[CI(3);CI(4)]; % x,y of obj 2
r3=[CI(5);CI(6)]; % x,y of obj 3
cm=(r1*m1+r2*m2+r3*m3)./3;

d_cm1=norm(cm-r1); % initial distance between 1 & cm
d_cm2=norm(cm-r2); % initial distance between 1 & cm
d_cm3=norm(cm-r3); % initial distance between 3 & cm
d=(d_cm1+d_cm2+d_cm3)/3;
Cnorm=sqrt(G*m1/d^3); % this shld give a unit of 1/time

CI(7:end)=CI(7:end)./Cnorm; %perform appropriate normalization for initial velocity

to = 0;
tf = 30;
tspan=to:0.01:tf;
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[t,Y]=ode45(@(t,Y) de(t,Y,m2,m3,d),tspan,CI,options);

figure
hold on
plot(Y(:,1),Y(:,2),'r',LineWidth=0.5)
plot(Y(:,3),Y(:,4),'g',LineWidth=0.5)
plot(Y(:,5),Y(:,6),'b',LineWidth=0.5)
hold off
box on
ax=gca;

alpha = 0.997;
tic
[t3, Y3] = flmm2(alpha,@(t,Y) de(t,Y,m2,m3,d),@(t,Y) Jac(t,Y,m2,m3,d),to,tf,CI,0.0001,[],3,1e-7,500);
toc
Y3=Y3.';
figure
hold on
plot(Y3(:,1),Y3(:,2),'r',LineWidth=0.5)
plot(Y3(:,3),Y3(:,4),'g',LineWidth=0.5)
plot(Y3(:,5),Y3(:,6),'b',LineWidth=0.5)
hold off
box on

%%
figure
tfig=tiledlayout(6,1);
labels={'x1(t)','y1(t)','x2(t)','y2(t)','x2(t)','y2(t)'};
for i =1:6
    nexttile
    hold on
    plot(t,Y(:,i))
    plot(t3,Y3(:,i),'--')
    hold off
    ylabel(labels{i})
    box on
end

tfig.TileSpacing = 'compact';
tfig.Padding = 'compact';

%%

function Y=de(t,Y,m2,m3,d)
%% Implementation 
%  Implementation of the equations of the acceleration of the 3 bodies
%  for example, the 3 eqns in https://en.wikipedia.org/wiki/Three-body_problem
%  Assume 3 bodies moving on xy-plane, z=0 always, hence z omitted

r1=[Y(1);Y(2)]; % x,y of obj 1
r2=[Y(3);Y(4)]; % x,y of obj 2
r3=[Y(5);Y(6)]; % x,y of obj 3

d_12=(norm(r1-r2)/d)^3; %distance^3 between 1 & 2
d_13=(norm(r1-r3)/d)^3; %distance^3 between 1 & 3
d_23=(norm(r2-r3)/d)^3; %distance^3 between 2 & 3

% we have converted the 6 2nd order differential equations to
% 12 1st order equations, so that we can applied standard solver to solve
% the system of differential equations
Y = [Y(7:end); % 6 components 
     m2*(r2-r1)./d_12 + m3*(r3-r1)./d_13; %2 components
        (r1-r2)./d_12 + m3*(r3-r2)./d_23; %2 components
        (r1-r3)./d_13 + m2*(r2-r3)./d_23]; %2 components
end

function J=Jac(t,Y,m2,m3,d)
u1=Y(1);
u2=Y(2);
u3=Y(3);
u4=Y(4);
u5=Y(5);
u6=Y(6);

r1=[Y(1);Y(2)]; % x,y of obj 1
r2=[Y(3);Y(4)]; % x,y of obj 2
r3=[Y(5);Y(6)]; % x,y of obj 3

d_12=norm(r1-r2); %distance between 1 & 2
d_13=norm(r1-r3); %distance between 1 & 3
d_23=norm(r2-r3); %distance between 2 & 3

J=zeros(12,12);
J(1,7)=1;
J(2,8)=1;
J(3,9)=1;
J(4,10)=1;
J(5,11)=1;
J(6,12)=1;
J(7:12,1:6)=d^3*[3*m2*(u1 - u3)^2/d_12^5 - m3/d_13^3 - m2/d_12^3 + 3*m3*(u1 - u5)^2/d_13^5,               3*m2*(u1 - u3)*(u2 - u4)/d_12^5 + 3*m3*(u1 - u5)*(u2 - u6)/d_13^5,                                  (m2*(- 2*(u1 - u3)^2 + (u2 - u4)^2))/d_12^5,                                        -3*m2*(u1 - u3)*(u2 - u4)/d_12^5,                                 -(m3*( 2*(u1 - u5)^2 - (u2 - u6)^2))/d_13^5,                                          -3*m3*(u1 - u5)*(u2 - u6)/d_13^5;
                         3*m2*(u1 - u3)*(u2 - u4)/d_12^5 + 3*m3*(u1 - u5)*(u2 - u6)/d_13^5,       3*m2*(u2 - u4)^2/d_12^5 - m3/d_13^3 - m2/d_12^3 + 3*m3*(u2 - u6)^2/d_13^5,                                             -3*m2*(u1 - u3)*(u2 - u4)/d_12^5,                              -m2*( 2*(u2 - u4)^2 - (u1 - u3)^2 )/d_12^5,                                            -3*m3*(u1 - u5)*(u2 - u6)/d_13^5,                                -m3*( 2*(u2 - u6)^2 - (u1 - u5)^2 )/d_13^5;
                                                   ((-2*(u1 - u3)^2 + (u2 - u4)^2))/d_12^5,                                                   -3*(u1 - u3)*(u2 - u4)/d_12^5,     3*(u1 - u3)^2/d_12^5 - m3/d_23^3 - 1/d_12^3 + (3*m3*(u3 - u5)^2)/d_23^5,          3*(u1 - u3)*(u2 - u4)/d_12^5 + 3*m3*(u4 - u6)*(u3 - u5)/d_23^5,                                   -m3*( 2*(u3 - u5)^2 - (u4 - u6)^2)/d_23^5,                                           -3*m3*(u4 - u6)*(u3 - u5)/d_23^5;
                                                             -3*(u1 - u3)*(u2 - u4)/d_12^5,                                         (( -2*(u2 - u4)^2 + (u1- u3)^2))/d_12^5,               3*(u1 - u3)*(u2 - u4)/d_12^5 + 3*m3*(u3 - u5)*(u4 - u6)/d_23^5,   3*(u2 - u4)^2/d_12^5 - m3/d_23^3 - 1/d_12^3 + 3*m3*(u4 - u6)^2/d_23^5,                                          -(3*m3*(u3 - u5)*(u4 - u6))/d_23^5,                                -m3*( 2*(u4 - u6)^2 - (u3 - u5)^2 )/d_23^5;
                                                   ((-2*(u1 - u5)^2 + (u2 - u6)^2))/d_13^5,                                                   -3*(u1 - u5)*(u2 - u6)/d_13^5,                                    -m2*( 2*(u3 - u5)^2 - (u4 - u6)^2)/d_23^5,                                      -(3*m2*(u4 - u6)*(u3 - u5))/d_23^5,       3*(u1 - u5)^2/d_13^5 - m2/d_23^3 - 1/d_13^3 + 3*m2*(u3 - u5)^2/d_23^5,            3*(u1 - u5)*(u2 - u6)/d_13^5 + 3*m2*(u4 - u6)*(u3 - u5)/d_23^5;
                                                             -3*(u1 - u5)*(u2 - u6)/d_13^5,                                         (( -2*(u2 - u6)^2 + (u1- u5)^2))/d_13^5,                                             -3*m2*(u3 - u5)*(u4 - u6)/d_23^5,                              -(m2*( 2*(u4 - u6)^2)- (u3 - u5)^2)/d_23^5,          (3*(u1 - u5)*(u2 - u6))/d_13^5 + (3*m2*(u3 - u5)*(u4 - u6))/d_23^5,    3*(u2 - u6)^2/d_13^5 - m2/d_23^3 - 1/d_13^3 + 3*m2*(u4 - u6)^2/d_23^5];

end                                                                                                                                                         