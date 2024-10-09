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

r1=[CI(1);CI(2)]; % x,y of body 1
r2=[CI(3);CI(4)]; % x,y of body 2
r3=[CI(5);CI(6)]; % x,y of body 3

to = 0;  % start time
tf = 20; % end time
tspan=to:0.001:tf; % time points

%for convenience, masses for body 1,2,3 and Gravitaional constant
m1=1;m2=1;m3=1;G=1;  

dY=@(t,Y) F(t,Y,m1,m2,m3,G);

% solver settings
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,Y]=ode45(dY,tspan,CI,options); %solving system

%% plotting
figure
hold on
plot(Y(:,1),Y(:,2),'r',LineWidth=0.5);
plot(Y(:,3),Y(:,4),'g',LineWidth=0.5);
plot(Y(:,5),Y(:,6),'b',LineWidth=0.5);
hold off
box on

%% creating GIF
[xmin,xmax] = bounds([Y(:,1);Y(:,3);Y(:,5)]); % bounds of x coordinates
[ymin,ymax] = bounds([Y(:,2);Y(:,4);Y(:,6)]); % bounds of y coordinates

gifFile="eight.gif";

figure('Position',[500 500 1920/3 1080/2])
hold on
plot(Y(:,1),Y(:,2),'k-');
hr=scatter(Y(1,1),Y(1,2),30,'r','filled');
hb=scatter(Y(1,3),Y(1,4),30,'b','filled');
hg=scatter(Y(1,5),Y(1,6),30,'g','filled');
hold off
box on
ax=gca;
ax.YLim=[-0.6,0.6];
ax.XLim=[-1.2,1.2];
exportgraphics(gca,gifFile)
for i = 1:100:length(Y(:,1))
    hr.XData=Y(i,1);
    hr.YData=Y(i,2);
    hb.XData=Y(i,3);
    hb.YData=Y(i,4);
    hg.XData=Y(i,5);
    hg.YData=Y(i,6);
    title(sprintf('t=%6.4f',t(i)))        
    exportgraphics(gca,gifFile,"Append",true)
    shg    
end

%% Governing Eqns
function Y=F(t,Y,m1,m2,m3,G)
%  Implementation of the equations of the acceleration of the 3 bodies
%  for example, the 3 eqns in https://en.wikipedia.org/wiki/Three-body_problem

%  Assume 3 bodies moving on xy-plane, z=0 always, hence z omitted

r1=[Y(1);Y(2)]; % x,y of obj 1
r2=[Y(3);Y(4)]; % x,y of obj 2
r3=[Y(5);Y(6)]; % x,y of obj 3

d_12=norm(r1-r2)^3; %distance^3 between 1 & 2
d_13=norm(r1-r3)^3; %distance^3 between 1 & 3
d_23=norm(r2-r3)^3; %distance^3 between 2 & 3

% we have converted the 6 2nd order differential equations to
% 12 1st order equations, so that we can applied standard solver to solve
% the system of differential equations
Y = [Y(7:end); % 6 components 
     G*m2*(r2-r1)./d_12 + G*m3*(r3-r1)./d_13; %2 components
     G*m1*(r1-r2)./d_12 + G*m3*(r3-r2)./d_23; %2 components
     G*m1*(r1-r3)./d_13 + G*m2*(r2-r3)./d_23]; %2 components
end

                                                                                                                                                         