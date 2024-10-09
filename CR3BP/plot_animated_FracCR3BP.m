% need to run FracCR3BP.m first to generate the trajectory data

close all
clear all

load("FracCR3BP_trajectory.mat")
x=Y(:,1);
y=Y(:,2);

[xmin,xmax] = bounds([x;-mu;1-mu])
[ymin,ymax] = bounds(y)

gifFile="FracCR3BP.gif";

plotmax = 5000 ;
figure
box on
ax=gca;
ax.YLim=[-0.8,0.8];
ax.XLim=[-0.8,1.2];
exportgraphics(gca,gifFile)
for i = 1:2000:length(Y(:,1))
    if i > plotmax
        plt_rng=i-plotmax:i;        
    else
        plt_rng=1:i;
    end
    clf
    hold on
    scatter([-mu,1-mu],[0,0],50,'r','filled');
    plot(Y(plt_rng,1),Y(plt_rng,2),'b',Y(i,1),Y(i,2),'-ob','MarkerSize',4,MarkerFaceColor='b')
    hold off    
    ax=gca;
    ax.YLim=[-0.8,0.8];
    ax.XLim=[-0.8,1.2];
    box on
    title(sprintf('t=%6.4f',t(i)))        
    exportgraphics(gca,gifFile,"Append",true)
    shg    
end