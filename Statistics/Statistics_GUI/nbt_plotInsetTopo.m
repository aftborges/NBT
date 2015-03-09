function nbt_plotInsetTopo(chanData,plotMe,colLims)


data = chanData;
news = plotMe;
%data = import channel locations from GSN-HydroCel-129.sfp
%note .SFP not .LOC

[a,b] = sort(data(:,3),1,'descend');
nz = zeros(129,1);
nz(b(1:6)) = 1;
nz(b(7:15)) = 2;
nz(b(16:30)) = 3;
nz(b(31:50)) = 4;
nz(b(51:70)) = 5;
nz(b(71:90)) = 6;
nz(b(91:109)) = 7;
nz(b(110:124)) = 8;
nz(b(125:129)) = 9;
nz(17) = 8;

for i = 1:128
    ox = data(i,1);
    oy = data(i,2);
    oz = data(i,3);
    oz = nz(i);
    ang = atan(ox/oy);
    nx(i) = oz*sin(ang);
    ny(i) = oz*cos(ang);
    
    if oy < 0
        ny(i) = -ny(i);
        nx(i) = -nx(i);
    end

    %text(nx(i),ny(i),num2str(i));
end
ny(129) = 0;
nx(129) = 0;
%axis([-20 20 -20 20])

figure
set(gcf,'position',get(0,'ScreenSize'))
for i = 1:129
    i
    h(i) = axes;
    t = [news(1:i,i)' news(i,i+1:end)] ;
    t(isnan(t)) = 0;
    topoplot(t,'GSN-HydroCel-129.loc');
    set(h(i),'Position',[0.5 + (nx(i)/20) 0.5 + (ny(i)/20) 0.05 0.05]);
    title(i);
    drawnow;

        
end
    coolWarm = load('nbt_CoolWarm.mat','coolWarm');
        coolWarm = coolWarm.coolWarm;
        colormap(coolWarm);
        cbh = colorbar('SouthOutside');
linkprop(get(gcf,'children'),'clim')
set(gca,'clim',colLims);