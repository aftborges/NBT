function[] = nbt_plot_EEG_channels(varargin)

load('nbt_CoolWarm.mat');
c=varargin{1};
chanlocs = varargin{4};
if ~isempty(varargin{2})
    cmin=varargin{2};
    cmax=varargin{3};
else
    %     cmin=min(c);
    %     cmax=max(c);
    m=max(abs(min(c)),abs(max(c)));
    cmin=-m;
    cmax=m;
end

% if length(varargin)>3
%     color=color(500:end,:);
%     step=(cmax-cmin)/500;
% else
    step=(cmax-cmin)/256;
% end

for i=1:length(c)
    temp(i)= round((c(i)-cmin)/step)+1;
end
temp(temp>256)=256;
temp(temp<1) = 1;
[intx,inty]=nbt_loadintxinty(chanlocs);

for i=1:length(c)
    try
    plot(inty(i),intx(i),'.','color',coolWarm(temp(i),:),'markersize',15)
    hold on
    catch
    end
end

axis off
ew = max([abs(intx) abs(inty)]);
axis([-ew ew -ew ew]);
caxis([cmin cmax])
% colorbar
%axis off
hold off
%   axis equal

end
