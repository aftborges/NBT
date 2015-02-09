%  nbt_plot_2conditions_topoALL - plots an overview of all selected biomarkers - topoplot for statistics with two groups or conditions
%
% Usage:
% nbt_plot_2conditions_topoALL(Group1,Group2,chanloc,s,unit,biomarker,regions)
%
% Inputs:
%   Group1,
%   Group2,
%   chanloc,
%   s,
%   unit,
%   biomarker,
%   regions
%
% Outputs:
%
%
% Example:
%
%
% References:
%
% See also:
%  nbt_plot_EEG_channels, nbt_plot_stat

%------------------------------------------------------------------------------------
% Originally created by Simon-Shlomo Poil, see NBT website (http://www.nbtwiki.net) for current email address
%------------------------------------------------------------------------------------
%
% ChangeLog - see version control log at NBT website for details.
% modified from nbt_plot_2conditions_topo.m (2013)
%
% Copyright (C) 2013 Simon-Shlomo Poil (based on functions developed by Rick Jansen)
%
% Part of the Neurophysiological Biomarker Toolbox (NBT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% See Readme.txt for additional copyright information.
% -------------------------------------------------------------------------
% --------------

function nbt_plot_2conditions_topoAll(Group1,Group2,chanloc,s,unit,biomarker,subplotIndex, MaxSubplotIndex)
condition1 = Group1.selection.group_name;
condition2 = Group2.selection.group_name;

biomarker = strrep(biomarker,'_','.');

meanc1=s.meanc1; %mean / median per channel condition 1
meanc2=s.meanc2; %mean / median per channel condition 2

p = s.p;
c1 = s.c1;
c2 = s.c2;
Cdiff = s.C;
if size(c1,2) == size(c2,2)
    diffC2C1 = c2-c1;
else
    diffC2C1 = [];
end

%%%%special for POWER %%%%
% meanc1 = meanc1/5.12;
% meanc2 = meanc2/5.12;
% diffC2C1 = diffC2C1/5.12;
%%%%%%%%%

diffC2C1_2 = meanc2-meanc1;
statistic=s.statistic;
statfunc = s.statfunc;
statfuncname={s.statfuncname};
fontsize = 10;

% ---for color scale in following plots
vmax=max([meanc1,meanc2]);
vmin=min([meanc1,meanc2]);
cmax = max(vmax);
cmin = min(vmin);
%cmin = 0; %change here for lower limit
NumberOfContours = 6;
levs = linspace(cmin, cmax, NumberOfContours + 2);
MinLevelIndex1 = find(levs > min(meanc1),1,'first');
MinLevelIndex2 = find(levs > min(meanc2),1,'first');
NumberOfContours1 = 6-(MinLevelIndex1-1);
NumberOfContours2 = 6-(MinLevelIndex2-1);
%NumberOfContours1 = 5;
%NumberOfContours2 = 5;
%meanc1(meanc1 < levs(MinLevelIndex1)) = levs(MinLevelIndex1);
%meanc2(meanc2 < levs(MinLevelIndex2)) = levs(MinLevelIndex2);

xa=-2.5;
ya = 0.25;
maxline = 20;

% %---plot grand average condition 1 per channel(Interpolated Plot)
  subplot(4,MaxSubplotIndex,subplotIndex)
  plot_interpolatedTopo(1);
  text(0.1,0.7,biomarker,'horizontalalignment','center')
  axis off
  cbfreeze
  freezeColors
 % drawnow
% 
% % % %---plot grand average condition 2 per channel(Interpolated Plot)
    subplot(4,MaxSubplotIndex,MaxSubplotIndex+subplotIndex)
    plot_interpolatedTopo(2);
    cbfreeze
    freezeColors
%    drawnow

% %---plot grand average difference between conditions or difference between group means (Interpolated Plot)
  subplot(4,MaxSubplotIndex,2*MaxSubplotIndex+subplotIndex)
  plot_interpolatedTopo(3)
  freezeColors
  cbfreeze
 % drawnow

% 
% %---plot P-values for the test (log scaled colorbar)
% minPValue = -2;% Plot log10(P-Values) to trick colour bar -
% maxPValue = -0.5;
% % %red white blue color scale
minPValue = log10(0.0005);
maxPValue = -log10(0.0005);
pLog = log10(p); % to make it log scaled

 if strcmp(char(statfunc),'ttest') || strcmp(char(statfunc),'signrank')
     pLog = sign(statistic(diffC2C1,2))'.*pLog;
 else
     pLog = sign(statistic(diffC2C1_2,2))'.*pLog;
 end
pLog = -1*pLog;
pLog(pLog<minPValue) = minPValue;
pLog(pLog> maxPValue) = maxPValue;


   subplot(4,MaxSubplotIndex,3*MaxSubplotIndex+subplotIndex)
   plot_pTopo()
  cbfreeze
  drawnow


%% Nested functions part
    function plot_interpolatedTopo(ConditionNr)
        if(ConditionNr ==3)
%             CoolWarm = load('nbt_colormapContourBlueRed', 'nbt_colormapContourBlueRed');
%             coolWarm = CoolWarm.nbt_colormapContourBlueRed;
%             colormap(coolWarm);
%             CoolWarm = load('nbt_CoolWarm', 'coolWarm');
%             coolWarm = CoolWarm.coolWarm;
%             colormap(coolWarm);
            Exp_redbluecmap= load('exp_redbluecmap','exp_redbluecmap');
            exp_redbluecmap = Exp_redbluecmap.exp_redbluecmap;
            colormap(exp_redbluecmap);
            if strcmp(char(statfunc),'ttest') || strcmp(char(statfunc),'signrank')
                topoplot(statistic(diffC2C1,2),chanloc,'headrad','rim','numcontour',0,'electrodes','on');
                if(subplotIndex ==1)
                    textThis = sprintf('Grand average for condition %s minus incondition %s ',condition2,condition1);
                end
              
                %mudei isto 
                 
                cmax = max(statistic(diffC2C1,2));
                cmin = min(statistic(diffC2C1,2));
              
                cmax =  max(abs([cmin cmax]));
                cmin = -1.*cmax;
                 
                caxis([cmin cmax])
                
            else
                topoplot(statistic(diffC2C1_2,2),chanloc,'headrad','rim','numcontour',0,'electrodes','on');
                if(subplotIndex == 1)
                    textThis = sprintf('Grand average for group %s minus group %s',condition2,condition1);
                end
                cmax =  max(abs([min(diffC2C1_2) max(diffC2C1_2)]));
                cmin = -1.*cmax;
                
%                 cmax = max(statistic(diffC2C1,2));
%                 cmin = min(statistic(diffC2C1,2));
%                 
                caxis([cmin cmax])
            end
            
        else
%             nbt_redwhite = load('nbt_colormapContourWhiteRed', 'nbt_colormapContourWhiteRed');
%             nbt_redwhite = nbt_redwhite.nbt_colormapContourWhiteRed;
%             colormap(nbt_redwhite);
            Reds5 = load('Reds5','Reds5');
            Reds5 = Reds5.Reds5;
            colormap(Reds5);
            if(ConditionNr == 1)
%                 topoplot(meanc1',chanloc,'headrad','rim','numcontour',NumberOfContours1,'electrodes','on');
                topoplot(meanc1',chanloc,'headrad','rim','numcontour',0,'electrodes','on');
                cmin = min(min(meanc1),min(meanc2)); cmax = max(max(meanc1),max(meanc2));
            else
%                 topoplot(meanc2',chanloc,'headrad','rim','numcontour',NumberOfContours2,'electrodes','on');
                topoplot(meanc2',chanloc,'headrad','rim','numcontour',0,'electrodes','on');
              
                 cmin = min(min(meanc1),min(meanc2)); cmax = max(max(meanc1),max(meanc2));
            end
            caxis([cmin,cmax])
            if(subplotIndex == 1)
                if strcmp(char(statfunc),'ttest') || strcmp(char(statfunc),'signrank')
                    if(ConditionNr == 1)
                        textThis = sprintf('Grand average for  %s (n = %i)',condition1,size(c1,2));
                    else
                        textThis = sprintf('Grand average for condition: %s (n = %i)',condition2, size(c2,2));
                    end
                else
                    if(ConditionNr == 1)
                        textThis = sprintf('Grand average for group: %s (n = %i) ',condition1,size(c1,2));
                    else
                        textThis = sprintf('Grand average for group: %s (n = %i) ',condition2, size(c2,2));
                    end
                end
            end
        end
        cb = colorbar('westoutside');
        set(get(cb,'title'),'String',unit);
        cin = (cmax-cmin)/6;
        %set(cb,'YTick',round([cmin:cin:cmax]/0.01)*0.01);
       
        
        if(subplotIndex == 1)
         nbt_split_text(xa,ya,textThis,maxline,fontsize);
         axis off
         set(gca,'fontsize',fontsize)
        end
    end

    function plot_pTopo(varargin)
        if(~isempty(varargin))
            figure;
        end
        CoolWarm = load('nbt_DarkBlueWhiteDarkRedSharp', 'nbt_DarkBlueWhiteDarkRedSharp');
        CoolWarm = CoolWarm.nbt_DarkBlueWhiteDarkRedSharp;
%        CoolWarm = load('nbt_colortmap', 'nbt_colortmap');
 %       CoolWarm = CoolWarm.nbt_colortmap;
        colormap(CoolWarm);
        topoplot(pLog,chanloc,'headrad','rim','numcontour',0,'electrodes','on')
        cb = colorbar('westoutside');
        caxis([minPValue maxPValue])
        
        axis square
        set(cb,'YTick',[-2.3010 -1.3010 0 1.3010 2.3010])
        set(cb,'YTicklabel',[0.005 0.05 0 0.05 0.005])
    end
end
