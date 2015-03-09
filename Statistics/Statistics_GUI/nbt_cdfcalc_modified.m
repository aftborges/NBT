% nbt_selectrunstatistics - this function is part of the statistics GUI, it allows
% to select specific statistical test and run the test for
% biomarkers/groups/channels or regions
%
% Usage:
%  G = nbt_selectrunstatistics(G);
%
% Inputs:
%   G is the struct variable containing informations on the selected groups
%      i.e.:  G(1).fileslist contains information on the files of Group 1
%           G(1).biomarkerslist list of selected biomarkers for the
%           statistics 
%           G(1).chansregs list of the channels and the regions relected
% Outputs:
%  G updated version of the input, 
%          
%
% Example:
%   G = nbt_selectrunstatistics(G);
%
% References:
% 
% See also: 
%  nbt_load_analysis,
%  nbt_run_stat_group, nbt_run_stat_2groups_or_2conditions,
%  nbt_plot_2conditions_topo, nbt_statisticslog
  
%------------------------------------------------------------------------------------
% Originally created by Giuseppina Schiavone (2012), see NBT website (http://www.nbtwiki.net) for current email address
%------------------------------------------------------------------------------------
%
% ChangeLog - see version control log at NBT website for details.
%
% Copyright (C) <year>  <Main Author>  (Neuronal Oscillations and Cognition group, 
% Department of Integrative Neurophysiology, Center for Neurogenomics and Cognitive Research, 
% Neuroscience Campus Amsterdam, VU University Amsterdam)
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


function G = nbt_selectrunstatistics(G);
global rsqQuests
%----------------------
% Interface
%----------------------
StatSelection = figure('Units','pixels', 'name','NBT: Select Statistics' ,'numbertitle','off','Position',[200  200  410  860.000],...
    'MenuBar','none','NextPlot','new','Resize','off');
% fit figure to screen, adapt to screen resolution
set(StatSelection,'CreateFcn','movegui')
hgsave(StatSelection,'onscreenfig')
close(StatSelection)
StatSelection = hgload('onscreenfig');
currentFolder = pwd;
delete([currentFolder '/onscreenfig.fig']);
% tests
statList{1} ='Grand median plot';
statList{2} ='Normality: Lilliefors test';
statList{3} ='Normality: Shapiro-Wilk test';
statList{4} ='Parametric (Univariate): Student paired t-test';
statList{5} ='Non-Parametric (Univariate): Wilcoxon signed rank test';
statList{6} ='Parametric (Bi-variate): Student unpaired t-test';
statList{7} ='Non-Parametric (Bi-variate):  Wilcoxon rank sum test';
statList{8} ='Non-Parametric (Univariate):  Permutation for difference (means)';
statList{9} ='Non-Parametric (Univariate):  Permutation for difference (medians)';
statList{10} ='Non-Parametric (Bi-variate):  Permutation for correlation';
hp3 = uipanel(StatSelection,'Title','SELECT TEST','FontSize',10,'Units','pixels','Position',[5 600 360 250],'BackgroundColor','w','fontweight','bold');
ListStat = uicontrol(hp3,'Units', 'pixels','style','listbox','Max',1,'Units', 'pixels','Position',[5 5 350 220],'fontsize',10,'String',statList,'BackgroundColor','w');
% biomarkers
hp2 = uipanel(StatSelection,'Title','SELECT BIOMARKER(S)','FontSize',10,'Units','pixels','Position',[5 330 360 250],'BackgroundColor','w','fontweight','bold');
ListBiom = uicontrol(hp2,'Units', 'pixels','style','listbox','Max',length(G(1).biomarkerslist),'Units', 'pixels','Position',[5 5 350 220],'fontsize',10,'String',G(1).biomarkerslist,'BackgroundColor','w');
% regions or channels
reglist{1} = 'Channels';
reglist{2} = 'Regions';
hp = uipanel(StatSelection,'Title','SELECT CHANNELS OR REGIONS','FontSize',10,'Units','pixels','Position',[5 240 360 70],'BackgroundColor','w','fontweight','bold');
ListRegion = uicontrol(hp,'Units', 'pixels','style','listbox','Min',0,'Max',1,'Units', 'pixels','Position',[5 5 350 50],'fontsize',10,'String',reglist,'BackgroundColor','w');

% select Group
NBTelementState = nbt_determineNBTelementState;
for i = 1:length(G)
    if(NBTelementState)
        gname = G(i).selection.group_name;
    else
        gname = G(i).fileslist(1).group_name;
    end
    groupList{i} = ['Group ' num2str(i) ' : ' gname];
end
    
hp4 = uipanel(StatSelection,'Title','SELECT GROUP(S)','FontSize',10,'Units','pixels','Position',[5 110 360 100],'BackgroundColor','w','fontweight','bold');
ListGroup = uicontrol(hp4,'Units', 'pixels','style','listbox','Max',length(groupList),'Units', 'pixels','Position',[5 5 350 80],'fontsize',10,'String',groupList,'BackgroundColor','w');
% run test
RunButton = uicontrol(StatSelection,'Style','pushbutton','String','Run Test','Position',[300 5 100 50],'fontsize',10,'callback',@get_settings);
% join button
joinButton = uicontrol(StatSelection,'Style','pushbutton','String','Join Groups','Position',[5 70 70 30],'fontsize',8,'callback',@join_groups);
% create difference group button
groupdiffButton = uicontrol(StatSelection,'Style','pushbutton','String','Difference Group','Position',[280 70 100 30],'fontsize',8,'callback',@diff_group);
% create difference group button
definegroupButton = uicontrol(StatSelection,'Style','pushbutton','String','Define New Group(s)','Position',[5 40 120 30],'fontsize',8,'callback',@define_new_groups);
% create difference group button
addgroupButton = uicontrol(StatSelection,'Style','pushbutton','String','Update Groups List','Position',[130 40 110 30],'fontsize',8,'callback',@add_new_groups);
% remove group(s)
removeGroupButton = uicontrol(StatSelection,'Style','pushbutton','String','Remove Group(s)','Position',[80 70 100 30],'fontsize',8,'callback',@remove_groups);
% save group(s)
saveGroupButton = uicontrol(StatSelection,'Style','pushbutton','String','Save Group(s)','Position',[185 70 90 30],'fontsize',8,'callback',@save_groups);

% move up
upButton = uicontrol(StatSelection,'Style','pushbutton','String','/\','Position',[370 165 25 25],'fontsize',8,'callback',@up_group);
% move down
downButton = uicontrol(StatSelection,'Style','pushbutton','String','\/','Position',[370 140 25 25],'fontsize',8,'callback',@down_group);


function get_settings(d1,d2)
disp('Waiting for statistics ...')
%----------------------
% get settings
%----------------------
% --- get statistics test (one)
statTest = get(ListStat,'Value');
nameTest = get(ListStat,'String');
nameTest = nameTest(statTest); 
% --- get channels or regions (one)
regs_or_chans_index = get(ListRegion,'Value');
regs_or_chans_name = get(ListRegion,'String');
regs_or_chans_name = regs_or_chans_name(regs_or_chans_index); 
% --- get biomarkers (one or more)
bioms_ind = get(ListBiom,'Value');
bioms_name = get(ListBiom,'String');
bioms_name = bioms_name(bioms_ind);
% --- get group (one or more)
group_ind = get(ListGroup,'Value');
group_name = get(ListGroup,'String');
group_name = group_name(group_ind);
%----------------------
% Load RSQ Questions
%----------------------
% if rsq present load rsqQuestions
rsqpresent = 0;
for i = 1:length(bioms_name)
    if strcmp(bioms_name{i},'rsq.Answers') 
        rsqpresent = 1;
        %load rsq
        rs = load([G(1).fileslist(1).path '/' G(1).fileslist(1).name],'rsq');
        rsq = rs.rsq;
        rsqQuests = rsq.Questions;
        clear rs rsq        
        break;
    end
end




%----------------------
% within a group
%----------------------
if length(group_ind) == 1
    G = evalin('base','G');
    Group = G(group_ind);
    NCHANNELS = length(Group.chansregs.chanloc);
    if(NBTelementState)
%         for l = 1:length(bioms_name)
%             [NBTelementName, Biomarker] = strtok(bioms_name{l},'.');
%             [FreqBand, Biomarker] = strtok(Biomarker, '.');
%             Biomarker = Biomarker(2:end);
%             NBTelementCall = ['nbt_GetData(' NBTelementName ',{Subject,['];
%             SubjectList = nbt_expandCell(Group.selection.sub);
%             for m=1:length(SubjectList)
%                 NBTelementCall = [NBTelementCall ' ' num2str(SubjectList(m)) ];
%             end
%             NBTelementCall = [NBTelementCall '];Condition, ''' Group.selection.con{1,1} ''' ;FrequencyBand,''' FreqBand '''},''' Biomarker ''');' ];
%             B_tmp = evalin('base', NBTelementCall);
%             for m=1:length(SubjectList)
%                  B_values(:,m,l) = B_tmp{m,1};
%             end
%         end
    else
        % load biomarker data from analysis files
        path = Group.fileslist.path;
        n_files = length(Group.fileslist);
%----------------------
% Load biomarkers
%----------------------
       [B_values_cell,Sub,Proj,unit] = nbt_checkif_groupdiff(Group,G,n_files,bioms_name,path);
    end
%----------------------
% check biomarkers dimensionality
%----------------------
    [dimens_diff,biomPerChans,IndexbiomNotPerChans] = nbt_dimension_check(B_values_cell,NCHANNELS);
%----------------------
% biomarkers for channels or regions
%----------------------
    if ~isempty(biomPerChans)
        B_values = nbt_extractBiomPerChans(biomPerChans,B_values_cell);
        % select channels or regions
        if strcmp(regs_or_chans_name,'Channels')
            chans = Group.chansregs.channel_nr;
            B_gebruik(:,:,:) = B_values(chans,:,:);
            regs = [];
        elseif strcmp(regs_or_chans_name,'Regions')
            regions = Group.chansregs.listregdata;
            for j = 1:n_files % subject
                for l = 1:length(bioms_name) % biomarker
                    B = B_values(:,j,l);
                    B_gebruik(:,j,l) = get_regions(B,regions);
                end
            end 
            regs = regions;
        end
        clear B B_values
        % load statistical test
        s = nbt_statisticslog(statTest);
%----------------------
% Run Statistics
%----------------------    
        for i = 1:length(bioms_name(biomPerChans))
            bioms_name2 = bioms_name(biomPerChans);
            disp(['Run ', s.statfuncname, ' for ', bioms_name2{i}])
            B(:,:) = B_gebruik(:,:,i);
%             unit{1,i} = 'tmpset';
            [stat_results(biomPerChans(i))] = nbt_run_stat_group(Group,B,s,bioms_name2{i},regs,unit{1,i});
        end
        assignin('base','stat_results',stat_results)
    end
%----------------------
% biomarkers with different dimensionality
%----------------------
   
    if ~isempty(IndexbiomNotPerChans)
        clear B
%         stat_results = evalin('base','stat_results');
        for dim = 1:length(IndexbiomNotPerChans)
            dimBio = IndexbiomNotPerChans{dim};
            dim_B_values(:,:,:) = nbt_extractBiomNotPerChans(dimBio,B_values_cell);
             % load statistical test
            s = nbt_statisticslog(statTest);
%----------------------
% Run Statistics
%---------------------- 
            for dim2 = 1:length(bioms_name(dimBio))
                bioms_name2 = bioms_name(dimBio);
                disp(['Run ', s.statfuncname, ' for ', bioms_name2{dim2}])
                B(:,:) = dim_B_values(:,:,dim2);
%                 unit{1,dim2} = 'tmpset';
                [stat_results(dimBio(dim2))] = nbt_run_stat_noChansBiom(Group,B,s,bioms_name2{dim2},unit{1,dim2});
                if strcmp(bioms_name2{1},'rsq.Answers')
                    set(get(gca, 'XLabel'), 'String','')
                    set(gca, 'Xtick',1:length(rsqQuests))
                    set(gca, 'Xticklabel',[])
%                     xlim = get(gca,'xlim');
                    ylim = get(gca,'ylim');
                    for ticklab = 1:length(rsqQuests)
                        text(ticklab,ylim(1)-0.1,rsqQuests{ticklab},'rotation',-60,'fontsize',8)
                    end
                end
            end
        end
        assignin('base','stat_results',stat_results)
    end
%----------------------
% Between two groups
%----------------------
elseif length(group_ind) == 2
    G = evalin('base','G');
    Group1 = G(group_ind(1));
    Group2 = G(group_ind(2));
    NCHANNELS = length(Group1.chansregs.chanloc);
    if(NBTelementState)
%        for l = 1:length(bioms_name)
%             [NBTelementName, Biomarker] = strtok(bioms_name{l},'.');
%             [FreqBand, Biomarker] = strtok(Biomarker, '.');
%             Biomarker = Biomarker(2:end);
%             NBTelementCall = ['nbt_GetData(' NBTelementName ',{Subject,['];
%             SubjectList = nbt_expandCell(Group1.selection.sub);
%             for m=1:length(SubjectList)
%                 NBTelementCall = [NBTelementCall ' ' num2str(SubjectList(m)) ];
%             end
%             NBTelementCall = [NBTelementCall '];Condition, ''' Group1.selection.con{1,1} ''' ;FrequencyBand,''' FreqBand '''},''' Biomarker ''');' ];
%             B_tmp = evalin('base', NBTelementCall);
%             for m=1:length(SubjectList)
%                  B_values1(:,m,l) = B_tmp{m,1};
%             end
%        end
%        for l = 1:length(bioms_name)
%             [NBTelementName, Biomarker] = strtok(bioms_name{l},'.');
%             [FreqBand, Biomarker] = strtok(Biomarker, '.');
%             Biomarker = Biomarker(2:end);
%             NBTelementCall = ['nbt_GetData(' NBTelementName ',{Subject,['];
%             SubjectList = nbt_expandCell(Group2.selection.sub);
%             for m=1:length(SubjectList)
%                 NBTelementCall = [NBTelementCall ' ' num2str(SubjectList(m)) ];
%             end
%             NBTelementCall = [NBTelementCall '];Condition, ''' Group2.selection.con{1,1} ''' ;FrequencyBand,''' FreqBand '''},''' Biomarker ''');' ];
%             B_tmp = evalin('base', NBTelementCall);
%             for m=1:length(SubjectList)
%                  B_values2(:,m,l) = B_tmp{m,1};
%             end
%         end
    else
        nameG2 = Group2.fileslist.group_name;
        path2 = Group2.fileslist.path;
        n_files2 = length(Group2.fileslist);
        nameG1 = Group1.fileslist.group_name;
        path1 = Group1.fileslist.path;
        n_files1 = length(Group1.fileslist);
        % load biomarkers
%----------------------
% Load biomarkers
%----------------------
        [B_values1_cell,Sub1,Proj1,unit] = nbt_checkif_groupdiff(Group1,G,n_files1,bioms_name,path1);
        [B_values2_cell,Sub2,Proj2,unit] = nbt_checkif_groupdiff(Group2,G,n_files2,bioms_name,path2);
        
    end
    % check that all biomarkers have same dimensionality
%----------------------
% check biomarkers dimensionality
%----------------------    
    [dimens_diff,biomPerChans,IndexbiomNotPerChans] = nbt_dimension_check(B_values1_cell,NCHANNELS);
%----------------------
% biomarkers for channels or regions
%----------------------
%     stat_results = strcut([]);
%     assignin('base','stat_results',stat_results)
    % compute statistics and plot biomarkers which have values for all the channels
    if ~isempty(biomPerChans)
        B_values1 = nbt_extractBiomPerChans(biomPerChans,B_values1_cell);
        B_values2 = nbt_extractBiomPerChans(biomPerChans,B_values2_cell);
        % select channels or regions
        if strcmp(regs_or_chans_name,'Channels')
            chans = Group1.chansregs.channel_nr;
            B_gebruik1(:,:,:) = B_values1(chans,:,:);
            B_gebruik2(:,:,:) = B_values2(chans,:,:);
            regs = [];
        elseif strcmp(regs_or_chans_name,'Regions')
            regions = Group1.chansregs.listregdata;
            for j = 1:n_files1 % subject
                for l = 1:length(bioms_name) % biomarker
                    B1 = B_values1(:,j,l);
                    B_gebruik1(:,j,l) = get_regions(B1,regions);
                end
            end 
             for j = 1:n_files2 % subject
                for l = 1:length(bioms_name) % biomarker
                    B2 = B_values2(:,j,l);
                    B_gebruik2(:,j,l) = get_regions(B2,regions);
                end
            end 
            regs = regions;
        end
        clear B1 B_values1 B2 B_values2
        % load statistical test
        s = nbt_statisticslog(statTest);
if strcmp(s.statfuncname,'Permutation for correlation') || strcmp(s.statfuncname,'Permutation for difference (means)') || strcmp(s.statfuncname,'Permutation for difference (medians)')
 disp('Test not available for this NBT verison')
else
%----------------------
% Run Statistics
%----------------------
        % compute statistics for biomarker with values on each channels
        bioms_name2 = bioms_name(biomPerChans);
        for i = 1:length(bioms_name2)
            disp(['Run ', s.statfuncname, ' for ', bioms_name2{i}])
            B1(:,:) = B_gebruik1(:,:,i);
            B2(:,:) = B_gebruik2(:,:,i);
            [stat_results(biomPerChans(i))] = nbt_run_stat_2groups_or_2conditions(Group1,Group2,B1,B2,s,bioms_name2{i},regs,unit{1,i});
        end    
        assignin('base','stat_results',stat_results)
%----------------------
% Plot Statistics
%----------------------
        pvaluesmatrix(s,stat_results,regs_or_chans_name,bioms_name(biomPerChans),regs,Group1,Group2,nameG1,nameG2)
end
    end
    % compute statistics and plot biomarkers which have values not for all
    % channels (i.e. biomarker extracted from questionnaires)
%----------------------
% biomarkers with different dimensionality
%----------------------
    if ~isempty(IndexbiomNotPerChans)
        clear B1 B2 B
%         stat_results = evalin('base','stat_results');
    if length(IndexbiomNotPerChans) == 1
        stop = length(IndexbiomNotPerChans);
    else
        stop = length(IndexbiomNotPerChans)-1;
    end
        for dim = 1:stop
            dimBio = IndexbiomNotPerChans{dim};
            dim_B_values1(:,:,:) = nbt_extractBiomNotPerChans(dimBio,B_values1_cell);
            dim_B_values2(:,:,:) = nbt_extractBiomNotPerChans(dimBio,B_values2_cell);
             % load statistical test
            s = nbt_statisticslog(statTest);
if strcmp(s.statfuncname,'Permutation for correlation') || strcmp(s.statfuncname,'Permutation for difference (means)') || strcmp(s.statfuncname,'Permutation for difference (medians)')
 disp('Test not available for this NBT verison')
else
            
%----------------------
% Run Statistics
%----------------------
            for dim2 = 1:length(bioms_name(dimBio))
                bioms_name2 = bioms_name(dimBio);
                disp(['Run ', s.statfuncname, ' for ', bioms_name2{dim2}])
                B1(:,:) = dim_B_values1(:,:,dim2);
                B2(:,:) = dim_B_values2(:,:,dim2);
%                 unit{1,dim2} = 'tmpset';
                
               [stat_results(dimBio(dim2))] = nbt_run_stat2_noChansBiom(Group1,Group2,B1,B2,s,bioms_name2{dim2},unit{1,dim2});
%----------------------
% Plot Statistics
%----------------------
                pvaluesmatrix_noChansBiom(s,stat_results(dimBio),bioms_name(dimBio),Group1,Group2,nameG1,nameG2)
                if strcmp(bioms_name2{1},'rsq.Answers')
                    if isfield(stat_results,'p')
                        set(get(gca, 'YLabel'), 'String','')
                        set(gca, 'Ytick',1:length(rsqQuests))
                        set(gca, 'Yticklabel',[])
    %                     xlim = get(gca,'xlim');
                        xlim = get(gca,'xlim');
                        for ticklab = 1:length(rsqQuests)
                            text(xlim(1)+1,ticklab,rsqQuests{ticklab},'rotation',-60,'fontsize',8)
                        end
                    end
                end
               
            end
            
        end
        assignin('base','stat_results',stat_results)
        
    end
    
    
    end
end
end

%-------------------------------------------------
 function pvaluesmatrix_noChansBiom(s,stat_results,bioms_name,Group1,Group2,nameG1,nameG2)
    if isfield(stat_results,'p')
    % ---------compare p-values test results
        for k = 1:length(stat_results)
            tmp(k,:)  = log10(stat_results(k).p);
            tmp2(k,:) = log10(stat_results(k).p);
        end
        x = tmp';
        y = tmp2';
    
        h1 = figure('Visible','off');
        ah=bar3(y);
        h2 = figure('Visible','on','numbertitle','off','Name',['p-values of biomarkers for ', s.statfuncname],'position',[10          80       1700      500]);
        %--- adapt to screen resolution
        set(h2,'CreateFcn','movegui')
        hgsave(h2,'onscreenfig')
        close(h2)
        h2 = hgload('onscreenfig');
        currentFolder = pwd;
        delete([currentFolder '/onscreenfig.fig']);
        %---
        hh=uicontextmenu;
        hh2 = uicontextmenu;
        bh=bar3(x);
        for i=1:length(bh)
            zdata = get(ah(i),'Zdata');
            set(bh(i),'cdata',zdata);
        end
        axis tight
%     axis image
%     zlabel('mean difference') 
%     axis off
        grid off
        view(-90,-90)
        colorbar('off')
        coolWarm = load('nbt_CoolWarm.mat','coolWarm');
        coolWarm = coolWarm.coolWarm;
        colormap(coolWarm);
        cbh = colorbar('EastOutside');
        minPValue = -2.6;%min(zdata(:));% Plot log10(P-Values) to trick colour bar
        maxPValue = 0;%max(zdata(:));
        caxis([minPValue maxPValue])
        set(cbh,'YTick',[-2 -1.3010 -1 0])
        set(cbh,'YTicklabel',[0.01 0.05 0.1 1]) %(log scale)
        set(cbh,'XTick',[],'XTickLabel','')
%     set(get(cbh,'title'),'String','p-values','fontsize',8,'fontweight','b
%     old');
        DeltaC = (maxPValue-minPValue)/20;
        pos_a = get(gca,'Position');
    %     pos = get(cbh,'Position');
        set(cbh,'Position',[1.5*pos_a(1)+pos_a(3) pos_a(2)+pos_a(4)/3 0.01 0.3 ])
        Pos=get(cbh,'position'); 
        set(cbh,'units','normalized'); 
%         uic=uicontrol('style','slider','units','normalized','position',[Pos(1)-0.015*Pos(1) Pos(2) 0.01 0.3 ],...
%         'min',minPValue,'max',maxPValue-DeltaC,'value',minPValue,...
%         'callback',{@Slider_fun,DeltaC,Pos,maxPValue,bh});
    
           for i = 1:length(bioms_name)
                umenu = text(i,-5, regexprep(bioms_name{i},'_',' '),'horizontalalignment','left','fontsize',8,'fontweight','bold');
                set(umenu,'uicontextmenu',hh);
            end
            set(gca,'YTick',1:5:size(y,1))
            set(gca,'YTickLabel',1:5:size(y,1),'Fontsize',8)
            set(gca,'XTick',[])
            set(gca,'XTickLabel','','Fontsize',8)

        
        set(bh,'uicontextmenu',hh2);
        title(['p-values of biomarkers for ', s.statfuncname, ' for ''', regexprep(nameG2,'_',''),''' vs ''',regexprep(nameG1,'_',''),''''],'fontweight','bold','fontsize',12)

        uimenu(hh,'label','plottest','callback',{@plot_test2Groups2,x,stat_results,Group1,Group2});
        uimenu(hh2,'label','plottest','callback',{@plot_subj_vs_subj2,x,stat_results,Group1,Group2});
        close(h1)  
    else
        disp('Graphic visualization not available.')
    end
    end
% plot pvalue matrix
    function pvaluesmatrix(s,stat_results,regs_or_chans_name,bioms_name,regs,Group1,Group2,nameG1,nameG2)
    if isfield(stat_results,'p')
        
    % ---------compare p-values test results
        for k = 1:length(stat_results)
            tmp(k,:)  = log10(stat_results(k).p);
            tmp2(k,:) = log10(stat_results(k).p);
        end
        x = tmp';
        y = tmp2';
    
        h1 = figure('Visible','off');
        ah=bar3(y);
        h2 = figure('Visible','on','numbertitle','off','Name',['p-values of biomarkers for ', s.statfuncname],'position',[10          80       1700      500]);
        %--- adapt to screen resolution
        set(h2,'CreateFcn','movegui')
        hgsave(h2,'onscreenfig')
        close(h2)
        h2 = hgload('onscreenfig');
        currentFolder = pwd;
        delete([currentFolder '/onscreenfig.fig']);
        %---
        hh=uicontextmenu;
        hh2 = uicontextmenu;
        bh=bar3(x);
        for i=1:length(bh)
            zdata = get(ah(i),'Zdata');
            set(bh(i),'cdata',zdata);
        end
        axis tight
%     axis image
%     zlabel('mean difference') 
%     axis off
        grid off
        view(-90,-90)
        colorbar('off')
        coolWarm = load('nbt_CoolWarm.mat','coolWarm');
        coolWarm = coolWarm.coolWarm;
        colormap(coolWarm);
        cbh = colorbar('EastOutside');
        minPValue = -2.6;%min(zdata(:));% Plot log10(P-Values) to trick colour bar
        maxPValue = 0;%max(zdata(:));
        caxis([minPValue maxPValue])
        set(cbh,'YTick',[-2 -1.3010 -1 0])
        set(cbh,'YTicklabel',[0.01 0.05 0.1 1]) %(log scale)
        set(cbh,'XTick',[],'XTickLabel','')
%     set(get(cbh,'title'),'String','p-values','fontsize',8,'fontweight','b
%     old');
        DeltaC = (maxPValue-minPValue)/20;
        pos_a = get(gca,'Position');
    %     pos = get(cbh,'Position');
        set(cbh,'Position',[1.5*pos_a(1)+pos_a(3) pos_a(2)+pos_a(4)/3 0.01 0.3 ])
        Pos=get(cbh,'position'); 
        set(cbh,'units','normalized'); 
%         uic=uicontrol('style','slider','units','normalized','position',[Pos(1)-0.015*Pos(1) Pos(2) 0.01 0.3 ],...
%         'min',minPValue,'max',maxPValue-DeltaC,'value',minPValue,...
%         'callback',{@Slider_fun,DeltaC,Pos,maxPValue,bh});

        if strcmp(regs_or_chans_name,'Channels')
            for i = 1:length(bioms_name)
                umenu = text(i,-20, regexprep(bioms_name{i},'_',' '),'horizontalalignment','left','fontsize',8,'fontweight','bold');
                set(umenu,'uicontextmenu',hh);
            end
            set(gca,'YTick',1:5:size(y,1))
            set(gca,'YTickLabel',1:5:size(y,1),'Fontsize',8)
            set(gca,'XTick',[])
            set(gca,'XTickLabel','','Fontsize',8)
            ylabel('Channels')
        elseif strcmp(regs_or_chans_name,'Regions')
            for i = 1:length(bioms_name)
                umenu = text(i,-2.5,regexprep(bioms_name{i},'_',' '),'horizontalalignment','left','fontsize',8,'fontweight','bold');
                set(umenu,'uicontextmenu',hh);
            end
            for i= 1:size(x,1)
                text(size(x,2)+0.5,i, regexprep(regs(i).reg.name,'_',' '),'verticalalignment','base','fontsize',8,'rotation',-30,'fontweight','bold');
            end
            set(gca,'XTick',[])
            set(gca,'XTickLabel','','Fontsize',8)
            set(gca,'YTick',[])
            set(gca,'YTickLabel','','Fontsize',8)
            ylabel('Regions')
        end
        set(bh,'uicontextmenu',hh2);
        title(['p-values of biomarkers for ', s.statfuncname, ' for ''', regexprep(nameG2,'_',''),''' vs ''',regexprep(nameG1,'_',''),''''],'fontweight','bold','fontsize',12)

        uimenu(hh,'label','plottest','callback',{@plot_test2Groups,x,stat_results,regs,Group1,Group2});
        uimenu(hh2,'label','plottest','callback',{@plot_subj_vs_subj,x,stat_results,regs,Group1,Group2,regs_or_chans_name});
        close(h1)  
    else
         disp('Graphic visualization not available.')
    
    end
    end
% slider
%     function Slider_fun(d1,d2,DeltaC,Pos,maxPValue,bh)
%         Level=get(gcbo,'Value');
%        col = get(gcf,'color');
%        uicontrol('style','text','units','normalized',...
%            'position',[Pos(1) Pos(2)+Pos(4)+0.1*Pos(4) 0.05 0.02],...
%            'string',['p <= ',sprintf('%4f',10^Level)],'fontsize',8,'BackgroundColor',col,'fontweight','bold');
%          set(findobj(gcf,'tag','Colorbar'),'YLim',[get(gcbo,'Value') maxPValue]);
%          set(findobj(gcf,'tag','Colorbar'),'YTick',linspace(get(gcbo,'Value'),maxPValue,3));
%          set(findobj(gcf,'tag','Colorbar'),'YTicklabel',sprintf('%4f',(10.^(linspace(get(gcbo,'Value'),maxPValue,3)))));
%          for i=1:length(bh)
%             zdata = get(bh(i),'Zdata');
%             ind = find(zdata(:)<=Level);
%             zdata(ind) = Level;
%             set(bh(i),'cdata',zdata);
%          end
%     end

%%
% difference group
function diff_group(d1,d2)
    notavail = 0;
    if notavail==1
        warning('This function is not yet available!')
    else
        G = evalin('base','G');
        % --- get group (one or more)
        group_ind = get(ListGroup,'Value');
        if length(group_ind) >2 
            warning('Select only two for difference group creation!')
        elseif length(group_ind) == 1
            warning('Two groups are necessary for difference group creation!')
        elseif length(group_ind) == 2
            if length(G(group_ind(1)).fileslist) == length(G(group_ind(2)).fileslist)            
                group_name = get(ListGroup,'String');
                group_name = group_name(group_ind);
                G(end+1).fileslist = struct([]);
                for g = 1: length(group_ind)
                    G(end).fileslist = [G(end).fileslist G(group_ind(g)).fileslist];
                    G(end).biomarkerslist = G(group_ind(g)).biomarkerslist;
                    G(end).chansregs = G(group_ind(g)).chansregs;
                end
                scrsz = get(0,'ScreenSize');
                % fit figure to screen, adapt to screen resolution
                hh2 = figure('Units','pixels', 'name','Define group difference' ,'numbertitle','off','Position',[scrsz(3)/4  scrsz(4)/2  250  120],...
                'MenuBar','none','NextPlot','new','Resize','off');
                col =  get(hh2,'Color' );
                set(hh2,'CreateFcn','movegui')
                hgsave(hh2,'onscreenfig')
                close(hh2)
                hh2 = hgload('onscreenfig');
                currentFolder = pwd;
                delete([currentFolder '/onscreenfig.fig']);
                step = 40;

                nameg1 = group_name{1};
                sep = findstr(nameg1,':');
                nameg1 = nameg1(sep+1:end);
                nameg2 = group_name{2};
                sep = findstr(nameg2,':');
                nameg2 = nameg2(sep+1:end);

                text_diff1= uicontrol(hh2,'Style','text','Position',[25 45+step 200 20],'string','Group 1     minus     Group 2','fontsize',10,'fontweight','Bold','BackgroundColor',col);
                text_diff2= uicontrol(hh2,'Style','edit','Position',[25 10+step 80 30],'string',nameg1,'fontsize',10);
                text_diff3= uicontrol(hh2,'Style','text','Position',[115 20+step 20 20],'string',' - ','fontsize',15,'fontweight','Bold','BackgroundColor',col);
                text_diff4= uicontrol(hh2,'Style','edit','Position',[150 10+step 80 30],'string',nameg2,'fontsize',10,'BackgroundColor',col);
                OkButton = uicontrol(hh2,'Style','pushbutton','String','OK','Position',[25 10 200 30],'fontsize',10,'callback',{@confirm_diff_group,G,text_diff2,text_diff4});
            else
                warning('The two groups must have same numebr of subjects!')
            end
        
        end
    end
        
end
    function confirm_diff_group(d1,d2,G,text_diff2,text_diff4)
            nameg1 = get(text_diff2,'string');
%             sep = findstr(nameg1,':');
%             nameg1 = nameg1(sep+1:end);
            nameg2 = get(text_diff4,'string');
%             sep = findstr(nameg2,':');
%             nameg2 = nameg2(sep+1:end);
            
            new_group_name = [nameg1 ' minus ' nameg2] ;
            for g = 1:length(G(end).fileslist)
                G(end).fileslist(g).group_name = new_group_name;
                G(end).group_difference = [nameg1 '-' nameg2];
            end
            groupList{end+1} = ['Group ' num2str(length(G)) ' : ' new_group_name];
            set(ListGroup,'Max',length(groupList),'fontsize',10,'String',groupList,'BackgroundColor','w');

            assignin('base','G',G)
            h = get(0,'CurrentFigure');
            close(h)
        end
%%
% move up
    function up_group(d1,d2)
        G = evalin('base','G');
        group_ind = get(ListGroup,'Value');
        group_name = get(ListGroup,'String');
        group_name = group_name(group_ind(1));
        groupList = get(ListGroup,'String');
        ngroups = length(G);
        groupList2 = groupList;
        G2 = G;
        if group_ind(1)~=1
            G2(group_ind(1)-1) = G(group_ind(1));
            G2(group_ind(1)) = G(group_ind(1)-1);
            groupList2{group_ind(1)-1} = groupList{group_ind(1)};
            groupList2{group_ind(1)} = groupList{group_ind(1)-1};
            G = G2;
            groupList = groupList2;
            assignin('base','G',G)
            set(ListGroup,'Max',length(groupList),'fontsize',10,'String',groupList,'BackgroundColor','w');
        end
        
        
    end
% move down
    function down_group(d1,d2)
        G = evalin('base','G');
        group_ind = get(ListGroup,'Value');
        group_name = get(ListGroup,'String');
        group_name = group_name(group_ind(1));
        groupList = get(ListGroup,'String');
        ngroups = length(G);
        groupList2 = groupList;
        G2 = G;
        if group_ind(1)~=ngroups
            G2(group_ind(1)+1) = G(group_ind(1));
            G2(group_ind(1)) = G(group_ind(1)+1);
            groupList2{group_ind(1)+1} = groupList{group_ind(1)};
            groupList2{group_ind(1)} = groupList{group_ind(1)+1};
            G = G2;
            groupList = groupList2;
            assignin('base','G',G)
            set(ListGroup,'Max',length(groupList),'fontsize',10,'String',groupList,'BackgroundColor','w');
        end
    end
% join groups
    function join_groups(d1,d2)
        G = evalin('base','G');
        % --- get group (one or more)
        group_ind = get(ListGroup,'Value');
        group_name = get(ListGroup,'String');
        group_name = group_name(group_ind);
        G(end+1).fileslist = struct([]);
        for g = 1: length(group_ind)
            G(end).fileslist = [G(end).fileslist G(group_ind(g)).fileslist];
            G(end).biomarkerslist = G(group_ind(g)).biomarkerslist;
            G(end).chansregs = G(group_ind(g)).chansregs;
        end
        new_group_name = (cell2mat(inputdlg('Assign a name to the joined group: ' )));
        for g = 1:length(G(end).fileslist)
            G(end).fileslist(g).group_name = new_group_name;
        end
        groupList{end+1} = ['Group ' num2str(length(G)) ' : ' new_group_name];
        set(ListGroup,'Max',length(groupList),'fontsize',10,'String',groupList,'BackgroundColor','w');
        
        assignin('base','G',G)
        
    end
% add group(s)
function define_new_groups(d1,d2)
    G = evalin('base','G'); 
    G = nbt_definegroups;
end
function add_new_groups(d1,d2)
    G = evalin('base','G');
    newG = length(G);
    groupList = get(ListGroup,'String');
    oldG = length(groupList);
    for i = 1:newG-oldG
        groupList{oldG+i} = ['Group ' num2str(oldG + i) ' : ' G(oldG+i).fileslist(1).group_name];
        G(oldG+i).biomarkerslist = G(oldG).biomarkerslist;
        G(oldG+i).chansregs = G(oldG).chansregs;        
    end
    
    set(ListGroup,'Max',length(groupList),'fontsize',10,'String',groupList,'BackgroundColor','w');
    assignin('base','G',G)
    eval(['evalin(''caller'',''clear SelectedFiles'');']);
end
% remove groups
    function remove_groups(d1,d2)
        G = evalin('base','G');
        group_ind = get(ListGroup,'Value');
        group_name = get(ListGroup,'String');
        if length(group_name)>1
        g3 = 1;
        for g = 1: length(G)
            for g2= 1: length(group_ind)
                if g ~= group_ind(g2)
                    newG(g3) = G(g);
                    groupList2{g3} = ['Group ' num2str(g3) ' : ' G(g).fileslist(1).group_name];
                    g3 = g3+1;
                end
            end
        end
        else
            newG =[];
            groupList2 = {''};
        end
        G = newG;
        groupList = groupList2; 
        set(ListGroup,'String','');
        set(ListGroup,'Value',length(groupList));
        set(ListGroup,'Max',length(groupList),'fontsize',10,'String',groupList,'BackgroundColor','w');
        assignin('base','G',G)
        
    end
% save groups
    function save_groups(d1,d2)
        G = evalin('base','G');
        group_ind = get(ListGroup,'Value');
        group_name = get(ListGroup,'String');
        group_name = group_name(group_ind);
        for g = 1:length(group_ind)
            saveG(g) = G(group_ind(g));
        end
        G = saveG;
        [FileName,PathName,FilterIndex] = uiputfile;
        save([PathName '/' FileName],'G');
        
    end
% get_regions
function B_regions = get_regions(B,regions) 
%         figure
        for i = 1:length(regions);
            chans_in_reg = regions(i).reg.channel_nr;
%             subplot(1,length(regions),i)
%             hist(B(chans_in_reg))
            B_regions(i) = nanmedian(B(chans_in_reg));
%             hold on
%             plot(B_regions(i),20,'r*')
        end
    end
   
%-------------------------------------------------
% plots
function plot_test2Groups(d1,d2,x,stat_results,regs,Group1,Group2)
    pos_cursor_unitfig = get(gca,'currentpoint');
    biomarker = round(pos_cursor_unitfig(1,1));
    if  biomarker>0 && biomarker<= size(x,2)
        s = stat_results(biomarker);
        chanloc = Group1.chansregs.chanloc;
        biom = s.biom_name;
        unit = s.unit;
        nbt_plot_2conditions_topo(Group1,Group2,chanloc,s,unit,biom,regs);
    end
%     
end
%-------------------------------------------------
% plots
function plot_subj_vs_subj(d1,d2,x,stat_results,regs,Group1,Group2,regs_or_chans_name)
    % this plot is valid when using same subject with different
    % conditions--> this implies that groups have same sumber of subject
    if length(Group1.fileslist) == length(Group2.fileslist)
        G1name = Group1.fileslist.group_name;
        G2name = Group2.fileslist.group_name;
        pos_cursor_unitfig = get(gca,'currentpoint');
        chan_or_reg = round(pos_cursor_unitfig(1,2));
        biomarker = round(pos_cursor_unitfig(1,1));
   
    if  biomarker>0 && biomarker<= size(x,2) && chan_or_reg>0 && chan_or_reg<= size(x,1)
        for i = 1:length(Group1.fileslist)
            subname1 = Group1.fileslist(i).name;
            tmp = findstr(subname1,'.');
            sub1{i} = subname1(tmp(1)+1:tmp(2)-1);
            subname2 = Group2.fileslist(i).name;
            tmp = findstr(subname2,'.');
            sub2{i} = subname2(tmp(1)+1:tmp(2)-1);
        end
        s = stat_results(biomarker);
        biom = s.biom_name;
        g(1,:) = s.c1(chan_or_reg,:);
        g(2,:) = s.c2(chan_or_reg,:);
        pval = sprintf('%.4f',s.p(chan_or_reg));
   
        if strcmp(regs_or_chans_name,'Channels')
            h4 = figure('Visible','on','numbertitle','off','Name',[biom ' values for channel ' num2str(chan_or_reg) ' for each subjects'],'Position',[1000   200   350   700]);
        elseif strcmp(regs_or_chans_name,'Regions')
            regname = regexprep(regs(chan_or_reg).reg.name,'_',' ');
            h4 = figure('Visible','on','numbertitle','off','Name',[biom ' values for reagion ' regname ' for each subjects'],'Position',[1000   200   350   700]);
        end
        set(h4,'CreateFcn','movegui')
        hgsave(h4,'onscreenfig')
        close(h4)
        h4= hgload('onscreenfig');
        currentFolder = pwd;
        delete([currentFolder '/onscreenfig.fig']);

        hold on
        plot([1.2 1.8],g,'g')
        for i = 1:length(g)
            text(1.2,g(1,i),sub1{i},'fontsize',5,'horizontalalignment','right')
            text(1.8,g(2,i),sub2{i},'fontsize',5)
        end
        boxplot(g')
        hold on
        plot(1,mean(g(1,:)),'s','Markerfacecolor','k')
        plot(2,mean(g(2,:)),'s','Markerfacecolor','k')
        text(1.02,mean(g(1,:)),'Mean','fontsize',8)
        text(2.02,mean(g(2,:)),'Mean','fontsize',8)
        xlim([0.8 2.2])
        ylim([min(g(:))-0.1*min(g(:)) max(g(:))+0.1*max(g(:))])
         set(gca,'Xtick', [1 2],'XtickLabel',{[G1name,'(n = ',num2str(length(g(1,:))),')'];[G2name,'(n = ',num2str(length(g(1,:))),')']},'fontsize',8,'fontweight','bold')
        xlabel('')
        ylabel([regexprep(biom,'_',' ')],'fontsize',8,'fontweight','bold')
        if strcmp(regs_or_chans_name,'Channels')
             title({[regexprep(biom,'_',' ') ' values for '],...
                 [' channel ''' num2str(chan_or_reg) ''' for each subjects (p = ', pval,')']},'fontweight','bold','fontsize',10)
        elseif strcmp(regs_or_chans_name,'Regions')
            regname = regexprep(regs(chan_or_reg).reg.name,'_',' ');
             title({[regexprep(biom,'_',' ') ' values for '] ,...
                 ['region ''' regname ''' for each subjects (p = ', pval,')']},'fontweight','bold','fontsize',10)
        end
       
       
    end
    else
        G1name = Group1.fileslist.group_name;
        G2name = Group2.fileslist.group_name;
        pos_cursor_unitfig = get(gca,'currentpoint');
        chan_or_reg = round(pos_cursor_unitfig(1,2));
        biomarker = round(pos_cursor_unitfig(1,1));
   
    if  biomarker>0 && biomarker<= size(x,2) && chan_or_reg>0 && chan_or_reg<= size(x,1)
        for i = 1:length(Group1.fileslist)
            subname1 = Group1.fileslist(i).name;
            tmp = findstr(subname1,'.');
            sub1{i} = subname1(tmp(1)+1:tmp(2)-1);
        end
        for i = 1:length(Group2.fileslist)
            subname2 = Group2.fileslist(i).name;
            tmp = findstr(subname2,'.');
            sub2{i} = subname2(tmp(1)+1:tmp(2)-1);
        end
        %----find subject correspondence
        k = 1;
        equalsub{k} = [];
        for i = 1:length(sub1)
            ss = strfind(sub2,sub1{i});
            for j = 1:length(ss)
                if ~isempty(ss{j})
                    equalsub{k}= [i j]; 
                    k = k+1;
                    break;
                end
            end
            clear ss
        end
        
        %----
        s = stat_results(biomarker);
        biom = s.biom_name;
        g1 = s.c1(chan_or_reg,:);
        g2 = s.c2(chan_or_reg,:);
        g = [g1 g2];
        z = [zeros(length(g1),1); ones(length(g2),1)];
        
        pval = sprintf('%.4f',s.p(chan_or_reg));
   
        if strcmp(regs_or_chans_name,'Channels')
            h4 = figure('Visible','on','numbertitle','off','Name',...
                [biom ' values for channel ' num2str(chan_or_reg) ' for each subjects'],...
                'Position',[1000   200   350   700]);
        elseif strcmp(regs_or_chans_name,'Regions')
            regname = regexprep(regs(chan_or_reg).reg.name,'_',' ');
            h4 = figure('Visible','on','numbertitle','off','Name',[biom ' values for reagion ' regname ' for each subjects'],'Position',[1000   200   350   700]);
        end
        set(h4,'CreateFcn','movegui')
        hgsave(h4,'onscreenfig')
        close(h4)
        h4= hgload('onscreenfig');
        currentFolder = pwd;
        delete([currentFolder '/onscreenfig.fig']);

        hold on
        plot(1.2,g1,'g')
        plot(1.8,g2,'g')
        %----plot subject correspondence
        if sum(cellfun('isempty', equalsub))==0
            for nn = 1:length(equalsub)
                egualsubind = equalsub{nn};
                hold on
                plot([1.2 1.8],[g1(egualsubind(1)) g2(egualsubind(2))],'g')
            end
        end
         %----
        for i = 1:length(g1)
            text(1.2,g1(i),sub1{i},'fontsize',5,'horizontalalignment','right')
        end
        for i = 1:length(g2)
            text(1.8,g2(i),sub2{i},'fontsize',5)
        end
        boxplot(g,z);
        hold on
        plot(1,mean(g1),'s','Markerfacecolor','k')
        plot(2,mean(g2),'s','Markerfacecolor','k')
        text(1.02,mean(g1),'Mean','fontsize',8)
        text(2.02,mean(g2),'Mean','fontsize',8)
        xlim([0.8 2.2])
        ylim([min(g(:))-0.1*min(g(:)) max(g(:))+0.1*max(g(:))])
         set(gca,'Xtick', [1 2],'XtickLabel',{[G1name, '(n = ',num2str(length(g1)),')' ];[G2name, '(n = ',num2str(length(g2)),')' ]},'fontsize',8,'fontweight','bold')
        xlabel('')
        ylabel([regexprep(biom,'_',' ')],'fontsize',8,'fontweight','bold')
        if strcmp(regs_or_chans_name,'Channels')
             title({[regexprep(biom,'_',' ') ' values for '],...
                 [' channel ''' num2str(chan_or_reg) ''' for each subjects (p = ', pval,')']},'fontweight','bold','fontsize',10)
        elseif strcmp(regs_or_chans_name,'Regions')
            regname = regexprep(regs(chan_or_reg).reg.name,'_',' ');
             title({[regexprep(biom,'_',' ') ' values for '] ,...
                 ['region ''' regname ''' for each subjects (p = ', pval,')']},'fontweight','bold','fontsize',10)
        end
       
       
    end
    end
%     
end


% plots
function plot_test2Groups2(d1,d2,x,stat_results,Group1,Group2)
    pos_cursor_unitfig = get(gca,'currentpoint');
    biomarker = round(pos_cursor_unitfig(1,1));
    G1name = Group1.fileslist.group_name;
    G2name = Group2.fileslist.group_name;
    if  biomarker>0 && biomarker<= size(x,2)
        s = stat_results(biomarker);
        biom = s.biom_name;
        unit = s.unit;
        h = nbt_plot_2conditions_NoChansBiom(Group1,Group2,s,unit,biom);
    end
    if strcmp(biom,'rsq.Answers')
       set(get(h, 'XLabel'), 'String','')
        set(h, 'Xtick',1:length(rsqQuests))
        set(h, 'Xticklabel',[])
%                     xlim = get(gca,'xlim');
        set(h,'ylim',[-3 3]);
        for ticklab = 1:length(rsqQuests)
            text(ticklab,-3,rsqQuests{ticklab},'rotation',-30,'fontsize',8, 'fontweight','bold')
        end
        if length(Group1.fileslist) == length(Group2.fileslist)
        h5 = figure('Visible','on','numbertitle','off','resize','off','Menubar','none',...
    'Name',['Table with all items and their numbering that have a p-value lower than 0.05' ],...
    'Position',[100   200  1000   500]);
        set(h5,'CreateFcn','movegui')
        hgsave(h5,'onscreenfig')
        close(h5)
        h4= hgload('onscreenfig');
        currentFolder = pwd;
        delete([currentFolder '/onscreenfig.fig']);
        Notsign = find(s.p>=0.05);
        for i = 1:length(s.p)
            strp{i} = (sprintf('%.4f',s.p(i)));
        end
        
        for i = 1:length(Notsign)
        strp{Notsign(i)} = 'N.S.';
        end
        i = 1;
        data = {i,rsqQuests{i}, s.meanc1(i),s.meanc2(i),s.meanc2(i)-s.meanc1(i),strp{i},size(s.c1,2)};
        for i = 2:length(s.p) 
            data{i,1} = i;
            data{i,2} = rsqQuests{i};
            data{i,3} = s.meanc1(i);
            data{i,4} = s.meanc2(i);
            data{i,5} = s.meanc2(i)-s.meanc1(i);
            data{i,6} = strp{i};
            data{i,7} = size(s.c1,2);
        end

        cnames = {'Question','Item',[s.statname '(' regexprep(G1name,'_',' ')  ')'],...
            [s.statname '(' regexprep(G2name,'_',' ') ')'],...
            [s.statname '(' regexprep(G2name,'_',' ')  ') - ' s.statname '(' regexprep(G1name,'_',' ') ')'],...
            [s.statfuncname ' p-value'],'n'};
        columnformat = {'numeric', 'char', 'numeric','numeric','numeric','char','numeric'};
        t = uitable('Parent',h5,'ColumnName',cnames,'ColumnFormat', columnformat,... 
                    'Data',data,'Position',[20 20 950 470],'RowName',[]);
        end
    end     
end
%-------------------------------------------------
% plots
function plot_subj_vs_subj2(d1,d2,x,stat_results,Group1,Group2)
    % this plot is valid when using same subject with different
    % conditions--> this implies that groups have same sumber of subject
    if length(Group1.fileslist) == length(Group2.fileslist)
        G1name = Group1.fileslist.group_name;
        G2name = Group2.fileslist.group_name;
        pos_cursor_unitfig = get(gca,'currentpoint');
        chan_or_reg = round(pos_cursor_unitfig(1,2));
        biomarker = round(pos_cursor_unitfig(1,1));
   
    if  biomarker>0 && biomarker<= size(x,2) && chan_or_reg>0 && chan_or_reg<= size(x,1)
        for i = 1:length(Group1.fileslist)
            subname1 = Group1.fileslist(i).name;
            tmp = findstr(subname1,'.');
            sub1{i} = subname1(tmp(1)+1:tmp(2)-1);
            subname2 = Group2.fileslist(i).name;
            tmp = findstr(subname2,'.');
            sub2{i} = subname2(tmp(1)+1:tmp(2)-1);
        end
        s = stat_results(biomarker);
        biom = s.biom_name;
        g(1,:) = s.c1(chan_or_reg,:);
        g(2,:) = s.c2(chan_or_reg,:);
        pval = sprintf('%.4f', s.p(chan_or_reg));
        h4 = figure('Visible','on','numbertitle','off','Name',[biom ' values for item ' num2str(chan_or_reg) ' for each subjects'],'Position',[1000   200   350   700]);
       
        set(h4,'CreateFcn','movegui')
        hgsave(h4,'onscreenfig')
        close(h4)
        h4= hgload('onscreenfig');
        currentFolder = pwd;
        delete([currentFolder '/onscreenfig.fig']);

        hold on
        plot([1.2 1.8],g,'g')
        for i = 1:length(g)
            text(1.2,g(1,i),sub1{i},'fontsize',5,'horizontalalignment','right')
            text(1.8,g(2,i),sub2{i},'fontsize',5)
        end
        boxplot(g')
        hold on
        plot(1,mean(g(1,:)),'s','Markerfacecolor','k')
        plot(2,mean(g(2,:)),'s','Markerfacecolor','k')
        text(1.02,mean(g(1,:)),'Mean','fontsize',8)
        text(2.02,mean(g(2,:)),'Mean','fontsize',8)
        xlim([0.8 2.2])
        ylim([min(g(:))-0.1*min(g(:)) max(g(:))+0.1*max(g(:))])
        set(gca,'Xtick', [1 2],'XtickLabel',{[regexprep(G1name,'_',' ') , '(n = ',num2str(length(g(1,:))),')' ];[regexprep(G2name,'_',' ') , '(n = ',num2str(length(g(2,:))),')' ]},'fontsize',8,'fontweight','bold')
        xlabel('')
        ylabel([regexprep(biom,'_',' ')],'fontsize',8,'fontweight','bold')
        if ~isempty(rsqQuests) & strcmp(biom,'rsq.Answers')
            title({[regexprep(biom,'_',' ') ' values for '],[' Question ' num2str(chan_or_reg) '.'''  rsqQuests{chan_or_reg} ''''],[' for each subjects (p = ',pval,')']},'fontweight','bold','fontsize',10)
        %----------------------
        answ(:,1) = s.c1(chan_or_reg,:);
        answ(:,2) = s.c2(chan_or_reg,:);
        
        h5 = figure('Visible','on','numbertitle','off','Name',...
            [biom ' for question ' num2str(chan_or_reg) ': Relative Frequencies of Responses'],...
            'Position',[1000   500   450   300]);
        
        set(h5,'CreateFcn','movegui')
        hgsave(h5,'onscreenfig')
        close(h5)
        h4= hgload('onscreenfig');
        currentFolder = pwd;
        delete([currentFolder '/onscreenfig.fig']);
        c_hist = hist(answ,[1 2 3 4 5]);
        c_hist = c_hist/size(answ,1)*100;
        bar(c_hist )
%         set(gca,'xtick',1:6)
        set(gca,'xticklabel',1:5)
        set(gca,'ylim',[0 100])
        xlabel('Scores')
        ylabel('Frequency [%]')
        legend([ regexprep(G1name,'_',' ')  ' (n = ' num2str(size(answ,1)) ')' ],[regexprep(G2name,'_',' ')  ' (n = ' num2str(size(answ,1)) ')' ])
        title({['Relative Frequencies of Responses for Question ', num2str(chan_or_reg) '.' ] ,...
            [''''  rsqQuests{chan_or_reg} '''']},'fontweight', 'bold' )
        h6 = figure('Visible','on','numbertitle','off','Name',...
            [biom ' for question ' num2str(chan_or_reg) ' Difference Distribution'],...
            'Position',[1000   500   450   300]);
       
        set(h6,'CreateFcn','movegui')
        hgsave(h6,'onscreenfig')
        close(h6)
        h6= hgload('onscreenfig');
        currentFolder = pwd;
        delete([currentFolder '/onscreenfig.fig']);
        ansdiff = s.c2(chan_or_reg,:)-s.c1(chan_or_reg,:);
        pdiff = s.p(chan_or_reg);
        d_hist = hist(ansdiff,-4:1:4);
        d_hist = d_hist/size(answ,1)*100;
        bar(d_hist)
        legend([ regexprep(G2name,'_',' ') ' (n = ' num2str(size(answ,1)) ') - ' regexprep(G1name,'_',' ') ' (n = ' num2str(size(answ,1)) ') , p = ' num2str(sprintf('%.4f',pdiff))])
        set(gca,'ylim',[0 100])
        set(gca,'xticklabel',-4:1:4)
        xlabel('Score Difference')
        ylabel('Frequency [%]')
        title({['Difference Distribution for Question ', num2str(chan_or_reg) '.' ] ,...
            [''''  rsqQuests{chan_or_reg} '''']},'fontweight', 'bold' )
        %----------------------
        else
        title({[regexprep(biom,'_',' ') ' values for '],[' item ''' num2str(chan_or_reg) ''' for each subjects (p = ', sprintf('%.4f',pval),')']},'fontweight','bold','fontsize',10)
        end   
        
      
       
    end
   else
        G1name = Group1.fileslist.group_name;
        G2name = Group2.fileslist.group_name;
        pos_cursor_unitfig = get(gca,'currentpoint');
        chan_or_reg = round(pos_cursor_unitfig(1,2));
        biomarker = round(pos_cursor_unitfig(1,1));
   
    if  biomarker>0 && biomarker<= size(x,2) && chan_or_reg>0 && chan_or_reg<= size(x,1)
        for i = 1:length(Group1.fileslist)
            subname1 = Group1.fileslist(i).name;
            tmp = findstr(subname1,'.');
            sub1{i} = subname1(tmp(1)+1:tmp(2)-1);
        end
        for i = 1:length(Group2.fileslist)
            subname2 = Group2.fileslist(i).name;
            tmp = findstr(subname2,'.');
            sub2{i} = subname2(tmp(1)+1:tmp(2)-1);
        end
         %----find subject correspondence
        k = 1;
        equalsub{k} = [];
        for i = 1:length(sub1)
            ss = strfind(sub2,sub1{i});
            for j = 1:length(ss)
                if ~isempty(ss{j})
                    equalsub{k}= [i j]; 
                    k = k+1;
                    break;
                end
            end
            clear ss
        end
        
        %----
        s = stat_results(biomarker);
        biom = s.biom_name;
        g1 = s.c1(chan_or_reg,:);
        g2 = s.c2(chan_or_reg,:);
        g = [g1 g2];
        z = [zeros(length(g1),1); ones(length(g2),1)];
        
        pval = sprintf('%.4f',s.p(chan_or_reg));
   
        h4 = figure('Visible','on','numbertitle','off','Name',[biom ' values for item ' num2str(chan_or_reg) ' for each subjects'],'Position',[1000   200   350   700]);
       
        set(h4,'CreateFcn','movegui')
        hgsave(h4,'onscreenfig')
        close(h4)
        h4= hgload('onscreenfig');
        currentFolder = pwd;
        delete([currentFolder '/onscreenfig.fig']);

        hold on
        plot(1.2 ,g1,'g')
        plot(1.8,g2,'g')
         %----plot subject correspondence
        if sum(cellfun('isempty', equalsub))==0
            for nn = 1:length(equalsub)
                egualsubind = equalsub{nn};
                hold on
                plot([1.2 1.8],[g1(egualsubind(1)) g2(egualsubind(2))],'g')
            end
        end
         %----
        for i = 1:length(g1)
            text(1.2,g1(i),sub1{i},'fontsize',5,'horizontalalignment','right')
        end
        for i = 1:length(g2)
            text(1.8,g2(i),sub2{i},'fontsize',5)
        end
        boxplot(g,z);
        hold on
        plot(1,mean(g1),'s','Markerfacecolor','k')
        plot(2,mean(g2),'s','Markerfacecolor','k')
        text(1.02,mean(g1),'Mean','fontsize',8)
        text(2.02,mean(g2),'Mean','fontsize',8)
        xlim([0.8 2.2])
        ylim([min(g(:))-0.1*min(g(:)) max(g(:))+0.1*max(g(:))])
        set(gca,'Xtick', [1 2],'XtickLabel',{[G1name, '(n = ',num2str(length(g1)),')' ];[G2name, '(n = ',num2str(length(g2)),')' ]},'fontsize',8,'fontweight','bold')
        xlabel('')
        ylabel([regexprep(biom,'_',' ')],'fontsize',8,'fontweight','bold')
        if ~isempty(rsqQuests) & strcmp(biom,'rsq.Answers')
            title({[regexprep(biom,'_',' ') ' values for '],...
                [' Question ' num2str(chan_or_reg) '.'''  rsqQuests{chan_or_reg} ''''],...
                [' for each subjects (p = ', pval,')']},'fontweight','bold','fontsize',10)
        answ1 = s.c1(chan_or_reg,:);
        answ2 = s.c2(chan_or_reg,:);
        
        h5 = figure('Visible','on','numbertitle','off','Name',...
            [biom ' for question ' num2str(chan_or_reg) ': Relative Frequencies of Responses'],...
            'Position',[1000   500   450   300]);
        
        set(h5,'CreateFcn','movegui')
        hgsave(h5,'onscreenfig')
        close(h5)
        h4= hgload('onscreenfig');
        currentFolder = pwd;
        delete([currentFolder '/onscreenfig.fig']);
        c_hist1 = hist(answ1,[1 2 3 4 5]);
        c_hist2 = hist(answ2,[1 2 3 4 5]);
        c_hist1 = c_hist1/size(answ1,2)*100;
        c_hist2 = c_hist2/size(answ2,2)*100;
        tothist = [c_hist1; c_hist2];
        bar(tothist');
       
%         set(gca,'xtick',1:6)
        set(gca,'xtick',1:5)
        set(gca,'xticklabel',1:5)
        set(gca,'ylim',[0 100])
        xlabel('Scores')
        ylabel('Frequency [%]')
        legend([ regexprep(G1name,'_',' ') ' (n = ' num2str(size(answ1,2)) ')' ],[regexprep(G2name,'_',' ') ' (n = ' num2str(size(answ2,2)) ')' ])
        title({['Relative Frequencies of Responses for Question ', num2str(chan_or_reg) '.' ] ,...
            [''''  rsqQuests{chan_or_reg} '''']},'fontweight', 'bold' )
        else
        title({[regexprep(biom,'_',' ') ' values for '],[' item ''' num2str(chan_or_reg) ''' for each subjects (p = ', sprintf('%.4f',pval),')']},'fontweight','bold','fontsize',10)
            
        end
       
    end
    end
%     
end


end