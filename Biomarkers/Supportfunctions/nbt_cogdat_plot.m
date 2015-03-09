function [] = nbt_cogdat_plot(PATH,COND1,COND2,ITEM,statistic)
%
%   Title: nbt_cogdat_plot.m
%   Authors: B.A. Diaz & Rick Jansen
%   Date: 2010-05-30
%   Version: 1.3
%   Description: This function will plot a relative frequency histogram of
%   the response over two conditions of a particular RSQ ITEM as well as a
%   histogran of the differences (changes) between both conditions.
%
% Inputs (can also be called without inputs):

% PATH=path to folder with NBT files to be analyzed (string)
% COND1 = first condition (string)
% COND2 = second condition (string)
% ITEM = RSQ question (string)
% STATISTIC = 1 (ttest) or 0 (default, wilcoxon test)

% For example:
% nbt_cogdat_plot('E:\NBT','ECRN1','ECRN2','Ik voelde me gelukkig.',1)
%
%%                Some flow control & assigning

if nargin<1 || isempty(PATH) disp('Select folder with NBT files to analyze');[path]=uigetdir([],'Select folder with NBT files to analyze');
    if path==0
        return
    end
else path = PATH; end;

%%           get conditions from file names

if nargin<3 || isempty(COND2)
%     d=dir(path);
%     d=d(~[d.isdir]);
%     j=1;
%     for i=1:length(d)
%         if ~isempty(strfind(d(i).name,'mat'));
%             if isempty(strfind(d(i).name,'info'));         %Skip Info files
%                 if isempty(strfind(d(i).name,'analysis'))
%                     FileNameIndex = strfind(d(i).name,'.');
%                     conditions{j} = d(i).name((FileNameIndex(3)+1):(FileNameIndex(4)-1));
%                     j=j+1;
%                 end
%             end
%         end
%     end
%     conditions=unique(conditions);
conditions = nbt_getInfo(path);

    if length(conditions)>2
        [selection,ok]=listdlg('liststring',conditions, 'SelectionMode','single', 'ListSize',[250 300],'PromptString','Select condition 1');
        condition1=conditions{selection};
        [selection,ok]=listdlg('liststring',conditions, 'SelectionMode','single', 'ListSize',[250 300],'PromptString','Select condition 2');
        condition2=conditions{selection};
    else
        condition1=conditions{1};
        condition2=conditions{2};
    end
else
    condition1=COND1;
    condition2=COND2;
end

%%         get rsq data
[c1,s1,p1]=nbt_load_analysis(path,condition1,[],@nbt_get_rsq,cell(1,2),[],[]);
questions1=c1{1};
answers1=c1{2};

[c2,s2,p2]=nbt_load_analysis(path,condition2,[],@nbt_get_rsq,cell(1,2),[],[]);
questions2=c2{1};
answers2=c2{2};

%check if c1 & c2 contain the same subjects
if size(c1,2)~=size(c2,2)
    error('ERROR: Not all subjects are present for both conditions, please add files');
    return;
else
    for i=1:length(s1)
        if ~strcmp(s1{i},s2{i})
            %     if size(find((s1 == s2) == 0),1) > 0
            error('ERROR: The same subjects are not being compared. Check that you only have analysis files for the same subjects in the directory');
            return;
        end
    end
end

if nargin < 4
    [selection,ok]=listdlg('liststring',questions1, 'SelectionMode','single', 'ListSize',[350 850],'PromptString','Select RSQ item');
    if isempty(selection)
        return
    end
    question=questions1{selection};
else
    question=ITEM;
    selection=[];
    for i=1:length(questions1)
        if strcmp(question,questions1{i})
            selection=i;
        end
    end
    if isempty(selection)
        disp('RSQ question is not present in analysis files')
        return
    end
end

if (nargin < 5)
    stat={'mean','median'};
    [selectionS,ok]=listdlg('liststring',stat, 'SelectionMode','single','PromptString','Select statistic');
    statistic=stat{selectionS};
end

disp(' ')
disp('Command window code:')
disp(['nbt_cogdat_plot(',char(39),path,char(39),',',char(39),condition1,char(39),',',char(39),condition2,char(39),',',char(39),question,char(39),',',char(39),statistic,char(39),')'])
disp(' ')

if strcmp(statistic,'mean')
    STATISTIC=1;
else
    STATISTIC=0;
end

%%        get figure handle
done=0;
figHandles = findobj('Type','figure');
for i=1:length(figHandles)
    if strcmp(get(figHandles(i),'Tag'),'RSQ_plot')
        figure(figHandles(i))
        done=1;
    end
end
if ~done
    figure()
    set(gcf,'Tag','RSQ_plot');
end
clf

%Plot the data fullscreen
fullscreen = get(0,'ScreenSize');
if(fullscreen(3) ~= 1 | fullscreen(4) ~= 1)
set(gcf,'Position',[0 0 fullscreen(3) fullscreen(4)])
else
    figure('Position', [1 1 1300 900])
end

%% gather data for comparison  & calculate statistics
var1 = answers1(:,selection);
var2 = answers2(:,selection);

%differences var2-var1
diff = var2-var1;

%n subjects
nSubs = length(var1);

%p-value
if(STATISTIC == 1)
    [h,p] = ttest(var1,var2);
else
    [p,h] = signrank(var1,var2);
end

% if( p < 0.001)
%     p = 'p < .001';
% elseif ( p < 0.01)
%     p = 'p < .01';
% elseif ( p < 0.05)
%     p = 'p < .05';
% end

%frequencies plot A
nb_A = [0.0:5.0];
[freq_A, bins_A] = hist([var1,var2],nb_A);

%frequencies plot B
nb_B = [-5:5];
[freq_B, bins_B] = hist(diff,nb_B);

%relative frequencies
freq_A = (freq_A/nSubs)*100;
freq_B = (freq_B/nSubs)*100;

%histogram A
sbp_A = subplot(1,2,1);
hist_A = bar(bins_A,freq_A,1,'hist');
p1=get(sbp_A,'pos');
p1(4) = p1(4) - 0.2;
set(sbp_A,'pos',p1)
legend(hist_A,{condition1,condition2});
set(sbp_A,'Title',text('string','Relative frequencies of responses','FontSize',14));
% set(sbp_A,'YLim',[0 YMAX1],'FontSize',12)
set(sbp_A,'XLim',[-0.5 5.5],'FontSize',12)
xlabel('Response categories','FontSize',14);
ylabel('Frequency [%]','FontSize',14);
colormap summer

%histogram B
sbp_B = subplot(1,2,2);
hist_B = bar(bins_B,freq_B);
p2=get(sbp_B,'pos');
p2(4) = p2(4) - 0.2;
set(sbp_B,'pos',p2)
legend(hist_B,[condition2,' - ',condition1]);
set(sbp_B,'Title',text('string','Relative frequencies of difference in responses','FontSize',14));
% set(sbp_B,'YLim',[0 YMAX2],'FontSize',12)
set(sbp_B,'XLim',[-5.5 5.5],'FontSize',12)
xlabel('{\Delta} Response categories','FontSize',14);
ylabel('Frequency [%]','FontSize',14);
colormap summer

%Main title
st_l1=['Comparison of RSQ question: ',question,', between condition ',condition2, ' and ',condition1,', N = ',num2str(nSubs),', P = ',num2str(p)];
spT = nbt_suptitle(st_l1);
spT_pos = get(spT,'position');
spT_pos(2) = spT_pos(2) - 0.05;
set(spT,'position',spT_pos);
end
