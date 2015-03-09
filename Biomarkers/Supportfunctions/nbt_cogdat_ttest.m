function [] = nbt_cogdat_ttest(PATH,COND1,COND2,THRESHOLD,statistic)
%
%   Title: nbt_cogdat_ttest.m
%   Authors: B.A. Diaz & Rick Jansen
%   Date: 2010-05-30
%   Version: 1.6
%   Description: This function will print a table of items and corresponding
%   p-values for the comparison of the items between two conditions
%   using t-test or wilcoxon test (paired samples)
%
% Inputs (can also be called without inputs):

% PATH=path to folder with NBT files to be analyzed (string)
% COND1 = first condition (string)
% COND2 = second condition (string)
% THRESHOLD = significance level
% STATISTIC = 'mean' (ttest) or 'median' (default, wilcoxon test)

% For example:
% nbt_cogdat_ttest('E:\NBT','ECRN1','ECRN2',0.05,1)
%
%%                Some flow control & assigning

if nargin<1 || isempty(PATH) disp('Select folder with NBT files to analyze');[path]=uigetdir([],'Select folder with NBT files to analyze');
    if path==0
        return
    end
else path = PATH; end;

%          get conditions from file names

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin < 4)
    THRESHOLD = 0.05;
    disp('ATTENTION: significance level set to 0.05.');
end;

if (nargin < 5)
    stat={'mean','median'};
    [selection,ok]=listdlg('liststring',stat, 'SelectionMode','single','PromptString','Select statistic');
    statistic=stat{selection};
end

disp(' ')
disp('Command window code:')
disp(['nbt_cogdat_ttest(',char(39),path,char(39),',',char(39),condition1,char(39),',',char(39),condition2,char(39),',',num2str(THRESHOLD),',',char(39),statistic,char(39),')'])
disp(' ')

if strcmp(statistic,'mean')
    STATISTIC=1;
else
    STATISTIC=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

n_subs = size(answers1,1);
n_items = length(questions1);

%%
%Table headers
% header1 = sprintf('=================================================================================================================================\n');
if(STATISTIC == 1)
    stat_lbl = 'Student paired samples t-test';
    stat_ev  = 'Mean';
    %     header2 = sprintf('[#]     [Item]\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t  [MEAN]\t   [MEAN]\t  [DIFF MEAN]  [p-value]\n');
    %     header3 = sprintf('\t\t\t\t(%s, N = %d)  \t\t\t\t\t      [Cond.1]    [Cond.2]   [Con.2-Con.1]          \n',stat_lbl,n_subs);
elseif(STATISTIC == 0)
    stat_lbl = 'Wilcoxon sign-rank test';
    stat_ev  = 'Median';
    %     header2 = sprintf('[#]     [Item]\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t  [MEDIAN]\t   [MEDIAN]\t [DIFF MEDIAN]  [p-value]\n');
    %     header3 = sprintf('\t\t\t\t(%s, N = %d)  \t\t\t\t\t\t\t      [Cond.1]     [Cond.2]  [Con.2-Con.1]          \n',stat_lbl,n_subs);
end
% disp([header1 header2 header3 header1]);

scrsz = get(0,'ScreenSize');
if(scrsz(3) ~= 1 | scrsz(4) ~= 1) 
figure('Position',[1 1 scrsz(3) scrsz(4)]);
else
    figure('Position', [1 1 1300 900])
end
cnames = {'Question#','Item',['     ' stat_ev ' (',condition1,')     '],['     ' stat_ev ' (',condition2,')     '],['     (' stat_ev '(',condition2,') - ' stat_ev '(',condition1,')     '],['     ',stat_lbl,' p-value     '],'n'};
try
t = uitable('ColumnName',cnames,'Position',[200 200 scrsz(3)-400 scrsz(4)-400]);
catch
    t = uitable('ColumnName',cnames,'Position',[200 200 1300-400 900-400]);
end 
dat = cell(n_items,7);
%%
%start the comparison loop
for i = 1:n_items
    
    %which statistic to use
    if(STATISTIC == 1)
        nans = isnan(answers1(:,i))|isnan(answers2(:,i));
        if nnz(nans) == n_subs
            H = nan;
            P = nan;
            CI = nan;
            STATS = nan;
            
            mean1 = nan;
            mean2 = nan;
            statVar1 = nan;
            statVar2 = nan;
            numTest = 0;
        else
            [H, P, CI, STATS] = ttest(answers1(find(~nans),i),answers2(find(~nans),i));
            mean1 = mean(answers1(find(~nans),i));
            mean2 = mean(answers2(find(~nans),i));
            
            statVar1 = mean1;
            statVar2 = mean2;
            numTest = nnz(~nans);
        end
        
    elseif(STATISTIC == 0)
        %         error('Non-parametric tests do not work on Student computers.')
        nans = isnan(answers1(:,i))|isnan(answers2(:,i));
        if nnz(nans) == n_subs
            P = nan;
            H = nan;
            STATS = nan;
            median1 = nan;
            median2 = nan;
            statVar1 = nan;
            statVar2 = nan;
            numTest = 0;
        else
            [P, H, STATS] = signrank(answers1(:,i),answers2(:,i));
            median1 = median(answers1(:,i));
            median2 = median(answers2(:,i));
            
            statVar1 = median1;
            statVar2 = median2;
            numTest = nnz(~nans);
        end
    end
    
    item_id = questions1(i);
    lbl_tmp = char(item_id);
    %lbl_tmp(end+(fix - end)) = ' ';
    
    if(~isnan(P))
        sigass = num2str(P);
    else
        sigass = 'failed';
    end
    
    %     %print the stuff
    %     value = sprintf('(%2d)\t%-75s\t%3.1f\t\t\t%3.1f\t\t\t%4.1f\t\t%s\n',i,lbl_tmp,statVar1,statVar2,statVar2-statVar1,sigass);
    %     disp(value);
    
    dat{i,1} = i;
    dat{i,2} = lbl_tmp;
    dat{i,3} = statVar1;
    dat{i,4} = statVar2;
    dat{i,5} = statVar2-statVar1;
    dat{i,6} = sigass;
    dat{i,7} = numTest;
    clearvars P H lbl_tmp
end
set(t,'Data',dat);
set(t,'ColumnWidth',{'auto',500,'auto','auto','auto','auto'});
set(t,'rowname',[]);
set(t,'columnformat',{'char','char','char','char','char','char'});
xlswrite(['RSQStats',condition1,condition2,'.xls'],[cnames;dat]);

end
