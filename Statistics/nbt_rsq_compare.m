function [] = rsq_compare(DATA,CONS,LBLS,THRESHOLD)
%RSQ_TTEST Perform paired t-tests on RSQ data
%
%   Title: rsq_ttest.m
%   Author: B.A. Diaz
%   Date: 2010-05-20
%   Version: 1.3
%   Description: This function will print a table of items and correponding
%   p-values for the comparison of the items between two conidtions /
%   trials of choice using the ttest function (paired samples)
%
%   <Input>:
%
%   DATA = 3D-array containing RSQ data (x = subjects, y = items, z =
%   conditions)
%
%   CONS = Vector of the two conditions to compare between. The number of
%   conditions is equal to size(DATA,3). Only TWO conditions can be
%   compaired (for now).
%
%   LBLS = Cell array containing the RSQ items.
%
%   THRESHOLD = Significance threshold (alpha). Use alpha/ nSubs to give
%   Bonferroni correction. If left empty, the table will contain all
%   p-values, regardless of significance of comparison.
%
%   <Ouput>:
%
%   statTable = Table with item names, p-values and assessment of significance.
%

%Some flow control
if (nargin < 4)
    THRESHOLD = 0.05; 
    test = 0;
else
    test = 1;
end;

if (nargin < 3)
    disp('ERROR#1: Need more input arguments');
    return;
end;
if(length(CONS) > size(DATA,3) || length(CONS) < 2)
    disp('ERROR#2: Need at least 2 conditions or more than 2 conditions specified');
    return;
end

%get dims
[n_items, n_subs, n_cons] = size(DATA);
cs = CONS;


%Longest entry in list
%l_entry = size(char(LBLS),2);

%start the comparison loop
header1 = sprintf('=================================================================================================================================\n');
%header2 = sprintf('[Item]\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t   [df] [t]\t\t[p-value]\n');
header2 = sprintf('[#]     [Item]\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t  [MEAN]\t   [MEAN]\t  [DIFF MEAN]  [p-value]\n');
header3 = sprintf('             \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t      [Cond.%d]    [Cond.%d]   [Con.%d-Con.%d]          \n',cs(1),cs(2),cs(2),cs(1));
disp([header1 header2 header3 header1]);

for i = 1:n_items
    [H, P, CI, STATS] = ttest(DATA(i,:,cs(1)),DATA(i,:,cs(2)));
    
    
    item_id = LBLS(i);
    lbl_tmp = char(item_id);
    %lbl_tmp(end+(fix - end)) = ' ';
    
    if(P < THRESHOLD)
        if(P<0.001)
            sigass = 'p < .001';
        elseif(P<0.01)
            sigass = 'p < .01 ';
        elseif(P<0.05)
            sigass = 'p < .05 ';
        end
    elseif(P > THRESHOLD)
        sigass = 'N.S.    ';
    end
    
    %print the stuff
    if(test == 0)
        value = sprintf('(%2d)\t%-75s\t%3.1f\t\t\t%3.1f\t\t\t%4.1f\t\t%s\n',i,lbl_tmp,mean(DATA(i,:,cs(1))),mean(DATA(i,:,cs(2))),mean(DATA(i,:,cs(2)))-mean(DATA(i,:,cs(1))),sigass);
        disp(value);
    elseif(P < THRESHOLD)
        value = sprintf('(%2d)\t%-75s\t%3.1f\t\t\t%3.1f\t\t\t%4.1f\t\t%s\n',i,lbl_tmp,mean(DATA(i,:,cs(1))),mean(DATA(i,:,cs(2))),mean(DATA(i,:,cs(2)))-mean(DATA(i,:,cs(1))),sigass);
        disp(value);
    end
    clearvars P H lbl_tmp
    
    
end