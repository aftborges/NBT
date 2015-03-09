function allConds= nbt_getBiomarkerSeveralConditions(noConditions)

%usage:  nbt_statistics_conditions(directory,condition1,condition2,biomarker,statistic,plotting)

% Inputs:
% directory = path to directoy with NBT files
% condition1 = name of first condition
% condition2 = name of second condition
% biomarker = name of biomarker, as it is written in the analysis file
% statistic = 'mean' (default) or 'median'
% plotting = 0 (default) or 1. If plotting =1, a pdf will be generated with a plot of the statistics

%  for example:
%  nbt_statistics_conditions('C:\NBT','E0R2','ECR2','amplitude_3_7_Hz.MarkerValues','mean',1)

% This function visualizes statistics of the biomarker within and between the
% two conditions

%------------------------------------------------------------------------------------
% Originally created by Rick Jansen (2010), see NBT website for current email address
%------------------------------------------------------------------------------------

% Copyright (C) 2008  Neuronal Oscillations and Cognition group, Department of Integrative Neurophysiology, Center for Neurogenomics and Cognitive Research, Neuroscience Campus Amsterdam, VU University Amsterdam.
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

%% assign fields


disp('select folder with NBT files');[path]=uigetdir();

% get conditions from file names in directory %this does not work if signals are missing!
d=dir(path);
d=d(~[d.isdir]);
j=1;
for i=1:length(d)
    if ~isempty(strfind(d(i).name,'mat'));
        if isempty(strfind(d(i).name,'info'));         %Skip Info files
                
                if ~isempty(strfind(d(i).name,'_'))
                    FileNameIndex = strfind(d(i).name,'_');
                    dotIndex = strfind(d(i).name,'.');
                    conditions{j} = d(i).name((dotIndex(3)+1):(FileNameIndex-1));
                else
                    FileNameIndex = strfind(d(i).name,'.');
                    conditions{j} = d(i).name((FileNameIndex(3)+1):(FileNameIndex(4)-1));
                end
                
                j=j+1;
            
        end
    end
end
conditions=unique(conditions);

% conditions = nbt_getInfo(path);

if length(conditions)>=noConditions
    for i = 1:noConditions
        [selection,ok]=listdlg('liststring',conditions, 'SelectionMode','single', 'ListSize',[250 300],'PromptString',['Select condition ' num2str(i)]);
        eval(['condition' num2str(i) '= char(conditions{selection})']);
    end
else
    error('Not enough conditions found');
end



% get biomarker
d=dir(path);
for i=1:length(d)
    if ~isempty(findstr('analysis',d(i).name))
        [biomarker_objects,biomarkers] = nbt_ExtractBiomarkers([path,'/', d(i).name]);
        break
    end
end
[selection,ok]=listdlg('liststring',biomarker_objects, 'SelectionMode','single', 'ListSize',[250 300],'PromptString','Select biomarker object');
biomarker_object=biomarker_objects{selection};
[selection1,ok]=listdlg('liststring',biomarkers{selection},'SelectionMode','single', 'ListSize',[250 300],'PromptString','Select biomarker');
biomarker=char(strcat(biomarker_object,'.',biomarkers{selection}(selection1)));




SignalInfo=nbt_getInfoObject(path,'Signal'); % get one Info file for locations channels, in case of EEG data
disp(['Computing statistics for ',regexprep(biomarker,'_',' ')])

disp(['for condition ',condition1,' and condition ',condition2])

%% get biomarkers from analysis files
for i = 1:noConditions
    eval(['[c' num2str(i) ',s' num2str(i) ',p' num2str(i) ']=nbt_load_analysis(path,condition' num2str(i) ',biomarker,@nbt_get_biomarker,[],[],[])']);
    eval(['allConds.c' num2str(i) ' = c' num2str(i)]);
end

