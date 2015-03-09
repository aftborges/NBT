function[] =nbt_import_rsq(varargin)

%usage:  nbt_import_rsq(RSQdirectory,NBTdirectory2)

%  for example:
%  nbt_import_rsq('C:\RSQ','C:\NBT')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% assigning directories
nargs=length(varargin);
if nargs==0
    rsq_path=(uigetdir([],'Select directory with RSQ files to be imported'));
    if rsq_path==0
        return
    end
    NBT_path=(uigetdir(rsq_path,'Select directory with NBT files'));
    if NBT_path==0
        return
    end
end
if nargs==1
    rsq_path=varargin{1};
    NBT_path=(uigetdir([],'Select directory with NBT files'));
end
if nargs==2
    rsq_path=varargin{1};
    NBT_path=varargin{2};
end

disp(' ')
disp('Command window code:')
disp(['nbt_import_rsq(',char(39),rsq_path,char(39),',',char(39),NBT_path,char(39),')'])
disp(' ')

rsqdir=dir(rsq_path);
rsqdir=rsqdir(~[rsqdir.isdir]);% remove directories

NBTdir=dir(NBT_path);
NBTdir=NBTdir(~[NBTdir.isdir]);% remove directories

%%
questions=importdata([rsq_path,'/','Questions_hn2011.txt']); % read in questions

for i=1:length(rsqdir)
    answers=[];
    rsqname=rsqdir(i).name; % name of rsq file
    if strfind(rsqname,'xls')
        disp(rsqname)
        [dummy,TXT,RAW]=xlsread([rsq_path,'/',rsqname]); % read in rsq file
        answers = zeros(length(questions),1);
        for ii = 1:length(questions)
            itm = cell2mat(RAW(3+ii,2));
            if isnumeric(itm)        
                answers(ii) = itm;
            else
                if strcmp(itm,'nan')
                    answers(ii) = nan;
                else
                    disp('CORRUPT RSQ FILE')
                end
            end
        end
        
        rsq=nbt_rsq(questions,abs(answers)); % make rsq object
        
        % find corresponding analysis file in NBT dir and append rsq object
        nbt_file=[NBT_path,'/',rsqname(1:end-4),'.mat'];
        if(exist(nbt_file,'file') == 0)
            disp('No matching NBT file found')
        else
            an_file=[NBT_path,'/',rsqname(1:end-4),'_analysis.mat'];
            if(exist(an_file,'file') == 2)
                save(an_file,'rsq','-append')
                disp('ATTENTION: Analysis file already exists. Appending to existing file!');
            elseif(exist(an_file,'file') == 0)
                save(an_file,'rsq')
                disp('ATTENTION: There is no Analysis file in the current folder, a new one will be made.');
            end
        end
    end
end






