function[] =nbt_import_rsq(rsqfile,rsqQfile, NBTpath)

%usage:

%  for example:
%


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



%%
questions=importdata([rsqQfile]); % read in questions
[dummy,TXT,RAW]=xlsread([rsqfile]); % read in rsq file

disp('Adding 1');
for fileID = 2:136 %75 = This should be number of files
    answers = zeros(length(questions),1);
    rsqname = [RAW(fileID,1) '.mat'];
    for ii = 1:length(questions)
        itm = cell2mat(RAW(fileID,3+ii));  % this should be changed to the column that the rsq questions start in - 1
        if isnumeric(itm)
            answers(ii) = itm+1;
            if isnan(itm)
                answers(ii) = nan;
            end
        else
            if strcmp(itm, 'M')
                answers(ii) = 4;
            elseif strcmp(itm, 'F')
                answers(ii) = 0;
            elseif strcmp(itm,'NaN')
                answers(ii) = nan;
            else
                disp('CORRUPT RSQ FILE')
            end
        end
    end
    rsq=nbt_rsq(questions,abs(answers)); % make rsq object
    
    rsqname{1,1}
    % find corresponding analysis file in NBT dir and append rsq object
    nbt_file=[NBTpath,'/',rsqname{1,1},'.mat'];
    if(exist(nbt_file,'file') == 0)
        disp('No matching NBT file found')
    else
        an_file=[NBTpath,'/',rsqname{1,1},'_analysis.mat'];
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
    





