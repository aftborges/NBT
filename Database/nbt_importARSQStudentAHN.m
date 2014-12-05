function nbt_importARSQStudentAHN(filename, SignalInfo, SaveDir)
idxFilename = strfind(filename,'.');
filename = filename(1:idxFilename(3)-1);
MP14 = MP14importfile([filename '.csv']);
rsq = nbt_ARSQ(125);

ARSQData = importdata([filename '.csv']);
%% Populate the ARSQ 
%First we import the normal ARSQ
    %Questions
    for i=1:55
        IndQ=strfind(MP14{i,1},'"');
        rsq.Questions{i,1} = MP14{i,1}(IndQ(1)+1:IndQ(2)-1);
    end
    %Answers 
    for i=1:55
        rsq.Answers(i,1) = str2double(MP14{57,i});
    end
        

%Add factors
tmpQ = rsq.Answers(1:55);
rsq.Questions{56,1} = 'Discontinuity of Mind';
rsq.Questions{57,1} = 'Theory of Mind';
rsq.Questions{58,1} = 'Self';
rsq.Questions{59,1} = 'Planning';
rsq.Questions{60,1} = 'Sleepiness';
rsq.Questions{61,1} = 'Comfort';
rsq.Questions{62,1} = 'Somatic Awareness';
rsq.Questions{63,1} = 'Health Concern';
rsq.Questions{64,1} = 'Visual Thought';
rsq.Questions{65,1} = 'Verbal Thought';
rsq.Answers(56,1) = nanmean(double(tmpQ([1,11,21])));
rsq.Answers(57,1) = nanmean(double(tmpQ([2,12,22])));
rsq.Answers(58,1) = nanmean(double(tmpQ([3,13,23])));
rsq.Answers(59,1) = nanmean(double(tmpQ([4,14,24])));
rsq.Answers(60,1) = nanmean(double(tmpQ([5,15,25])));
rsq.Answers(61,1) = nanmean(double(tmpQ([6,16,26])));
rsq.Answers(62,1) = nanmean(double(tmpQ([7,17,27])));
rsq.Answers(63,1) = nanmean(double(tmpQ([10,20,30])));
rsq.Answers(64,1) = nanmean(double(tmpQ([8,18,28])));
rsq.Answers(65,1) = nanmean(double(tmpQ([9,19,29])));

%% Then we import additional questionaire data.
%First we identify the sorting variable (music name)
MusicNameTemplate = { 'BH.1';'BH.2';'BR.1';'BR.2'; 'CH.1'; 'CH.2'; 'GR.1' ...
    ; 'GR.2'; 'HD.1';'HD.2'; 'MZ.1';'MZ.2'};
%Then we read the music names
 %here it becomes a bit more complicated because first row is 8 down, next
 %are 4 down.
 MusicName{1,1} = MP14{67,5}(9:12);
 MusicAnswers{1,1} = 67;
for i=2:15
   MusicName{i,1} = MP14{67+5*(i-1),5}(9:12);
   MusicAnswers{i,1} = 67+5*(i-1);
end
%Okay.. then we sort - 
idxMusicName = nbt_searchvector(MusicName,MusicNameTemplate); %this is the sorted index of music names

%additional complication
%  'CH.1' nr. 6, HD.2 nr 12 and MZ.2 nr 15 are repeated. 

%add questions
MusicNameTemplate = { 'BH1';'BH2';'BR1';'BR2'; 'CH1'; 'CH1R'; 'CH2'; 'GR1' ...
    ; 'GR2'; 'HD1';'HD2'; 'HD2R'; 'MZ1';'MZ2';'MZ2R'};
startQ = 66;
for m=1:15
    t = 0;
    for i=68:71
        IndQ=strfind(MP14{i,1},'"');
        rsq.Questions{startQ+t,1} = [MusicNameTemplate{m,1} '.' MP14{i,1}(IndQ(1)+1:IndQ(2)-1)];
        t=t+1;
    end
    startQ = startQ+4;
end
%and then finally we add the answers
startQ = 65;
for m=1:15
    for t=1:4
        rsq.Answers(startQ+t,1) = str2double(MP14{MusicAnswers{idxMusicName(m),1},t}(2));
    end
startQ = startQ+4;
end

rsq = nbt_UpdateBiomarkerInfo(rsq, SignalInfo);
nbt_SaveClearObject('rsq', SignalInfo, SaveDir)

end


function MP14 = MP14importfile(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   MP14 = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   MP14 = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   MP14 = importfile('MP14.S0020.20141203.csv', 1, 151);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2014/12/05 14:18:12

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: text (%s)
%	column4: text (%s)
%   column5: text (%s)
%	column6: text (%s)
%   column7: text (%s)
%	column8: text (%s)
%   column9: text (%s)
%	column10: text (%s)
%   column11: text (%s)
%	column12: text (%s)
%   column13: text (%s)
%	column14: text (%s)
%   column15: text (%s)
%	column16: text (%s)
%   column17: text (%s)
%	column18: text (%s)
%   column19: text (%s)
%	column20: text (%s)
%   column21: text (%s)
%	column22: text (%s)
%   column23: text (%s)
%	column24: text (%s)
%   column25: text (%s)
%	column26: text (%s)
%   column27: text (%s)
%	column28: text (%s)
%   column29: text (%s)
%	column30: text (%s)
%   column31: text (%s)
%	column32: text (%s)
%   column33: text (%s)
%	column34: text (%s)
%   column35: text (%s)
%	column36: text (%s)
%   column37: text (%s)
%	column38: text (%s)
%   column39: text (%s)
%	column40: text (%s)
%   column41: text (%s)
%	column42: text (%s)
%   column43: text (%s)
%	column44: text (%s)
%   column45: text (%s)
%	column46: text (%s)
%   column47: text (%s)
%	column48: text (%s)
%   column49: text (%s)
%	column50: text (%s)
%   column51: text (%s)
%	column52: text (%s)
%   column53: text (%s)
%	column54: text (%s)
%   column55: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
MP14 = [dataArray{1:end-1}];
end
