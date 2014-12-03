function QBiomarker = nbt_importXLStoQBiomarker(fileName, SubjectColumn, SignalInfo)

%First we load the xls file
[dummy, dummy, raw] = xlsread(fileName);
QBiomarker = nbt_ARSQ(size(raw,2)-1);

%% Match subject with xls
% Generate subject list
SubjectList = cell(size(raw,1)-1,1);
for i=2:size(raw,1)
    SubjectList{i,1} = raw{i,SubjectColumn};
end
if(ischar(SignalInfo.subjectID))
    SubjectIndex = find(strcmp(SubjectList, SignalInfo.subjectID));
else
    SubjectIndex = find(SubjectList == SignalInfo.subjectID);
end
QIndex = nbt_negSearchVector(1:size(raw,2),SubjectColumn);
%% Insert 'questions' and Answers
for m=1:length(QIndex)
    QBiomarker.Questions{m,1} = raw{1,QIndex(m)};
    QBiomarker.Answers(m) = raw{SubjectIndex,QIndex(m)};
end

end