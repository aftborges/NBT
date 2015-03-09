function RSQObject=importRSQdata(RSQObject,fileToRead)

% Import the file xls file
xlsData = importdata(fileToRead);

lastSubject=xlsData.textdata{1,size(xlsData.textdata,2)};
NumSubjects =str2double(lastSubject(2:end));
NumQuestions = size(xlsData.data,1);

RSQmatrix = nan(NumQuestions,NumSubjects);

for i=1:size(xlsData.data,2)
   SubjectID = xlsData.textdata{1,i+1};
   RSQmatrix(:,str2double(SubjectID(2:end))) = xlsData.data(:,i);
end

RSQObject.MarkerValue = RSQmatrix;
RSQObject.LastUpdate = datestr(now);

end



