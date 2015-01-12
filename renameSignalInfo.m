function renameSignalInfo(SignalInfo)


SignalInfo.file_name=SignalInfo.file_name(1,1:23);

filename=strcat('SignalInfo.file_name','_info.mat')

save(filename,SignalInfo)