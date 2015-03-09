path_faster = 'C:\Users\Giusi\Desktop\gNBT\data\human_NeuroPhysiology\HN2012\Hum2011\FasterClean';
path_student = 'C:\Users\Giusi\Desktop\gNBT\data\human_NeuroPhysiology\HN2012\Hum2011\StudentClean'; 
d1 = dir(path_faster);
d2 = dir(path_student);
k = 1;
for i = 1:length(d1)
    if ~isempty(findstr(d1(i).name,'_info.mat')) 
        newd1(k) = d1(i);
        k = k+1;
    end
end

StudentSignalName = 'ICASignal';
k = 1;
for i = 2:length(newd1)
   FileName1 = newd1(i).name;
   [Biom_error,chans_faster,additional_student_chans,faster_student_chans,eye_student_chans,duration_faster_percent,duration_student_percent] = nbt_verify_cleaning(path_faster,path_student,FileName1,StudentSignalName);
   CompareFasterStudent(k).Biom_error = Biom_error;
   CompareFasterStudent(k).chans_faster = chans_faster;
   CompareFasterStudent(k).additional_student_chans = additional_student_chans;
   CompareFasterStudent(k).faster_student_chans = faster_student_chans;
   CompareFasterStudent(k).eye_student_chans = eye_student_chans;
   CompareFasterStudent(k).duration_faster_percent = duration_faster_percent;
   CompareFasterStudent(k).duration_student_percent = duration_student_percent;
   CompareFasterStudent(k).file_name = newd1(i).name;
%    pause(0.1);
%    close all
   k = k+1;
end
   