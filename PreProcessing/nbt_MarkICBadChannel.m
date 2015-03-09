function EEG=nbt_MarkICBadChannel(EEG)

answer = inputdlg('ICA components to mark?');
answer = str2num(answer{1});

for i=1:length(answer)
    [dummy BadChannel(i)] = max(abs( EEG.icawinv(:,answer(i))));
    BadChannel(i) = EEG.icachansind(BadChannel(i));
end

indelec =zeros(EEG.nbchan,1);
indelec(BadChannel) = 1;

if(~isempty(EEG.NBTinfo.BadChannels))
    EEG.NBTinfo.BadChannels(find(indelec)) = 1;
else
    EEG.NBTinfo.BadChannels = indelec;
end
disp(' ');
disp('You should now rerun ICA, no components have been rejected, but there are is now a new bad channel list');
find(EEG.NBTinfo.BadChannels)
end