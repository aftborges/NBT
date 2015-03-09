function FASTER(option_wrapper)

% Copyright (C) 2010 Hugh Nolan, Robert Whelan and Richard Reilly, Trinity College Dublin,
% Ireland
% nolanhu@tcd.ie, robert.whelan@tcd.ie
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
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

startDir=option_wrapper.options.file_options.folder_name;
chan_locs=option_wrapper.options.file_options.channel_locations;
is_bdf=option_wrapper.options.file_options.is_bdf==1;

%chan_locs = check_chan_locs(chan_locs,num_data_chans,num_ext_chans);

if ~ischar(chan_locs)
    fprintf('Invalid channel location file.\n');
    return;
end

if (option_wrapper.options.file_options.resume)
    try
        first_file_num=option_wrapper.options.file_options.current_file_num;
    catch
        first_file_num=1;
    end
    if ~isempty(option_wrapper.options.file_options.plist)
        plist=option_wrapper.options.file_options.plist;
        nlist=option_wrapper.options.file_options.nlist;
    else
        if (is_bdf)
            [plist nlist] = extsearchc(startDir,'.bdf',0);
        else
            % Assume .set if not .bdf
            % Later versions may support other file formats
            [plist nlist] = extsearchc(startDir,'.set',0);
        end

        x=true(size(plist));
        for i=1:length(plist)
            if length(plist{i})>12
                if strcmp(plist{i}(end-11:end),'Intermediate')
                    x(i)=0;
                end
            end
        end
        plist={plist{x}};
        nlist={nlist{x}};

        option_wrapper.options.file_options.plist=plist;
        option_wrapper.options.file_options.nlist=nlist;
    end
else
    first_file_num=1;

    if (is_bdf)
        [plist nlist] = extsearchc(startDir,'.bdf',0);
    else
        % Assume .set if not .bdf
        % Later versions may support other file formats
        [plist nlist] = extsearchc(startDir,'.set',0);
    end
    x=true(size(plist));
    for i=1:length(plist)
        if length(plist{i})>12
            if strcmp(plist{i}(end-11:end),'Intermediate')
                x(i)=0;
            end
        end
    end
    plist={plist{x}};
    nlist={nlist{x}};
    x = findrepeats(plist);
    for i = 1:length(x)
        try
            mkdir(plist{x(i)},nlist{x(i)}(1:end-4));
            movefile([plist{x(i)} filesep nlist{x(i)}],[plist{x(i)} filesep nlist{x(i)}(1:end-4)],'f');
            plist{x(i)} = [plist{x(i)} filesep nlist{x(i)}(1:end-4)];
        catch
            error('Error in organising files in %s\n',plist{x(i)});
            return;
        end
    end

    option_wrapper.options.file_options.plist=plist;
    option_wrapper.options.file_options.nlist=nlist;
end

top_log = fopen([startDir filesep option_wrapper.options.file_options.file_prefix 'FASTER.log'],'a');
n_errors=0;
n_procd=0;
n_skipd=0;

for i=first_file_num:length(plist)
    fprintf(top_log,'%s%s%s:\n',plist{i},filesep,nlist{i});
    if (~isempty(option_wrapper.options.file_options.searchstring))
        searchstring2=option_wrapper.options.file_options.searchstring;
    else
        searchstring2=nlist{i};
    end
    if (~isempty(strfind(nlist{i},searchstring2)) || ~isempty(strfind(plist{i},searchstring2)))
        tic
        fprintf('******************\n');
        fprintf('* File %.3d / %.3d *\n',i,length(nlist));
        fprintf('******************\n');

        log_file = fopen([plist{i} filesep nlist{i}(1:end-4) '.log'],'a');
        option_wrapper.options.file_options.current_file = [plist{i} filesep nlist{i}];
        option_wrapper.options.file_options.current_file_num=i;
        if (~isempty(option_wrapper.options.job_filename))
            option_wrapper.options.job_filename
            save(option_wrapper.options.job_filename,'option_wrapper','-mat');
        end

        % Version number check because try-catch syntax is different for
        % different versions of MATLAB

        try
            FASTER_process(option_wrapper,log_file);
            n_procd=n_procd+1;
        catch
            m=lasterror;
            fprintf('Error - %s.\n',m.message);
            fprintf(top_log,'Error - %s.\n',m.message);
            n_errors=n_errors+1;
            try
                fclose(log_file);
            catch
            end
        end
    else
        n_skipd=n_skipd+1;
        fprintf('Skipped file.\n');
        fprintf(top_log,'Skipped due to filename filter.\n');
    end
end

if (option_wrapper.options.averaging_options.make_GA)
    p=1;
    for i=1:length(plist)
        if (~isempty(option_wrapper.options.file_options.searchstring))
            searchstring2=option_wrapper.options.file_options.searchstring;
        else
            searchstring2=nlist{i};
        end
        if (~isempty(strfind(nlist{i},searchstring2)) || ~isempty(strfind(plist{i},searchstring2))) && exist([plist{i} filesep nlist{i}(1:end-4) '.set'],'file')
            EEGt=pop_loadset('filepath',plist{i},'filename',[nlist{i}(1:end-4) '.set']);
            if ~isempty(option_wrapper.options.averaging_options.GA_markers) && ~isempty(EEGt.data)
                GA_markers=option_wrapper.options.averaging_options.GA_markers;
                for v=1:length(GA_markers)
                    EEGt1=h_epoch(EEGt,GA_markers(v),[EEGt.xmin EEGt.xmax]);
                    if ~isempty(EEGt1.data)
                        GAs{v}(:,:,p)=mean(EEGt1.data,3);
                    end
                end
                p=p+1;
            elseif ~isempty(EEGt.data)
                GAs(:,:,p)=mean(EEGt.data,3);p=p+1;
            end
        end
    end
    if ~isempty(option_wrapper.options.averaging_options.GA_markers) && exist('GAs','var') && ~isempty(GAs)
        for v=1:length(GA_markers)
            if ~isempty(GAs{v})
                cl=EEGt.chanlocs; ci=EEGt.chaninfo;
                EEGt=make_EEG(GAs{v},GA_markers(v),EEGt.srate,[EEGt.xmin EEGt.xmax]);
                EEGt.chanlocs=cl; EEGt.chaninfo=ci;
                if (option_wrapper.options.averaging_options.subject_removal_on) && size(EEGt.data,3)>1
                    list_properties=GA_properties(EEGt,option_wrapper.options.channel_options.eeg_chans,option_wrapper.options.ica_options.EOG_channels);
                    lengths=min_z(list_properties,option_wrapper.options.averaging_options.rejection_options);
                    EEGt=pop_rejepoch(EEGt,find(lengths),0);
                end
                EEGt=eeg_checkset(EEGt);
                pop_saveset(EEGt,'filename',sprintf('GA_%d.set',GA_markers(v)),'filepath',startDir);
            end
        end
    elseif exist('GAs','var') && ~isempty(GAs)
        cl=EEGt.chanlocs; ci=EEGt.chaninfo;
        EEGt=make_EEG(GAs,1,EEGt.srate,[EEGt.xmin EEGt.xmax]);
        EEGt.chanlocs=cl; EEGt.chaninfo=ci;
        if (option_wrapper.options.averaging_options.subject_removal_on) && size(EEGt.data,3)>1
            list_properties=GA_properties(EEGt,option_wrapper.options.channel_options.eeg_chans,option_wrapper.options.ica_options.EOG_channels);
            lengths=min_z(list_properties,option_wrapper.options.averaging_options.rejection_options);
            EEGt=pop_rejepoch(EEGt,find(lengths),0);
        end
        EEGt=eeg_checkset(EEGt);
        pop_saveset(EEGt,'filename','GA.set','filepath',startDir);
    end
end

fprintf('*******************\n');
fprintf('* FASTER Finished *\n');
fprintf('*  %.3d processed  *\n',n_procd);
fprintf('*   %.3d errors    *\n',n_errors);
fprintf('*   %.3d skipped   *\n',n_skipd);
fprintf('*******************\n');
fprintf(top_log,'\nFinished. %d processed, %d errors, %d skipped.\n',n_procd,n_errors,n_skipd);
fclose(top_log);

    function out_chan_locs = check_chan_locs(chan_locs,num_chans,num_exts)
        fid = fopen(chan_locs,'r');
        if fid==-1
            out_chan_locs = -1;
            return;
        end
        count=0;
        x='';
        while ~feof(fid)
            x = fgetl(fid);
            if ~isempty(x)
                count = count+1;
            end
        end
        if count < num_chans + num_exts
            out_chan_locs = -1;
            fclose(fid);
            return;
        end
        if count == num_chans + num_exts
            out_chan_locs = chan_locs;
            fclose(fid);
            return;
        end
        frewind(fid);
        fid2 = fopen(['tempchanloc' chan_locs(end-4:end)],'w');
        count2 = 0;
        while ~feof(fid)
            x = fgets(fid);
            if (count2 < num_chans || count - count2 <= num_exts)
                fprintf(fid2,'%s',x);
            end
            count2 = count2 + 1;
        end
        fclose(fid);
        fclose(fid2);
        out_chan_locs = ['tempchanloc' chan_locs(end-4:end)];
    end

    function indices = findrepeats(input)
        indices=zeros(size(input));
        for u=1:length(input)
            indices(u)=sum(strcmp(input{u},input));
        end
        indices=find(indices>1);
    end

end