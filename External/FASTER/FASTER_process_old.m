function EEG=FASTER_process(option_wrapper,log_file)

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

tic;
o=option_wrapper.options;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File options %
%%%%%%%%%%%%%%%%
% 1 File name including full path (string)
% 2 Reference channel (integer > 0)
% 3 Number of data channels (integer > 0)
% 4 Number of extra channels (integer > 0)
% 5 Channel locations file including full path (string)
% 6 Save options (cell)
%%%%%%%%%%%%%%%%
fullfilename = o.file_options.current_file;
ref_chan = o.channel_options.ref_chan;
eeg_chans = o.channel_options.eeg_chans;
if (eeg_chans==0)
    eeg_chans=[];
end
ext_chans = o.channel_options.ext_chans;
if (ext_chans==0)
    ext_chans=[];
end
channel_locations_file = o.file_options.channel_locations;
save_options = o.save_options;
cutoff_markers = o.file_options.cutoff_markers;
do_reref = o.channel_options.do_reref;

if (~do_reref)
    ref_chan=[];
end
[filepath,filename,extension,version] = fileparts(fullfilename);

mkdir([filepath filesep 'Intermediate']);

%log_file = fopen([filepath filesep filename '.log'],'a');

c=clock;
months{1}='Jan';months{2}='Feb';months{3}='Mar';months{4}='Apr';months{5}='May';months{6}='Jun';
months{7}='Jul';months{8}='Aug';months{9}='Sep';months{10}='Oct';months{11}='Nov';months{12}='Dec';
fprintf(log_file,'\n%d/%s/%d %d:%d:%d\n',c(3),months{c(2)},c(1),c(4),c(5),round(c(6)));
fprintf(log_file,'%.2f - Opened log file.\n',toc);

%%%%%%%%%%%%%%%%%%%%%%
% File setup section %
%%%%%%%%%%%%%%%%%%%%%%

% Import .bdf file or load .set file
% Note: import all channels and then remove the unnecessary ones, as
% otherwise the event channel gets removed and we have no event data.
if strcmpi(extension,'.bdf')
    fprintf('Importing %s.\n',fullfilename);
    EEG = pop_biosig(fullfilename);
    EEG.setname = filename;
    EEG = pop_select(EEG, 'nochannel',length(eeg_chans)+length(ext_chans)+1:size(EEG.data,1));
    if (do_reref)
        EEG = h_pop_reref( EEG, ref_chan, 'exclude',ext_chans, 'keepref', 'on');
    end
    EEG = pop_saveset(EEG,'filename',[filename '.set'],'filepath',filepath,'savemode','onefile');
    fprintf(log_file,'%.2f - Imported and converted file %s.\n',toc,fullfilename);
elseif strcmpi(extension,'.set')
    fprintf('Loading %s.\n',fullfilename);
    EEG = pop_loadset('filename',[filename '.set'],'filepath',filepath);
    fprintf(log_file,'%.2f - Loaded file %s.\n',toc,fullfilename);
    pop_saveset(EEG,'filename',['Original_' filename '.set'],'filepath',[filepath filesep 'Intermediate']);
    EEG.filename = [filename '.set'];
else
    EEG=[];
    fprintf('Unknown file format.\n');
    fprintf(log_file,'%.2f - Unknown file format. Cannot process.\n',toc);
    return;
end
EEG = eeg_checkset(EEG);

% Check if channel locations exist, and if not load them from disk.
if (~isfield(EEG.chanlocs,'X') || ~isfield(EEG.chanlocs,'Y') || ~isfield(EEG.chanlocs,'Z') || isempty(EEG.chanlocs))
    EEG = pop_chanedit(EEG, 'load', {channel_locations_file});
    EEG.saved='no';
    fprintf(log_file,'%.2f - Loaded channel locations file from %s.\n',toc,channel_locations_file);
end
EEG = pop_saveset(EEG,'savemode','resave');

%%%%%%%%%%%%%%%%
% Save options %
%%%%%%%%%%%%%%%%
save_before_filter = save_options(1);
save_before_interp = save_options(2);
save_before_epoch = save_options(3);
save_before_ica_rej = save_options(4);
save_before_epoch_interp = save_options(5);

if save_before_filter
    EEGBAK=EEG;
    EEGBAK.setname = ['pre_filt_' EEG.setname];
    pop_saveset(EEGBAK,'filename',['1_pre_filt_' EEG.filename],'filepath',[filepath filesep 'Intermediate'],'savemode','onefile');
    clear EEGBAK;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter options %
%%%%%%%%%%%%%%%%%%
% 1 Hipass approx -3dB freq (integer > 0) (default: 1)
% 2 Hipass filter order (integer > 0) (default: 3)
% 3 Lopass -3dB freq (integer > 0) (default: 95)
% 4 Lopass filter order (default: 6) (integer > 0) (default: 6)
% 5 Notch -3dB freq (array of 2 integers > 0) (default: [47 53])
% 6 Notch filter order (integer > 0) (default: 2)
% 7 Downsample (integer > 0) (default: off)
%
% If any of the above values are set to zero, the corresponding filter is
% turned off.
%%%%%%%%%%%%%%%%%%
w_h=o.filter_options.hpf_freq;
t_h=o.filter_options.hpf_bandwidth;
r_h=o.filter_options.hpf_ripple;
a_h=o.filter_options.hpf_attenuation;

w_l=o.filter_options.lpf_freq;
t_l=o.filter_options.lpf_bandwidth;
r_l=o.filter_options.lpf_ripple;
a_l=o.filter_options.lpf_attenuation;

w_n=[o.filter_options.notch_freq-o.filter_options.notch_bandwidth1/2 o.filter_options.notch_freq+o.filter_options.notch_bandwidth1/2];
t_n=o.filter_options.notch_bandwidth2;
r_n=o.filter_options.notch_ripple;
a_n=o.filter_options.notch_attenuation;
resample_frequency=o.filter_options.resample_freq;
do_resample=o.filter_options.resample_on;
% Downsampling is done later (shouldn't really be done at all).

do_hipass=o.filter_options.hpf_on;
do_lopass=o.filter_options.lpf_on;
do_notch=o.filter_options.notch_on;

if any(any(isnan(EEG.data)))
    fprintf('NaN in EEG data before filtering.\n');
end

testing=0;
if (~testing)
    if do_hipass
        [m, wtpass, wtstop] = pop_firpmord([w_h-(t_h/2) w_h+(t_h/2)], [0 1], [10^(-1*abs(a_h)/20) (10^(r_h/20)-1)/(10^(r_h/20)+1)], EEG.srate);
        if mod(m,2);m=m+1;end;
        EEG = pop_firpm(EEG, 'fcutoff', w_h, 'ftrans', t_h, 'ftype', 'highpass', 'wtpass', wtpass, 'wtstop', wtstop, 'forder', m);
        EEG.saved='no';
        fprintf(log_file,'%.2f - Highpass filter: %.3fHz, transition band: %.2f, order: %d.\n',toc,w_h,t_h,m);
    end

    if do_lopass
        [m, wtpass, wtstop] = pop_firpmord([w_l-(t_l/2) w_l+(t_l/2)], [1 0], [(10^(r_l/20)-1)/(10^(r_l/20)+1) 10^(-1*abs(a_l)/20)], EEG.srate);
        if mod(m,2);m=m+1;end;
        EEG = pop_firpm(EEG, 'fcutoff', w_l, 'ftrans', t_l, 'ftype', 'lowpass', 'wtpass', wtpass, 'wtstop', wtstop, 'forder', m);
        EEG.saved='no';
        fprintf(log_file,'%.2f - Lowpass filter: %.3fHz, transition band: %.2f, order: %d.\n',toc,w_l,t_l,m);
    end

    if do_notch
        [m, wtpass, wtstop] = pop_firpmord([w_n(1)-(t_n/2) w_n(1)+(t_n/2) w_n(2)-(t_n/2) w_n(2)+(t_n/2)], [0 1 0], [10^(-1*abs(a_n)/20) (10^(r_n/20)-1)/(10^(r_n/20)+1) 10^(-1*abs(a_n)/20)], EEG.srate);
        if mod(m,2);m=m+1;end;
        EEG = pop_firpm(EEG, 'fcutoff', w_n, 'ftrans', t_n, 'ftype', 'bandstop', 'wtpass', wtpass, 'wtstop', wtstop, 'forder', m);
        EEG.saved='no';
        fprintf(log_file,'%.2f - Notch filter: %.3f to %.3fHz, transition band: %.2f, order: %d.\n',toc,w_n(1),w_n(2),t_n,m);
    end
else
    load('butterworth_filters.mat');
    EEG.data=EEG.data(:,:)';
    EEG.data = filtfilthd(Hd_HP,EEG.data(:,:));
    EEG.data = filtfilthd(Hd_LP,EEG.data(:,:));
    EEG.data = filtfilthd(Hd_Notch,EEG.data(:,:));
    EEG.data=EEG.data(:,:)';
end

EEG = pop_saveset(EEG,'savemode','resave');

if save_before_interp
    EEGBAK=EEG;
    EEGBAK.setname = ['pre_interp_' EEG.setname];
    pop_saveset(EEGBAK,'filename',['2_pre_interp_' EEG.filename],'filepath',[filepath filesep 'Intermediate'],'savemode','onefile');
    clear EEGBAK;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data cutoff point      %
% Will be re-implemented %
%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ~isempty(cutoff_markers) && any(cutoff_markers)
% 	cutoff_point=[0 size(EEG.data,2)+1];
% 	for u=1:length(EEG.event)
% 		if EEG.event(u).type == cutoff_markers(1) || strcmp(EEG.event(u).type,cutoff_markers(1))
% 			cutoff_point(1)=EEG.event(u).latency; % Finds the last 255 (check this one)
% 		end
% 		if EEG.event(u).type == cutoff_markers(2) || strcmp(EEG.event(u).type,cutoff_markers(2))
% 			cutoff_point(2)=EEG.event(u).latency; % Finds the last 255 (check this one)
% 		end
% 	end
% 	if cutoff_point(1) > 1
% 		EEG = pop_select( EEG, 'nopoint',[1 cutoff_point(1)] );
% 	end
% 	if cutoff_point(2) < size(EEG.data,2)
% 		EEG = pop_select( EEG, 'nopoint',[cutoff_point(2) size(EEG.data(:,:),2)] );
% 	end
% end

% %New cutoff points for VESPA
%
% EEG = remevent(EEG,768);EEG=remevent(EEG,33536);
%
%
% first_real_event = -1;
% last_real_event = -1;
%
% for u=1:length(EEG.event)-2
%
% 	if ((EEG.event(u).latency - EEG.event(u+1).latency) * (1000/EEG.srate) < 100 && (EEG.event(u+1).latency - EEG.event(u+2).latency) * (1000/EEG.srate) < 100 && first_real_event == -1)
% 		first_real_event = u;
% 	end
%
% 	if (first_real_event ~= -1 && (EEG.event(u).latency - EEG.event(u+1).latency) * (1000/EEG.srate) > 100 && (EEG.event(u+1).latency - EEG.event(u+2).latency) * (1000/EEG.srate) > 100 && last_real_event == -1)
% 		last_real_event = u;
% 	end
%
% end
%
% first_real_time=max(EEG.event(first_real_event).latency - EEG.srate,1);
%
% if (last_real_event==-1)
% 	last_real_time=min(EEG.event(end).latency + EEG.srate,size(EEG.data(:,:),2));
% else
% 	last_real_time=min(EEG.event(last_real_event).latency + EEG.srate,1);
% end
%
% EEG = pop_select( EEG, 'point',[first_real_time:last_real_time] );
% EEG.saved='no';
% EEG = pop_saveset(EEG,'savemode','resave');
%
% fprintf(log_file,'Cropped between %.2f and %.2f seconds.\n',first_real_time/EEG.srate,last_real_time/EEG.srate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel interpolation options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 Automatic interpolation of bad channels on or off (1 / 0)
% 2 Radius for channel interpolation hypersphere (integer > 0)
% 3 Automatic interpolation of channels per single epoch at end of process (1 / 0)
% 4 Radius for epoch interpolation hypersphere (integer > 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_auto_interp = o.channel_options.channel_rejection_on;

if do_auto_interp
    list_properties = channel_properties(EEG,eeg_chans,ref_chan);
    lengths = min_z(list_properties,o.channel_options.rejection_options); % Need to edit to make rejection_options.measure a vector, instead of multiple fields
    chans_to_interp = union(find(lengths),o.channel_options.bad_channels);
    chans_to_interp = setdiff(chans_to_interp,ref_chan); % Ref chan may appear bad, but we shouldn't interpolate it!
    if ~isempty(chans_to_interp)
        fprintf('Interpolating channel(s)');
        fprintf(' %d',chans_to_interp);
        fprintf('.\n');
        EEG = h_eeg_interp_spl(EEG,chans_to_interp,'spherical',ext_chans);
        EEG.saved='no';
        fprintf(log_file,'%.2f - Interpolated channels',toc); fprintf(log_file,' %d',chans_to_interp); fprintf(log_file,'.\n');
    end
end

EEG = pop_saveset(EEG,'savemode','resave');

if save_before_epoch
    EEGBAK=EEG;
    EEGBAK.setname = ['pre_epoch_' EEG.setname];
    pop_saveset(EEGBAK,'filename',['3_pre_epoch_' EEG.filename],'filepath',[filepath filesep 'Intermediate'],'savemode','onefile');
    clear EEGBAK;
end

%%% Do resampling here (if done pre-filtering, it creates problems). %%%
%%% It does anyway, it seems. %%%
if do_resample
    old_name = EEG.setname;
    old_srate = EEG.srate;
    EEG = pop_resample( EEG, resample_frequency);
    EEG.setname = old_name;
    fprintf(log_file,'%.2f - Resampled from %dHz to %dHz.\n',toc,old_srate,resample_frequency);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Epoch options %
%%%%%%%%%%%%%%%%%
% 1 Epoching on or off (1 / 0)
% 2 Markers to epoch from (array of integers or cell of strings)
% 3 Epoch length (vector of 2 floats, 1 negative, 1 positive) - seconds
% 4 Baseline length for mean subtraction (vector of 2 integers) (0 => baseline subtraction off) - milliseconds
% 5 Auto epoch rejection on or off (1 / 0)
% 6 Radius for epoch rejection hypersphere (integer > 0)
%%%%%%%%%%%%%%%%%
do_epoching = ~isempty(o.epoch_options.epoch_markers) && any(o.epoch_options.epoch_limits);
markers = o.epoch_options.epoch_markers;
epoch_length = o.epoch_options.epoch_limits;
baseline_time = o.epoch_options.baseline_sub * 1000;
do_epoch_rejection = o.epoch_options.epoch_rejection_on;

%%%%%%%%%%%%%%
% Epoch data %
%%%%%%%%%%%%%%
if do_epoching
    oldname = EEG.setname;
    EEG = h_epoch(EEG,markers,epoch_length);
    EEG.setname = oldname;
    EEG.saved='no';
    if isnumeric(markers)
        fprintf(log_file,'%.2f - Epoched data on markers',toc);
        fprintf(log_file,' %d',markers);
        fprintf(log_file,'.\n');
    else
        fprintf(log_file,'%.2f - Epoched data on markers',toc);
        fprintf(log_file,' %s',markers{:});
        fprintf(log_file,'.\n');
    end

    % Remove epoch baselines after epoching:
    if ~any(baseline_time)
        EEG = pop_rmbase( EEG, baseline_time);
    end
    
    fprintf(log_file,'Initial baseline variance: %.2f.\n',var(mean(EEG.data(1:round(EEG.srate*-1*EEG.xmin),:,:),3),[],2));
end

EEG = pop_saveset(EEG,'savemode','resave');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Epoch rejection section %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_epoch_rejection
    list_properties = epoch_properties(EEG,eeg_chans);
    [lengths] = min_z(list_properties,o.epoch_options.rejection_options);
    EEG=pop_rejepoch(EEG, find(lengths),0);
    fprintf(log_file,'%.2f - Rejected %d epochs',toc,length(find(lengths)));
    fprintf(log_file,' %d',find(lengths));
    fprintf(log_file,'.\n');
    EEG.saved='no';
end

EEG = pop_saveset(EEG,'savemode','resave');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average reference %
%%%%%%%%%%%%%%%%%%%%%
if (do_reref)
    EEG = h_pop_reref(EEG, [], 'exclude',ext_chans, 'refstate', ref_chan);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ICA options %
%%%%%%%%%%%%%%%
% 1 ICA on or off (1 / 0)
% 2 Auto component rejection on or off (1 / 0)
% 3 Radius for component rejection hypersphere (integer > 0)
% 4 EOG channels (vector of integers)
%%%%%%%%%%%%%%%
do_ica = o.ica_options.run_ica;
k_value = o.ica_options.k_value;
do_component_rejection = o.ica_options.component_rejection_on;
EOG_chans = o.ica_options.EOG_channels;

%%%%%%%%%%
% Do ICA %
%%%%%%%%%%
if do_ica
    num_pca = min(floor(sqrt(size(EEG.data(:,:),2) / k_value)),(size(EEG.data,1) - length(chans_to_interp) - 1));
    EEG = pop_runica(EEG,  'icatype', 'runica', 'dataset',1, 'options',{'extended',1,'pca',num_pca});
    EEG.saved='no';
    fprintf(log_file,'%.2f - Ran ICA.\n',toc);
end

EEG = pop_saveset(EEG,'savemode','resave');

if save_before_ica_rej
    EEGBAK=EEG;
    EEGBAK.setname = ['pre_comp_rej_' EEG.setname];
    pop_saveset(EEGBAK,'filename',['4_pre_comp_rej_' EEG.filename],'filepath',[filepath filesep 'Intermediate'],'savemode','onefile');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Component rejection section %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_component_rejection && ~isempty(EEG.icaweights)
    EEG = eeg_checkset(EEG);
    original_name=EEG.setname;
    list_properties = component_properties(EEG,EOG_chans,[w_l-(t_l/2) w_l+(t_l/2)]);
    [lengths] = min_z(list_properties,o.ica_options.rejection_options);
    fprintf('Rejecting components');
    fprintf(' %d',find(lengths));
    fprintf('.\n');
    EEG = pop_subcomp(EEG, find(lengths), 0);
    fprintf(log_file,'%.2f - Rejected %d components',toc,length(find(lengths)));
    fprintf(log_file,' %d',find(lengths));
    fprintf(log_file,'.\n');
    EEG.setname=original_name;
    EEG.saved='no';
end

EEG = pop_saveset(EEG,'savemode','resave');

if save_before_epoch_interp
    EEGBAK=EEG;
    EEGBAK.setname = ['pre_epoch_interp_' EEG.setname];
    pop_saveset(EEGBAK,'filename',['5_pre_epoch_interp_' EEG.filename],'filepath',[filepath filesep 'Intermediate'],'savemode','onefile');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Epoch interpolation section %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_epoch_interp=o.epoch_interp_options.epoch_interpolation_on;
if do_epoch_interp && length(size(EEG.data)) > 2
    status = '';
    lengths_ep=cell(1,size(EEG.data,3));
    for v=1:size(EEG.data,3)
        list_properties = single_epoch_channel_properties(EEG,v,eeg_chans);
        lengths_ep{v}=find(min_z(list_properties,o.epoch_interp_options.rejection_options));
        status = [status sprintf('%d: ',v) sprintf('%d ',lengths_ep{v}) sprintf('\n')];
    end
    EEG=h_epoch_interp_spl(EEG,lengths_ep,ext_chans);
    EEG.saved='no';
    epoch_interps_log_file=fopen([filepath filesep filename(1:end-4) '_epoch_interpolations.txt'],'a');
    fprintf(epoch_interps_log_file,'%s',status);
    fclose(epoch_interps_log_file);
    fprintf(log_file,'%.2f - Did per-epoch interpolation cleanup.\n',toc);
    fprintf(log_file,['See ' filename(1:end-4) '_epoch_interpolations.txt for details.\n']);
end
EEG = pop_saveset(EEG,'savemode','resave');

fprintf('Done with file %s.\nTook %d seconds.\n',[filepath filesep filename],toc);

fprintf(log_file,'%.2f - Finished.\n',toc);
if (size(EEG.data,3>1))
    fprintf(log_file,'Final baseline variance: %.2f.\n',var(mean(EEG.data(1:round(EEG.srate*-1*EEG.xmin),:,:),3),[],2));
    % More stats here!
end
fclose(log_file);
end