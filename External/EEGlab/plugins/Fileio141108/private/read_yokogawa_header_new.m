function hdr = read_yokogawa_header_new(filename)


% READ_YOKOGAWA_HEADER_NEW reads the header information from continuous,
% epoched or averaged MEG data that has been generated by the Yokogawa
% MEG system and software and allows that data to be used in combination
% with FieldTrip.
%
% Use as
%  [hdr] = read_yokogawa_header_new(filename)
%
% This is a wrapper function around the functions
% getYkgwHdrSystem
% getYkgwHdrChannel
% getYkgwHdrAcqCond
% getYkgwHdrCoregist
% getYkgwHdrDigitize
% getYkgwHdrSource
%
% See also READ_YOKOGAWA_DATA_NEW, READ_YOKOGAWA_EVENT

% ** 
% Copyright (C) 2005, Robert Oostenveld and 2010, Tilmann Sander-Thoemmes
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: read_yokogawa_header_new.m 7123 2012-12-06 21:21:38Z roboos $

% FIXED
%  txt -> m
%  fopen iee-le

if ~ft_hastoolbox('yokogawa_meg_reader')
    error('cannot determine whether Yokogawa toolbox is present');
end

handles = definehandles;

sys_info = getYkgwHdrSystem(filename);
id = sys_info.system_id;
ver = sys_info.version;
rev = sys_info.revision;
sys_name = sys_info.system_name;
model_name = sys_info.model_name;
clear('sys_info'); % remove structure as local variables are collected in the end

channel_info = getYkgwHdrChannel(filename);
channel_count = channel_info.channel_count;

acq_cond = getYkgwHdrAcqCond(filename);
acq_type = acq_cond.acq_type;

% these depend on the data type
sample_rate        = [];
sample_count       = [];
pretrigger_length  = [];
averaged_count     = [];
actual_epoch_count = [];

switch acq_type
  case handles.AcqTypeContinuousRaw
    sample_rate = acq_cond.sample_rate;
    sample_count = acq_cond.sample_count;
    if isempty(sample_rate) | isempty(sample_count)
      error('invalid sample rate or sample count in ', filename);
      return;
    end
    pretrigger_length = 0;
    averaged_count = 1;

  case handles.AcqTypeEvokedAve
    sample_rate = acq_cond.sample_rate;
    sample_count = acq_cond.frame_length;
    pretrigger_length = acq_cond.pretrigger_length;
    averaged_count = acq_cond.average_count; 
    if isempty(sample_rate) | isempty(sample_count) | isempty(pretrigger_length) | isempty(averaged_count)
      error('invalid sample rate or sample count or pretrigger length or average count in ', filename);
      return;
    end
    if acq_cond.multi_trigger.enable 
      error('multi trigger mode not supported for ', filename);
      return;
    end
  case handles.AcqTypeEvokedRaw
    sample_rate = acq_cond.sample_rate;
    sample_count = acq_cond.frame_length;
    pretrigger_length = acq_cond.pretrigger_length;
    actual_epoch_count = acq_cond.average_count; 
    if isempty(sample_rate) | isempty(sample_count) | isempty(pretrigger_length) | isempty(actual_epoch_count)
      error('invalid sample rate or sample count or pretrigger length or epoch count in ', filename);
      return;
    end
    if acq_cond.multi_trigger.enable 
      error('multi trigger mode not supported for ', filename);
      return;
    end

  otherwise
    error('unknown data type');
end
clear('acq_cond'); % remove structure as local variables are collected in the end

coregist = getYkgwHdrCoregist(filename);
digitize = getYkgwHdrDigitize(filename);
source = getYkgwHdrSource(filename);

% put all local variables into a structure, this is a bit unusual matlab programming style
tmp = whos;
orig = [];
for i=1:length(tmp)
  if isempty(strmatch(tmp(i).name, {'tmp', 'ans', 'handles'}))
    orig = setfield(orig, tmp(i).name, eval(tmp(i).name));
  end
end

% convert the original header information into something that FieldTrip understands
hdr = [];
hdr.orig         = orig;                % also store the original full header information
hdr.Fs           = orig.sample_rate;    % sampling frequency
hdr.nChans       = orig.channel_count;  % number of channels
hdr.nSamples     = [];                  % number of samples per trial
hdr.nSamplesPre  = [];                  % number of pre-trigger samples in each trial
hdr.nTrials      = [];                  % number of trials

switch orig.acq_type
  case handles.AcqTypeEvokedAve
    hdr.nSamples    = orig.sample_count;
    hdr.nSamplesPre = orig.pretrigger_length;
    hdr.nTrials     = 1;                % only the average, which can be considered as a single trial
  case handles.AcqTypeContinuousRaw
    hdr.nSamples    = orig.sample_count;
    hdr.nSamplesPre = 0;                % there is no fixed relation between triggers and data
    hdr.nTrials     = 1;                % the continuous data can be considered as a single very long trial
  case handles.AcqTypeEvokedRaw
    hdr.nSamples    = orig.sample_count;
    hdr.nSamplesPre = orig.pretrigger_length;
    hdr.nTrials     = orig.actual_epoch_count;
  otherwise
    error('unknown acquisition type');
end

% construct a cell-array with labels of each channel
for i=1:hdr.nChans
% this should be consistent with the predefined list in ft_senslabel,
% with yokogawa2grad_new and with ft_channelselection
  if     hdr.orig.channel_info.channel(i).type == handles.NullChannel
    prefix = '';
  elseif hdr.orig.channel_info.channel(i).type == handles.MagnetoMeter
    prefix = 'M';
  elseif hdr.orig.channel_info.channel(i).type == handles.AxialGradioMeter
    prefix = 'AG';
  elseif hdr.orig.channel_info.channel(i).type == handles.PlannerGradioMeter
    prefix = 'PG';
  elseif hdr.orig.channel_info.channel(i).type == handles.RefferenceMagnetoMeter
    prefix = 'RM';
  elseif hdr.orig.channel_info.channel(i).type == handles.RefferenceAxialGradioMeter
    prefix = 'RAG';
  elseif hdr.orig.channel_info.channel(i).type == handles.RefferencePlannerGradioMeter
    prefix = 'RPG';
  elseif hdr.orig.channel_info.channel(i).type == handles.TriggerChannel
    prefix = 'TRIG';
  elseif hdr.orig.channel_info.channel(i).type == handles.EegChannel
    prefix = 'EEG';
  elseif hdr.orig.channel_info.channel(i).type == handles.EcgChannel
    prefix = 'ECG';
  elseif hdr.orig.channel_info.channel(i).type == handles.EtcChannel
    prefix = 'ETC';
  end
  hdr.label{i} = sprintf('%s%03d', prefix, i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this defines some usefull constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = definehandles;
handles.output = [];
handles.sqd_load_flag = false;
handles.mri_load_flag = false;
handles.NullChannel         = 0;
handles.MagnetoMeter        = 1;
handles.AxialGradioMeter    = 2;
handles.PlannerGradioMeter  = 3;
handles.RefferenceChannelMark = hex2dec('0100');
handles.RefferenceMagnetoMeter       = bitor( handles.RefferenceChannelMark, handles.MagnetoMeter );
handles.RefferenceAxialGradioMeter   = bitor( handles.RefferenceChannelMark, handles.AxialGradioMeter );
handles.RefferencePlannerGradioMeter = bitor( handles.RefferenceChannelMark, handles.PlannerGradioMeter );
handles.TriggerChannel      = -1;
handles.EegChannel          = -2;
handles.EcgChannel          = -3;
handles.EtcChannel          = -4;
handles.NonMegChannelNameLength = 32;
handles.DefaultMagnetometerSize       = (4.0/1000.0);       % Square of 4.0mm in length
handles.DefaultAxialGradioMeterSize   = (15.5/1000.0);      % Circle of 15.5mm in diameter
handles.DefaultPlannerGradioMeterSize = (12.0/1000.0);      % Square of 12.0mm in length
handles.AcqTypeContinuousRaw = 1;
handles.AcqTypeEvokedAve     = 2;
handles.AcqTypeEvokedRaw     = 3;
handles.sqd = [];
handles.sqd.selected_start  = [];
handles.sqd.selected_end    = [];
handles.sqd.axialgradiometer_ch_no      = [];
handles.sqd.axialgradiometer_ch_info    = [];
handles.sqd.axialgradiometer_data       = [];
handles.sqd.plannergradiometer_ch_no    = [];
handles.sqd.plannergradiometer_ch_info  = [];
handles.sqd.plannergradiometer_data     = [];
handles.sqd.eegchannel_ch_no   = [];
handles.sqd.eegchannel_data    = [];
handles.sqd.nullchannel_ch_no   = [];
handles.sqd.nullchannel_data    = [];
handles.sqd.selected_time       = [];
handles.sqd.sample_rate         = [];
handles.sqd.sample_count        = [];
handles.sqd.pretrigger_length   = [];
handles.sqd.matching_info   = [];
handles.sqd.source_info     = [];
handles.sqd.mri_info        = [];
handles.mri                 = [];
