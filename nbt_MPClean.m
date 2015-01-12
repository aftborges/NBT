% nbt_MPClean(Signal,SignalInfo)
%
%
%
% Usage:
%
%
% Inputs:
% ICAswitch     : -1 for "extended ICA with no pca", 0 for "automatic pca
%                reduction", any other number = number of pca reduced ICA components
% NonEEGCh      : list of Non-EEG channels
% EyeCh         : list of eye-channels (used for cleaning eye artifacts)
% ResampleFS    : to resample set this to the resampling frequency.
%
% Outputs:
%
% Example:
%
%
% References:
%
% See also

function [Signal, SignalInfo] = nbt_MPClean(Signal, SignalInfo, SignalPath, ICAswitch, NonEEGCh, EyeCh)
narginchk(3,6) % took out the Resample part
if(~isempty(NonEEGCh))
    SignalInfo.NonEEGch = NonEEGCh;
end
if(~isempty(EyeCh))
    SignalInfo.EyeCh    = EyeCh;
end

if(isempty(SignalInfo.NonEEGch))
    SignalInfo.NonEEGch = input('Please specify Non-EEG channels: ');
end

if(isempty(SignalInfo.EyeCh))
    SignalInfo.EyeCh = input('Please specify eye channels: ');
end

NonEEGCh = SignalInfo.NonEEGch;
EyeCh = SignalInfo.EyeCh;


% Possivelmente n eh NECESSARIO AQUI
%. 0. Ref-ref to Cz
%first we find Cz
cznotfound = true;
% for CzID = 1:SignalInfo.Interface.number_of_channels
%     if(strcmpi(SignalInfo.Interface.EEG.chanlocs(CzID).labels,'Cz'))
%         cznotfound = false;
%         break;
%     end
% end
if(cznotfound)
    CzID = 129;
    %CzID = input('Please specify Cz channel number')
end



% 1. Filter Data
%Let's not bandpass filter the data for now! - check functioning of FASTER
%and ADJUST without it;
%[Signal] = nbt_filter_fir(Signal,0.5,45,SignalInfo.converted_sample_frequency,2/0.5,1);

% 2. Mark Bad Channels
[Signal, SignalInfo] = nbt_EEGLABwrp(@nbt_FindBadChannels, Signal, SignalInfo, [] , 0, 'f', NonEEGCh);
%decidir qual opcao  melhor nbt_findBadChannels 'f' - FASTER
SignalInfo.BadChannels(NonEEGCh) = 1;
% 3. Reject Transient artifacts -- Nao quero fazer epoching aqui
% provavelmente --remove muito esta funcao cerca de 25% surreal (mas no MFF) 
 [Signal, SignalInfo] = nbt_AutoRejectTransient(Signal,SignalInfo,NonEEGCh);
% % 4. Run ICA
% %Re-reference to Cz - because autoreject ICA expects Cz referenced
% %topomaps.
 %[Signal, SignalInfo] = nbt_EEGLABwrp(@nbt_ReRef, Signal,SignalInfo,[],0,CzID);
 [Signal, SignalInfo] = nbt_EEGLABwrp(@nbt_filterbeforeICA, Signal,SignalInfo, [], 0, '',1,ICAswitch); % I am going to decrease the offset
% to 1 because we already remove the first 5 seconds (1 instead of 4)
% % 5. Reject ICA compoents
 [Signal, SignalInfo] = nbt_EEGLABwrp(@nbt_AutoRejectICA,Signal, SignalInfo, [],0, EyeCh,0);
% % 6. Average Ref
 [Signal, SignalInfo] = nbt_EEGLABwrp(@nbt_ReRef,Signal, SignalInfo, [],0,[]);
nbt_SaveSignal(Signal, SignalInfo, SignalPath,1,'AutoICASignal')
end