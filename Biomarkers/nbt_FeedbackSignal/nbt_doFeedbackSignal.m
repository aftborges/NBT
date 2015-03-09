function feedbackBiomarker = nbt_doFeedbackSignal(Signal,Info,Save_dir)

disp(' ')
disp('Command window code:')
disp(['nbt_get_amplitude_freq_bands(Signal,SignalInfo,',char(39),Save_dir,char(39),')'])
disp(' ')

disp(['Computing amplitudes for ',Info.file_name])

%% Calculate feedbackSignals

[ECRfname,ECRpname] = uigetfile('*ECR*mat',['Select the ECR session for file ',Info.file_name]);
ECRSignal = load ([ECRpname,ECRfname],'RawSignal');
ECRSignal = ECRSignal.RawSignal;
FeedbackSignal = Signal;

ECRSignal = ECRSignal * 0.0240;
FeedbackSignal = FeedbackSignal * 0.0240;
spatialFilter = load('AlphaICAWeights.mat');
spatialFilter = spatialFilter.AlphaICAWeights;
icaFBD = FeedbackSignal * spatialFilter(:,3);
icaECR = ECRSignal * spatialFilter(:,3);

%% Calculate weights
k = 0;
fs = 250;
for i = 1:fs:length(icaECR)-(3*fs)
    k = k+1;
    [p,f] = pwelch(icaECR(i:i+(3*fs)),hamming(2^9),[],2^9,fs);
    powsECR(k) = sum(p(f>=8 & f<=13));
end


k = 0;
for i = 1:fs:length(icaFBD)-(3*fs)
    k = k+1;
    [p,f] = pwelch(icaFBD(i:i+(3*fs)),hamming(2^9),[],2^9,fs);
    powsFBD(k) = sum(p(f>=8 & f<=13));
end

%% Calculate range based on feedback
sortedSignal = sort(powsECR);
lengthSignal = length(powsECR);
min5 = log10(sortedSignal(floor(lengthSignal *0.05)));
max5 = log10(sortedSignal(floor(lengthSignal *0.95)));
rangeSignal = (max5-min5)/3;
min5 = min5 - rangeSignal;
max5 = max5 + rangeSignal;
med = min5+((max5-min5)/2);


powsECR = 100*((log10(powsECR) - min5)/ (max5-min5));
powsFBD = 100*((log10(powsFBD) - min5)/ (max5-min5));

% figure
% plot(powsECR,'b')
% hold on
% plot(powsFBD,'r')

ECRbelow50 = 100*nnz(powsECR<50)/length(powsECR);
ECRabove50 = 100*nnz(powsECR>50)/length(powsECR);
FBDbelow50 = 100*nnz(powsFBD<50)/length(powsFBD);
FBDabove50 = 100*nnz(powsFBD>50)/length(powsFBD);
%% create biomarker and store

feedbackBiomarker = nbt_FeedbackSignal(powsFBD,powsECR,FBDbelow50,FBDabove50,ECRbelow50,ECRabove50);
BiomarkerTemplate = nbt_UpdateBiomarkerInfo(feedbackBiomarker, Info);

end