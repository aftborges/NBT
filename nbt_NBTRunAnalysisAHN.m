% Copyright (C) 2012  Neuronal Oscillations and Cognition group, Department of Integrative Neurophysiology, Center for Neurogenomics and Cognitive Research, Neuroscience Campus Amsterdam, VU University Amsterdam.
%
% Part of the Neurophysiological Biomarker Toolbox (NBT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
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
%
% See Readme.txt for additional copyright information.
%

% ChangeLog - see version control log for details
% <date> - Version <#> - <text>

function nbt_NBTRunAnalysisAll(varargin)
script = NBTwrapper();
nbt_NBTcompute(script,'AutoICASignal',pwd,pwd)
end


function NBTfunction_handle = NBTwrapper()

% delta_hp = 1;
% delta_lp = 4;
% delta_fo = 2;
% 
% theta_hp = 4;
% theta_lp = 8;
% theta_fo = 0.5;
% 
% alpha_hp = 8;
% alpha_lp = 13;
% alpha_fo = 0.26;
% 
% beta_hp = 13;
% beta_lp = 30;
% beta_fo = 0.16;
% 
% gamma_hp = 30;
% gamma_lp = 45;
% gamma_fo = 2/30;
% 
% highGamma_hp = 55;
% highGamma_lp = 125;
% highGamma_fo = 2/55;
% 
% 
% DFAshortFit = 1;
% DFAlongFit = 20;
% DFAshortCalc = 0.1;
% DFAlongCalc = 80;
% DFA_Overlap = 0.5;
% DFA_Plot = 0;
% ChannelToPlot = 1;

    function NBTscript(Signal, SignalInfo, SaveDir)
       % [Signal, SignalInfo] = nbt_EEGLABwrp(@pop_resample, Signal, SignalInfo, [], 0, 200);
      %  SignalInfo.converted_sample_frequency = 200;
        
%         Hjorthpara = nbt_doHjorthPara(Signal, SignalInfo,0,0);
%         nbt_SaveClearObject('Hjorthpara', SignalInfo, SaveDir)
%         
%         Barlow = nbt_doBarlow(Signal,SignalInfo,0);
%         nbt_SaveClearObject('Barlow', SignalInfo, SaveDir)
%         
%         PeakFit = nbt_doPeakFit(Signal, SignalInfo);
%         nbt_SaveClearObject('PeakFit', SignalInfo, SaveDir)
%         
%         SignalInfo.frequencyRange = [1 45];
%         
%         AmpDist_signal = nbt_doAmpDist(Signal, SignalInfo);
%         nbt_SaveClearObject('AmpDist_signal', SignalInfo, SaveDir)

%% Delta
%         freqstab_delta=nbt_doFreqStability( Signal, SignalInfo, delta_hp, delta_lp);
%         nbt_SaveClearObject('freqstab_delta', SignalInfo, SaveDir)
        
%         [AmplitudeEnvelope, AmplitudeEnvelopeInfo] = nbt_GetAmplitudeEnvelope(Signal,SignalInfo,delta_hp,delta_lp,delta_fo);
%         
%         AmpDist_delta = nbt_doAmpDist(AmplitudeEnvelope, AmplitudeEnvelopeInfo);
%         nbt_SaveClearObject('AmpDist_delta', SignalInfo, SaveDir)
%         
%         OscB_delta =  nbt_doOscBursts(AmplitudeEnvelope, AmplitudeEnvelopeInfo, 1, 0.95, 0, 0, 0,[],1);
%         nbt_SaveClearObject('OscB_delta', SignalInfo, SaveDir)
%         
%         AmpCorr_delta = nbt_doAmplitudeCorr(AmplitudeEnvelope,AmplitudeEnvelopeInfo);
%         nbt_SaveClearObject('AmpCorr_delta', SignalInfo, SaveDir)
%         
%         PhaseLock_delta = nbt_doPhaseLocking(Signal,SignalInfo,[delta_hp delta_lp],[0 100],2/delta_hp,[],[],[1 1]);
%         nbt_SaveClearObject('PhaseLock_delta',SignalInfo, SaveDir)
        
%         [DFA_delta,DFA_exp] = nbt_doDFA(AmplitudeEnvelope, AmplitudeEnvelopeInfo, [DFAshortFit DFAlongFit], [DFAshortCalc DFAlongCalc], DFA_Overlap, DFA_Plot, ChannelToPlot, []);
%         nbt_SaveClearObject('DFA_delta', SignalInfo, SaveDir)
% %% theta frequency
%         freqstab_theta=nbt_doFreqStability(Signal, SignalInfo, theta_hp, theta_lp);
%         nbt_SaveClearObject('freqstab_theta', SignalInfo, SaveDir)
%      
%          [AmplitudeEnvelope, AmplitudeEnvelopeInfo] = nbt_GetAmplitudeEnvelope(Signal,SignalInfo,theta_hp,theta_lp,theta_fo);
% %      
%         AmpDist_theta = nbt_doAmpDist(AmplitudeEnvelope, AmplitudeEnvelopeInfo);
%         nbt_SaveClearObject('AmpDist_theta', SignalInfo, SaveDir)
%      
%         OscB_theta =  nbt_doOscBursts( AmplitudeEnvelope, AmplitudeEnvelopeInfo, 1, 0.95, 0, 0, 0,[],1);
%         nbt_SaveClearObject('OscB_theta', SignalInfo, SaveDir)
%         
%         AmpCorr_theta = nbt_doAmplitudeCorr(AmplitudeEnvelope,AmplitudeEnvelopeInfo);
%         nbt_SaveClearObject('AmpCorr_theta', SignalInfo, SaveDir)
%         
%         PhaseLock_theta = nbt_doPhaseLocking(Signal,SignalInfo,[theta_hp theta_lp],[0 100],2/theta_hp,[],[],[1 1]);
%         nbt_SaveClearObject('PhaseLock_theta',SignalInfo, SaveDir)
        
%         [DFA_theta,DFA_exp] = nbt_doDFA(AmplitudeEnvelope, AmplitudeEnvelopeInfo, [DFAshortFit DFAlongFit], [DFAshortCalc DFAlongCalc], DFA_Overlap, DFA_Plot, ChannelToPlot, []);
%         nbt_SaveClearObject('DFA_theta', SignalInfo, SaveDir)
% %         
% %% Alpha       
%         freqstab_alpha=nbt_doFreqStability(Signal, SignalInfo, alpha_hp, alpha_lp);
%         nbt_SaveClearObject('freqstab_alpha', SignalInfo, SaveDir)
%         
         %[AmplitudeEnvelope, AmplitudeEnvelopeInfo]= nbt_GetAmplitudeEnvelope(Signal,SignalInfo,alpha_hp,alpha_lp,alpha_fo);
%         
%         AmpDist_alpha = nbt_doAmpDist(AmplitudeEnvelope, AmplitudeEnvelopeInfo);
%         nbt_SaveClearObject('AmpDist_alpha', SignalInfo, SaveDir)
%        
%         OscB_alpha =  nbt_doOscBursts(AmplitudeEnvelope, AmplitudeEnvelopeInfo, 1, 0.95, 0, 0, 0,[],1);
%         nbt_SaveClearObject('OscB_alpha', SignalInfo, SaveDir)
%        
%         AmpCorr_alpha = nbt_doAmplitudeCorr(AmplitudeEnvelope,AmplitudeEnvelopeInfo);
%         nbt_SaveClearObject('AmpCorr_alpha', SignalInfo, SaveDir)
%         
%         PhaseLock_alpha = nbt_doPhaseLocking(Signal,SignalInfo,[alpha_hp alpha_lp],[0 100],2/alpha_hp,[],[],[1 1]);
%         nbt_SaveClearObject('PhaseLock_alpha',SignalInfo, SaveDir)
        
        %[DFA_alpha,DFA_exp] = nbt_doDFA(AmplitudeEnvelope, AmplitudeEnvelopeInfo, [DFAshortFit DFAlongFit], [DFAshortCalc DFAlongCalc], DFA_Overlap, DFA_Plot, ChannelToPlot, []);
%         nbt_SaveClearObject('DFA_alpha', SignalInfo, SaveDir)
%         
% %% Beta frequency
%         freqstab_beta = nbt_doFreqStability(Signal, SignalInfo, beta_hp, beta_lp);
%         nbt_SaveClearObject('freqstab_beta', SignalInfo, SaveDir)
% 
         %[AmplitudeEnvelope, AmplitudeEnvelopeInfo] = nbt_GetAmplitudeEnvelope(Signal,SignalInfo,beta_hp,beta_lp,beta_fo);
%         
%         AmpDist_beta = nbt_doAmpDist(AmplitudeEnvelope, AmplitudeEnvelopeInfo);
%         nbt_SaveClearObject('AmpDist_beta', SignalInfo, SaveDir)
%         
%         OscB_beta =  nbt_doOscBursts( AmplitudeEnvelope, AmplitudeEnvelopeInfo, 1, 0.95, 0, 0, 0,[],1);
%         nbt_SaveClearObject('OscB_beta', SignalInfo, SaveDir)
%         
%         AmpCorr_beta = nbt_doAmplitudeCorr( AmplitudeEnvelope,AmplitudeEnvelopeInfo);
%         nbt_SaveClearObject('AmpCorr_beta', SignalInfo, SaveDir)
%         
%         PhaseLock_beta = nbt_doPhaseLocking(Signal,SignalInfo,[beta_hp beta_lp],[0 100],2/beta_hp,[],[],[1 1]);
%         nbt_SaveClearObject('PhaseLock_beta',SignalInfo, SaveDir)
        
%         [DFA_beta,DFA_exp] = nbt_doDFA( AmplitudeEnvelope, AmplitudeEnvelopeInfo, [DFAshortFit DFAlongFit], [DFAshortCalc DFAlongCalc], DFA_Overlap, DFA_Plot, ChannelToPlot, []);
%         nbt_SaveClearObject('DFA_beta', SignalInfo, SaveDir,1)  
%         
%         [AmplitudeEnvelope, AmplitudeEnvelopeInfo] = nbt_GetAmplitudeEnvelope(Signal,SignalInfo,gamma_hp,gamma_lp,gamma_fo);
%          [DFA_gamma,DFA_exp] = nbt_doDFA( AmplitudeEnvelope, AmplitudeEnvelopeInfo, [DFAshortFit DFAlongFit], [DFAshortCalc DFAlongCalc], DFA_Overlap, DFA_Plot, ChannelToPlot, []);
%         nbt_SaveClearObject('DFA_gamma', SignalInfo, SaveDir,1)  
%         
%         [AmplitudeEnvelope, AmplitudeEnvelopeInfo] = nbt_GetAmplitudeEnvelope(Signal,SignalInfo,highGamma_hp,highGamma_lp,highGamma_fo);
%          [DFA_highGamma,DFA_exp] = nbt_doDFA( AmplitudeEnvelope, AmplitudeEnvelopeInfo, [DFAshortFit DFAlongFit], [DFAshortCalc DFAlongCalc], DFA_Overlap, DFA_Plot, ChannelToPlot, []);
%         nbt_SaveClearObject('DFA_highGamma', SignalInfo, SaveDir,1)  
        
        nbt_importARSQStudentAHN(SignalInfo.file_name,SignalInfo,SaveDir)

    end

NBTfunction_handle = @NBTscript;
end