
% nbt_runAmplitudeHF - Run computation of HF amplitude among channels and save the biomarker in the analysis.mat file; 
%
% Usage:
%   nbt_runAmplitudeHF(Signal, SignalInfo, SaveDir)
%
% Inputs:
%   Signal
%   SignalInfo
%   SaveDir - directory where you want to save the analysis file
%
% Outputs: 
%
% Example:
%
% References:
% 
% See also: 
%  
  
%------------------------------------------------------------------------------------
% Originally created by Rick Jansen(2012), see NBT website (http://www.nbtwiki.net) for current email address
%------------------------------------------------------------------------------------
%
% ChangeLog - see version control log at NBT website for details.
%
% Copyright (C) <year>  <Main Author>  (Neuronal Oscillations and Cognition group, 
% Department of Integrative Neurophysiology, Center for Neurogenomics and Cognitive Research, 
% Neuroscience Campus Amsterdam, VU University Amsterdam)
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
% -------------------------------------------------------------------------
function nbt_runAmplitudewHF(Signal, SignalInfo, SaveDir)
% compute biomarkers
% [amplitude_1_4_Hz amplitude_4_8_Hz amplitude_8_13_Hz amplitude_13_30_Hz amplitude_30_45_Hz ... 
%  amplitude_1_4_Hz_Normalized amplitude_4_8_Hz_Normalized amplitude_8_13_Hz_Normalized  ...
%  amplitude_13_30_Hz_Normalized amplitude_30_45_Hz_Normalized,amplitude_1_4_Hz_cmbgrad,amplitude_4_8_Hz_cmbgrad, ...
%     amplitude_8_13_Hz_cmbgrad,amplitude_13_30_Hz_cmbgrad,amplitude_30_45_Hz_cmbgrad,amplitude_1_4_Hz_cmbgradN, ...
%     amplitude_4_8_Hz_cmbgradN,amplitude_8_13_Hz_cmbgradN,amplitude_13_30_Hz_cmbgradN,amplitude_30_45_Hz_cmbgradN] = nbt_doAmplitude(Signal,SignalInfo);

[amplitude_1_4_Hz, ...
    amplitude_4_8_Hz,amplitude_8_13_Hz,amplitude_13_30_Hz,amplitude_30_45_Hz,...
    amplitude_55_125_Hz] = nbt_doAmplitudewHF(Signal,SignalInfo);






% save biomarkers


nbt_SaveClearObject('amplitude_1_4_Hz',SignalInfo,SaveDir);
nbt_SaveClearObject('amplitude_4_8_Hz',SignalInfo,SaveDir);
nbt_SaveClearObject('amplitude_8_13_Hz',SignalInfo,SaveDir);
nbt_SaveClearObject('amplitude_13_30_Hz',SignalInfo,SaveDir);
nbt_SaveClearObject('amplitude_30_45_Hz',SignalInfo,SaveDir);
nbt_SaveClearObject('amplitude_55_125_Hz' ,SignalInfo,SaveDir);

