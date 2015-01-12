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
nbt_NBTcompute(script,'RawSignal',pwd,pwd)
end


function NBTfunction_handle = NBTwrapper()
    function NBTscript(Signal, SignalInfo, SignalPath)        
        %nbt_importARSQStudentAHN(SignalInfo.file_name, SignalInfo, SaveDir)
        
        nbt_MPClean(Signal,SignalInfo,SignalPath,15,130:136,125:126);
        %nbt_epochingMP14(SignalInfo.file_name, SignalInfo, Signal)
        
        
    end

NBTfunction_handle = @NBTscript;
end


