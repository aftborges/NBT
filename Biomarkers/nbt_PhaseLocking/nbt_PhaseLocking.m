% nbt_PhaseLocking - Creates a Phase Locking biomarker object 
%
% Usage:
%   BiomarkerObject = nbt_PhaseLocking
%   or
%   BiomarkerObject = nbt_PhaseLocking(NumChannels)
% Inputs:
%   NumChannels
%
% Outputs:
%   PhaseLocking BiomarkerObject    
%
% Example:
%   
% References:
% 
% See also: 
%   nbt_CrossPhaseLocking
%  
  
%------------------------------------------------------------------------------------
% Originally created by Giuseppina Schiavone (2011), see NBT website (http://www.nbtwiki.net) for current email address
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

    
classdef nbt_PhaseLocking < nbt_Biomarker  
    properties 
        
    Ratio 
    PLV 
    Instphase
%     frequencyRange 
    filterorder 
    interval 
%     synchlag
    IndexE %index based on the Shannon entropy
    IndexF %based on the intensity of the first Fourier mode of the distribution
    IndexCP %based on the conditional probability
%     PLV_in_time
%     time_int
        
    end
    methods
       
        function BiomarkerObject = nbt_PhaseLocking(LengthSign,NumChannels)%,overlap,nw)
            if nargin == 0
                LengthSign = 1;
                NumChannels = 1;
            end
%                 overlap = 1;
%                 nw = 1;
%             elseif nargin == 1
%                 NumChannels = 1;
%                 overlap = 1;
%                 nw = 1;
%             elseif nargin == 2
%                 overlap = 1;
%                 nw = 1;
%             elseif nargin == 2
%                 overlap = 1;
%                 nw = 1;
%             end
            % assign values for this biomarker object:
            %% Define Phase Locking values
            BiomarkerObject.Ratio = nan(NumChannels,NumChannels);
            BiomarkerObject.PLV = nan(NumChannels,NumChannels);
            BiomarkerObject.filterorder =  nan(1);
            BiomarkerObject.interval =  nan(1,2); 
%             BiomarkerObject.synchlag =  nan(2,NumChannels,NumChannels); 
            BiomarkerObject.IndexE = nan(NumChannels,NumChannels); %index based on the Shannon entropy
            BiomarkerObject.IndexCP = nan(NumChannels,NumChannels);%based on the conditional probability
            BiomarkerObject.IndexF = nan(NumChannels,NumChannels);%based on the intensity of the first Fourier mode of the distribution     
            %% Define fields for additional information
            BiomarkerObject.DateLastUpdate = datestr(now);
% %             PhaseLogkingObject.PLV_in_time = nan(NumChannels,NumChannels,floor((LengthSign-nw)/(nw/overlap)));
% %             PhaseLogkingObject.time_int = nan(floor((LengthSign-nw)/(nw/overlap)),1);
            PhaseLogkingObject.Instphase = nan(LengthSign,NumChannels);
            BiomarkerObject.PrimaryBiomarker = 'PLV';
            BiomarkerObject.Biomarkers = {'PLV','Instphase'};
           
            
        end
        function plotPLV(obj)
            figure
            pcolor(obj.PLV)
            view(90,-90)
            colorbar
            xlabel('Channels')
            ylabel('Channels')
            axis tight
        end
    end

end

