classdef nbt_rsq < nbt_Biomarker
    %% Copyright (c) 2010,  B.A. Diaz (Center for Neurogenomics and Cognitive Research (CNCR), VU University Amsterdam)
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
    %% ChangeLog - remember to set NBTversion property
    %$ Version 1.3 - May 30th 2010 : Modified by B.A. Diaz, alexander.diazATcncrDOTvuDOTnl$
    % Implementing new matlab object structure.
    %%
    properties
        Questions
        Answers
    end
    methods
        function RSQobject = nbt_rsq(Questions, Answers)
            if nargin == 0
               Questions = [];
               Answers = [];
            end
            RSQobject.Questions = Questions;
            RSQobject.Answers = Answers;
        end
    end
end