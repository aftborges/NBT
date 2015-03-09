% nbt_run_stat2_noChansBiom - computes the select statistics for two
% groups or conditions and provides visualization of the results
%
% Usage:
%  [s] =
%  nbt_run_stat2_noChansBiom(Group1,Group2,B1,B2,s,biom,unit)
%
% Inputs:
%   Group1,
%   Group2,
%   B1,
%   B2
%   s structure containing statistics information
%   biom
%   regions
%   unit
%
% Outputs:
%  s updated version of the input structure
%
% Example:
%  
%
% References:
% 
% See also: 
%  nbt_plot_2conditions_topo
  
%------------------------------------------------------------------------------------
% Originally created by Giuseppina Schiavone (2012), see NBT website (http://www.nbtwiki.net) for current email address
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
% --------------

function  [s] = nbt_run_stat2_noChansBiom(Group1,Group2,B1,B2,s,biom,unit)
statistic=s.statistic;
statfunc =s.statfunc;
statfuncname=s.statfuncname;
statname=s.statname;
nchans_o_nregs = size(B1,1);
s.biom_name = biom;

if strcmp(char(statfunc),'nanmedian')
    warning('This test is not design for multiple groups');
elseif strcmp(char(statfunc),'lillietest')
    warning('This test is not design for multiple groups');
elseif strcmp(char(statfunc),'swtest')
    warning('This test is not design for multiple groups');
elseif strcmp(char(statfunc),'ttest')
    try
    for i = 1:nchans_o_nregs
        [h,p(i),C(i,:),stats] = ttest(B1(i,:),B2(i,:));
        statvalues(i) = stats.tstat;
    end    
    s = plot_group(Group1,Group2,B1,B2,C,p,s,biom,[],unit);
    catch
        warning('The two groups do not have the same number of subjects')
    end
elseif strcmp(char(statfunc),'signrank')
    try
    B = B2-B1;
    for i = 1:nchans_o_nregs
        [p(i),h,stats] = signrank(B1(i,:),B2(i,:));
        C(i,:)=bootci(1000,statistic,B(i,:));
        statvalues(i) = stats.signedrank;
    end
    s = plot_group(Group1,Group2,B1,B2,C,p,s,biom,[],unit);
    catch
        warning('The two groups do not have the same number of subjects')
    end
elseif strcmp(char(statfunc),'ttest2')
    for i = 1:nchans_o_nregs
        [h,p(i),C(i,:),stats] = ttest2(B1(i,:),B2(i,:));
        statvalues(i) = stats.tstat;
    end  
    s = plot_group(Group1,Group2,B1,B2,C,p,s,biom,[],unit);
elseif strcmp(char(statfunc),'ranksum')
     for i = 1:nchans_o_nregs
        [p(i),h,stats] = ranksum(B1(i,:),B2(i,:));
        C(i,:)=bootci(1000,{@median_diff,B1(i,:),B2(i,:)});
        statvalues(i) = stats.ranksum;
     end
    s = plot_group(Group1,Group2,B1,B2,C,p,s,biom,[],unit);
    
elseif strcmp(char(statfunc),'nbt_perm_group_diff')
    if strcmp(char(statname),'mean')
        for i = 1:nchans_o_nregs
            [p(i)] = nbt_perm_group_diff(B1(i,:),B2(i,:),'mean',5000,0);
        end
    elseif strcmp(char(statname),'median')
        for i = 1:nchans_o_nregs
            [p(i)] = nbt_perm_group_diff(B1(i,:),B2(i,:),'median',5000,0);
        end
    end
%     s.C = C;
    s.p = p;
    s.c1 = B1;
    s.c2 = B2;
    s.meanc1 = statistic(B1,2);
    s.meanc2 = statistic(B2,2);
    s.diff_biom = B2-B1;
    s.diff_mean_or_med_biom = s.meanc2-s.meanc1;
    
elseif strcmp(char(statfunc),'nbt_perm_corr')
    for i = 1:nchans_o_nregs
        [p(i)] = nbt_perm_corr(B1(i,:),B2(i,:),[],5000,0);
    end

    for i = 1:nchans_o_nregs
        corrB1(i) = corr(B1(i,:)','rows','pairwise');
        corrB2(i) = corr(B2(i,:)','rows','pairwise');
    end
    %   s.C = C;
    s.p = p;
    s.c1 = B1;
    s.c2 = B2;
    s.meanc1 =corrB1; %correlation
    s.meanc2 =corrB2; %correlation
    s.diff_biom = B2-B1;
    s.diff_mean_or_med_biom = s.meanc2-s.meanc1;% difference of correlation coeff
end

function[d]=median_diff(M,N)
        m1=nanmedian(M);
        m2=nanmedian(N);
        d=m2-m1;
end


%-------------------------------------------------
     function s = plot_group(Group1,Group2,B1,B2,C,p,s,biom,regions,unit)
            chanloc = Group1.chansregs.chanloc;
            unit = unit;
            s.C = C;
            s.p = p;
            s.c1 = B1;
            s.c2 = B2;
            s.unit = unit;
            s.meanc1 = statistic(B1,2);
            s.meanc2 = statistic(B2,2);
            if size(B1,2) == size(B2,2)
            s.diff_biom = B2-B1;
            end
            s.diff_mean_or_med_biom = s.meanc2-s.meanc1;
            % you can comment the next line, because you can also plot same
            % data using the p-value max figure
%             nbt_plot_2conditions_topo(Group1,Group2,chanloc,s,unit,biom,[]);
     end
end