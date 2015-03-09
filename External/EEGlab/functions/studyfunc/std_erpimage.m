<<<<<<< HEAD
% std_erpimage() - Compute ERP images and save them on disk.
%
% Usage:
%   >> std_erpimage( EEG, 'key', 'val', ...);
%
% Inputs:
%   EEG          - a loaded epoched EEG dataset structure. May be an array
%                  of such structure containing several datasets.
%
% Optional inputs:
%   'components' - [numeric vector] components of the EEG structure for which 
%                  the measure will be computed {default|[] -> all}
%   'channels'   - [cell array] channels of the EEG structure for which 
%                  activation ERPs will be computed {default|[] -> none}
%   'trialindices' - [cell array] indices of trials for each dataset.
%                  Default is all trials.
%   'recompute'  - ['on'|'off'] force recomputing data file even if it is 
%                  already on disk.
%   'rmcomps'    - [integer array] remove artifactual components (this entry
%                  is ignored when plotting components). This entry contains 
%                  the indices of the components to be removed. Default is none.
%   'interp'     - [struct] channel location structure containing electrode
%                  to interpolate ((this entry is ignored when plotting 
%                  components). Default is no interpolation.
%   'fileout'    - [string] name of the file to save on disk. The default
%                  is the same name (with a different extension) as the 
%                  dataset given as input.
%
% ERPimage options:
%   'concatenate' - ['on'|'off'] concatenate single trial of different
%                  subjects for plotting ERPimages ('on'). The default
%                  ('off') computes an ERPimage for each subject and then
%                  averages these ERPimages. This allows to perform
%                  statistics (the 'on' options does not allow statistics).
%   'smoothing'  - Smoothing parameter (number of trials). {Default: 10}
%                  erpimage() equivalent: 'avewidth'
%   'nlines'     - Number of lines for ERPimage. erpaimge() equivalent is 
%                  'decimate'. Note that this parameter must be larger than
%                  the minimum number of trials in each design cell 
%                  {Default: 10}
%   'sorttype'   - Sorting event type(s) ([int vector]; []=all). See Notes below.
%                  Either a string or an integer.
%   'sortwin'    - Sorting event window [start, end] in milliseconds ([]=whole epoch)
%   'sortfield'  - Sorting field name. {default: latency}.
%   'erpimageopt'  - erpimage() options, separated by commas (Ex: 'erp', 'cbar').
%                  {Default: none}. For further details see >> erpimage help
% Outputs:
%   erpimagestruct - structure containing ERPimage information that is
%                    been saved on disk.
%
%   Files are saved on disk.
%    [dataset_file].icaerpim     % component ERPimage file
% OR
%    [dataset_file].daterpim     % channel ERPimage file
%
% Author: Arnaud Delorme, SCCN & CERCO, CNRS, 2011-

% Copyright (C) 2011 Arnaud Delorme
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

function allerpimage = std_erpimage( EEG, varargin);

if nargin < 1
    help std_erpimage;
    return;
end

allerpimage = [];
[opt moreopts] = finputcheck( varargin, { ...
    'components'    'integer'     []                    [];
    'channels'      { 'cell','integer' }  { [] [] }     {};
    'trialindices' { 'integer','cell' }   []            [];
    'recompute'     'string'      { 'on','off' }        'off';
    'savefile'      'string'      { 'on','off' }        'on';
    'fileout'       'string'      []                    '';
    'rmcomps'       'cell'        []                    cell(1,length(EEG));
    'interp'        'struct'      { }                   struct([]);
    'nlines'        ''            []                    10;
    'smoothing'     ''            []                    10;
    'sorttype'       ''           {}                    '';
    'sortwin'        ''           {}                    [];
    'sortfield'      ''           {}                    'latency';
    'concatenate'   'string'      { 'on','off' }        'off';
    'erpimageopt'   'cell'        {}                    {}}, ...
    'std_erpimage', 'ignore');
if isstr(opt), error(opt); end;
if length(EEG) == 1 && isempty(opt.trialindices), opt.trialindices = { [1:EEG.trials] }; end;
if isempty(opt.trialindices), opt.trialindices = cell(length(EEG)); end;
if ~iscell(opt.trialindices), opt.trialindices = { opt.trialindices }; end;
if isfield(EEG,'icaweights')
    numc = size(EEG(1).icaweights,1);
else
    error('EEG.icaweights not found');
end
if isempty(opt.components)
    opt.components = 1:numc;
end

% filename
% --------
if isempty(opt.fileout), opt.fileout = fullfile(EEG(1).filepath, EEG(1).filename(1:end-4)); end;
if ~isempty(opt.channels)
    filenameshort = [ opt.fileout '.daterpim'];
    prefix = 'chan';
    if iscell(opt.channels)
        if ~isempty(opt.interp)
            opt.indices = eeg_chaninds(opt.interp, opt.channels, 0);
        else
            opt.indices = eeg_chaninds(EEG(1), opt.channels, 0);
            for ind = 2:length(EEG)
                if ~isequal(eeg_chaninds(EEG(ind), opt.channels, 0), opt.indices)
                    error([ 'Channel information must be consistant when ' 10 'several datasets are merged for a specific design' ]);
                end;
            end;
        end;
    else
        opt.indices = opt.channels;
    end;
else
    opt.indices = opt.components;
    filenameshort = [ opt.fileout '.icaerpim'];
    prefix = 'comp';
end;
filename = filenameshort;

% ERP information found in datasets
% ---------------------------------
if exist(filename) && strcmpi(opt.recompute, 'off')
    fprintf('File "%s" found on disk, no need to recompute\n', filenameshort);
    return;
end

allerpimage = [];
if strcmpi(opt.concatenate, 'off')
    % compute ERP images
    % ------------------
    if isempty(opt.channels)
         X = eeg_getdatact(EEG, 'component', opt.indices, 'trialindices', opt.trialindices );
    else X = eeg_getdatact(EEG, 'channel'  , opt.indices, 'trialindices', opt.trialindices, 'rmcomps', opt.rmcomps, 'interp', opt.interp);
    end;
    if ~isempty(opt.sorttype)
         events = eeg_getepochevent(EEG, 'type', opt.sorttype, 'timewin', opt.sortwin, 'fieldname', opt.sortfield, 'trials', opt.trialindices);
    else events = [];
    end;
        
    % reverse engeeneering the number of lines for ERPimage
    finallines = opt.nlines;
    if ~isempty(events)
         if all(isnan(events))
             error('Cannot sort trials for one of the dataset');
         end;
         lastx  = sum(~isnan(events));
    else lastx  = size(X,3);
    end;
    if lastx < finallines + floor((opt.smoothing-1)/2) + 3
        error('The default number of ERPimage lines is too large for one of the dataset');
    end;
    firstx = 1;
    xwidth = opt.smoothing;
    %xadv   = lastx/finallines;
    nout   = finallines; %floor(((lastx-firstx+xadv+1)-xwidth)/xadv);
    nlines = (lastx-xwidth)/(nout-0.5)*i; % make it imaginary
    %nlines = ceil(lastx/((lastx-firstx+1-xwidth)/(nout-1)));

    if 0
        % testing conversion back and forth
        % ---------------------------------
        for lastx = 20:300
            for xwidth = 1:19
                for nlines = (xwidth+1):100
                    
                    nout = floor(((lastx+nlines)-xwidth)/nlines);
                    realnlines = (lastx-xwidth)/(nout-0.5);
                    noutreal = floor(((lastx+realnlines)-xwidth)/realnlines);
                    
                    if nout ~= noutreal
                        error('Wrong conversion 2');
                    end;
                    
                end;
            end;
        end;
    end;
    
    clear tmperpimage eventvals;
    parfor index = 1:size(X,1)
        [tmpX tmpevents] = erpimage(squeeze(X(index,:,:)), events, EEG(1).times, '', opt.smoothing, nlines, 'noplot', 'on', opt.erpimageopt{:}, moreopts{:});
        if isempty(events), tmpevents = []; end;
        eventvals{index}   = tmpevents;
        tmperpimage{index} = tmpX';
    end;
    allerpimage.events = eventvals{1};
    for index = 1:size(X,1)
        allerpimage.([ prefix int2str(opt.indices(index)) ]) = tmperpimage{index};
    end;
else
    % generate dynamic loading commands
    % ---------------------------------
    for dat = 1:length(EEG)
        filenames{dat} = fullfile(EEG(1).filepath, EEG(1).filename);
    end;
    allerpimage.times = EEG(1).times;
    for index = 1:length(opt.indices)
        if ~isempty(opt.channels)
             com = sprintf('squeeze(eeg_getdatact(%s, ''interp'', chanlocsforinterp));', vararg2str( { filenames 'channel'  , opt.indices(index), 'rmcomps', opt.rmcomps, 'trialindices', opt.trialindices } ));
        else com = sprintf('squeeze(eeg_getdatact(%s));', vararg2str( { filenames 'component', opt.indices(index), 'trialindices', opt.trialindices } ));
        end;
        allerpimage = setfield(allerpimage, [ prefix int2str(opt.indices(index)) ], com);
    end;
    allerpimage = setfield(allerpimage, 'chanlocsforinterp', opt.interp);
    if ~isempty(opt.sorttype)
         events = eeg_getepochevent(EEG, 'type', opt.sorttype, 'timewin', opt.sortwin, 'fieldname', opt.sortfield, 'trials', opt.trialindices);
         %geteventcom = sprintf('eeg_getepochevent(%s);', vararg2str( { filenames 'type', opt.sorttype, 'timewin', opt.sortwin, 'fieldname', opt.sortfield } ));
    else events = [];
    end;
    allerpimage = setfield(allerpimage, 'events', events);
end;
allerpimage.times       = EEG(1).times;
allerpimage.parameters  = varargin;
allerpimage.datatype    = 'ERPIMAGE';
allerpimage.datafiles   = computeFullFileName( { EEG.filepath }, { EEG.filename });
allerpimage.datatrials  = opt.trialindices;

% Save ERPimages in file (all components or channels)
% ----------------------------------------------
if strcmpi(opt.savefile, 'on')
    if strcmpi(prefix, 'comp')
        std_savedat(filename, allerpimage);
    else
        tmpchanlocs = EEG(1).chanlocs;
        allerpimage.labels = opt.channels;
        std_savedat(filename, allerpimage);
    end;
end;

% compute full file names
% -----------------------
function res = computeFullFileName(filePaths, fileNames);
for index = 1:length(fileNames)
    res{index} = fullfile(filePaths{index}, fileNames{index});
end;
=======
% std_erpimage() - This is a legacy function. This function is not
%                  compatible yet with STUDY design. 
%
%                  Plot an erpimage using multiple subject data
%
% Usage:
%   >> std_erpimage( STUDY, ALLEEG, 'key', 'val', ...);
% Inputs:
%   STUDY    - STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG   - global vector of EEG structures for the datasets included 
%              in the STUDY. ALLEEG for a STUDY set is typically created 
%              using load_ALLEEG().  
% Optional inputs:
%   'clusters' - [numeric vector|'all'] indices of clusters to plot.
%                If no component indices ('comps' below) are given, the average 
%                ERSPs of the requested clusters are plotted in the same figure, 
%                with ERSPs for different conditions (and groups if any) plotted 
%                in different colors. In 'comps' (below) mode, ERSP for each 
%                specified cluster are plotted in separate figures (one per 
%                condition), each overplotting cluster component ERSP plus the
%                average cluster ERSP in bold. Note this parameter has no effect 
%                if the 'comps' option (below) is used. {default: 'all'}
%   'comps'    - [numeric vector|'all'] indices of the cluster components to plot.
%                Note that 'comps', 'all' is equivalent to 'plotsubjects', 'on'.
%
% Optional inputs for channel plotting:
%   'channels' - [numeric vector]  specific channel group to plot. By
%                default, the grand mean channel ERSP is plotted (using the 
%                same format as for the cluster component means described above)
%   'subject'  - [numeric vector]  In 'changrp' mode (above), index of 
%                the subject(s) to plot. Else by default, plot all components 
%                in the cluster.
%   'plotsubjects' - ['on'|'off'] When 'on', plot ERSP of all subjects.
%
% Commandline options:
%   'projchan' - Channel to back-project the selected component(s) to. 
%                If plotting channel activity, this argument is ignored. 
%                If [], the ICA component activation is plotted {default []}.
%   'title'      - ['string'] Plot title {default: []}
%   'smooth'     - Smoothing parameter (number of trials). {Default: 5} 
%                erpimage() equivalent: 'avewidth' 
%   'decimate'   - Decimate parameter (i.e. ratio of trials_in/trials_out).
%                erpaimge() equivalent: 'decimate' {Default: 0}
%   'sorttype'   - Sorting event type(s) ([int vector]; []=all). See Notes below.
%                Either a string or an integer.
%   'sortwin'    - Sorting event window [start, end] in seconds ([]=whole epoch)
%   'sortfield'  - Sorting field name. {default: none}. See Notes below.
%   'erpimageopt'  - erpimage() options, separated by commas (Ex: 'erp', 'cbar'). 
%                {Default: none}. For further details see >> erpimage help   
%
% Outputs:
%  STUDY       - STUDY structure 
%  allphases   - all phase from erpimage
%  allerpimage - all sorting variables from erpimage
%  subjamp     - value of amplitude contribution [0 to 1] for each trial.
%                1 indicates that all subject contribute equally to the
%                amplitude of each trial. 0 indicates that only one
%                subject dominates all trials' amplitude.
%
% Author: Arnaud Delorme, SCCN & CERCO, CNRS, 2004-

% Copyright (C) 2004 Arnaud Delorme
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

function [STUDY, allphases, allsortvar, subjamptime, subjamptrial, globalent ] = std_erpimage( STUDY, ALLEEG, varargin);
        
    error('std_erpimage is a legacy function. This function is not compatible yet with STUDY design scheme');

    if nargin < 3
        help std_erpimage;
        return;
    end;
    
    STUDY = pop_erpimparams(STUDY, 'default');

    [ opt moreparams ] = finputcheck( varargin, { ...
        'erpimageopt' 'cell'    [] STUDY.etc.erpimparams.erpimageopt;
        'sorttype'    { 'string' 'cell' } [] STUDY.etc.erpimparams.sorttype;
        'sortwin'     'real'    [] STUDY.etc.erpimparams.sortwin;
        'sortfield'   'string'  [] STUDY.etc.erpimparams.sortfield;
        'statistics'  'string'  [] STUDY.etc.erpimparams.statistics;
        'groupstats'  'string'  [] STUDY.etc.erpimparams.groupstats;
        'condstats'   'string'  [] STUDY.etc.erpimparams.condstats;
        'statmode'    'string'  [] STUDY.etc.erpimparams.statmode;
        'threshold'   'real'    [] STUDY.etc.erpimparams.threshold;
        'naccu'       'integer' [] STUDY.etc.erpimparams.naccu;
        'channels'    'cell'    []              {};
        'clusters'    'integer' []              [];
        'mode'        'string'  []              '';
        'comps'       {'integer','string'}  []              []; % for backward compatibility
        'plotsubjects' 'string' { 'on' 'off' }  'off';
        'normalize'    'string' { 'on' 'off' }  'off';
        'subject'     'string'  []              '' }, ...
                                      'std_erpimage', 'ignore');
    if isstr(opt), error(opt); end;
    
    if strcmpi(opt.mode, 'comps'), opt.plotsubjects = 'on'; end;

    if ~isempty(opt.subject), opt.groupstats = 'off'; disp('No group statistics for single subject'); end;
    if ~isempty(opt.subject), opt.condstats = 'off'; disp('No condition statistics for single subject'); end;
    
    if length(opt.comps) == 1
        opt.condstats = 'off'; opt.groupstats = 'off'; 
        disp('Statistics cannot be computed for single component');
    end;

    % read data from disk
    % -------------------
    if ~isempty(opt.channels)
         [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'channels', opt.channels, 'infotype', 'data');
         [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'channels', opt.channels, 'infotype', 'event', ...
                                            'type', opt.sorttype, 'timewin', opt.sortwin, 'fieldname', opt.sortfield);
    else [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'clusters', opt.clusters, 'infotype', 'data');
         [STUDY tmp allinds] = std_readdata(STUDY, ALLEEG, 'clusters', opt.clusters, 'infotype', 'event', ...
                                            'type', opt.sorttype, 'timewin', opt.sortwin, 'fieldname', opt.sortfield);
    end;
    opt.legend = 'off';
    
    for index = 1:length(allinds)

        if length(allinds) > 1, try, subplot(nr,nc,index, 'align'); catch, subplot(nr,nc,index); end; end;
        if ~isempty(opt.channels)
            eval( [ 'alldata     = STUDY.changrp(allinds(index)).data;' ]);
            eval( [ 'alltimes    = STUDY.changrp(allinds(index)).datatimes;' ]);
            eval( [ 'allvals     = STUDY.changrp(allinds(index)).datasortvals;' ]);
            eval( [ 'allcontinds = STUDY.changrp(allinds(index)).datacontinds;' ]);
            setinds  = STUDY.changrp(allinds(index)).setinds;
        else
            eval( [ 'alldata  = STUDY.cluster(allinds(index)).data;' ]);
            eval( [ 'alltimes    = STUDY.cluster(allinds(index)).datatimes;' ]);
            eval( [ 'allvals     = STUDY.cluster(allinds(index)).datasortvals;' ]);
            eval( [ 'allcontinds = STUDY.cluster(allinds(index)).datacontinds;' ]);
            compinds = STUDY.cluster(allinds(index)).allinds;
            setinds  = STUDY.cluster(allinds(index)).setinds;
        end;
        
        % plot specific subject
        % ---------------------
        if ~isempty(opt.subject) & isempty(opt.comps)
            for c = 1:size(alldata,1)
                for g = 1:size(alldata,2)
                    for l=length(setinds{c,g}):-1:1
                        if ~strcmpi(opt.subject, STUDY.datasetinfo(setinds{c,g}(l)).subject)
                            inds = find(allcontinds == l);
                            alldata{c,g}(:,inds) = [];
                            allvals{c,g}(inds) = [];
                            allcontinds{c,g}(inds) = [];
                        end;
                    end;
                end;
            end;
        end;
        
        % plot specific component
        % -----------------------
        if ~isempty(opt.comps)
            
            % find and select group
            % ---------------------
            sets   = STUDY.cluster(allinds(index)).sets(:,opt.comps);
            comps  = STUDY.cluster(allinds(index)).comps(opt.comps);
            grp    = STUDY.datasetinfo(sets(1)).group;
            grpind = strmatch( grp, STUDY.group );
            if isempty(grpind), grpind = 1; end;
            alldata = alldata(:,grpind);
            
            % find component
            % --------------
            for c = 1:size(alldata,1)
                for ind = length(compinds{1,grpind}):-1:1
                    if ~any(compinds{1,grpind}(ind) == comps) | ~any(setinds{1,grpind}(ind) == sets)
                        inds = find(allcontinds == ind);
                        alldata{c,g}(:,inds) = [];
                        allvals{c,g}(inds) = [];
                        allcontinds{c,g}(inds) = [];
                    else
                        comp_names{c,1} = comps;
                    end;
                end;
            end;
            opt.subject = STUDY.datasetinfo(sets(1)).subject;
        end;
        
        % normalize data
        % --------------
        if strcmpi(opt.normalize, 'on')
            indbef0 = find(alltimes < 0);
            for c = 1:size(alldata,1)
                for g = 1:size(alldata,2)
                    for tmpind = 1:length(STUDY.subject)
                        inds = find(allcontinds{c,g} == tmpind);
                        alldata{c,g}(:,inds) = alldata{c,g}(:,inds) - ...
                                                repmat(mean(alldata{c,g}(indbef0,inds),1), [size(alldata{c,g},1) 1]);
                        alldata{c,g}(:,inds) = alldata{c,g}(:,inds) / std(reshape(alldata{c,g}(indbef0,inds), ...
                                                                          [length(inds)*length(indbef0) 1] ));
                    end;
                end;
            end;
        end;
        
        % plot specific component
        % -----------------------
        if index == length(allinds), opt.legend = 'on'; end;
        figure('color', 'w');
        for c = 1:size(alldata,1)
            for g = 1:size(alldata,2)
                if ~isempty(alldata{c,g})
                    subplot(size(alldata,1), size(alldata,2),c+(g-1)*size(alldata,1));
                    if ~isempty(alldata{c,g})
                        if all(isnan(allvals{c,g}))
                            erpimage(alldata{c,g}, [], alltimes, opt.erpimageopt{:});
                        else
                            erpimage(alldata{c,g}, allvals{c,g}, alltimes, opt.erpimageopt{:});
                        end;
                    end;
                end;
            end;
        end;

        if isempty(opt.channels), %title(sprintf('Cluster %d', allinds(index)));
            title([ STUDY.cluster(allinds(index)).name ' (' num2str(length(STUDY.cluster(allinds(index)).comps)),' ICs, ' ...
                    num2str(length(unique(STUDY.cluster(allinds(index)).sets(1,:)))) ' Ss)' ]);            
            set(gcf,'name','Cluster ERPIMAGE')
        else                      
            title(sprintf('%s', opt.channels{index}));  
            set(gcf,'name','Channel grand ERPIMAGE')
        end;
    end;
    
>>>>>>> FETCH_HEAD
