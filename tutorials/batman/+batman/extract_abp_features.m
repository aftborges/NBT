function extract_abp_features
% EXTRACT_ABP_FEATURES - Extract ABP features from BATMAN data

% Start in a completely clean state
close all;
clear all;
clear classes;

% Add meegpipe to your path, and initialize it
switch lower(get_hostname),
    case {'somerenserver', 'nin389'},
        addpath(genpath('/data1/toolbox/meegpipe_v0.0.8'));
    case 'outolintulan',
        addpath(genpath('/Volumes/DATA/mlib/meegpipe_v0.0.8'));
    otherwise
        error('I don''t know where is meegpipe on %s', get_hostname);
end

meegpipe.initialize;

% Import some utilities
import mperl.file.find.finddepth_regex_match;

switch lower(get_hostname),
    case {'somerenserver', 'nin389'},
        % The directory where the split data files are located
        INPUT_DIR = ...
            '/data1/projects/meegpipe/batman_tut/gherrero/split_files_output';
        % The output directory where we want to store the features
        OUTPUT_DIR = ...
            '/data1/projects/meegpipe/batman_tut/gherrero/extract_abp_features_output';        
    case 'outolintulan'
        INPUT_DIR = '/Volumes/DATA/tutorial/batman/split_files_output';
        OUTPUT_DIR = '/Volumes/DATA/tutorial/batman/extract_abp_features_output';
end

% Some (optional) parameters that you may want to play with when experimenting
% with your processing pipeline
PARALELLIZE = true; % Should each file be processed in parallel?
DO_REPORT   = true; % Should full HTML reports be generated?

% Create an instance of the feature extraction pipeline
myPipe = batman.extract_abp_features_pipeline(...
    'GenerateReport', DO_REPORT, ...
    'Parallelize',    PARALELLIZE);

% Note that we have not yet written function extract_abp_feature_pipeline!

% Generate links to the relevant data files into the output directory. This
% step is equivalent to copying the relevant data files into the output
% directory but has the advantage of saving valuable disk space.
regex = 'split_files-.+_\d+\.pseth?';
splittedFiles = finddepth_regex_match(INPUT_DIR, regex, false);
somsds.link2files(splittedFiles, OUTPUT_DIR);
regex = '\.pseth$';
files = finddepth_regex_match(OUTPUT_DIR, regex);

% files should now be a cell array containing the full paths to the single
% sub-block .pseth files that were generated in the data splitting stage.

run(myPipe, files{:});

end

function name = get_hostname
[ret, name] = system('hostname');
if ret
    if ispc,
        name = getenv('COMPUTERNAME');
    else
        name = getenv('HOSTNAME');
    end
end
name = regexprep(name, '[^\w\d]', '');
end