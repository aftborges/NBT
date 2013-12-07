Splitting raw data files
===

The raw data files that we just [linked to in the previous step][getting_raw]
of this tutorial are very large: about 30 GB each. It is certainly possible to
use `meegpipe` to work with such large files directly, but it is generally
a good idea to try to work with smaller chunks of your data at a time, if that
makes sense for your particular analysis. Otherwise, you may need to wait a long
for every processing stage to complete on a given file.

In this tutorial we want to extract features for each experimental condition
separately. Thus, it makes sense to split our original data files into 12
single-block files, each containing just one experimental manipulation.
`meegpipe` allows you to process files into parallel jobs and thus breaking
your files into 12 smaller chunks has the potential of reducing computation
times considerably.

In the code snippets shown below we assume that all meegpipe's processing nodes
have been imported using:

````matlab
import meegpipe.node.*;
````



[getting_raw]: ./getting_raw_data.md


## Keeping your scripts organized

We are going to wrap all the scripts necessary to perform the file splitting
into a MATLAB package called `batman`. Open MATLAB and type:

````matlab
cd /data1/projects/meegpipe/batman_tut/gherrero
mkdir +batman
````

From now on we will save all scripts under `+batman`.


## Main processing script

Before writing our data processing pipeline we are going to write the scheleton
of our _main_ processing script where we perform the necessary preliminaries,
and where we run the pipeline (which we will write later) on the relevant data
files. Below you can see a profusely commented example of how such a
[split_files.m][split_files_m] script may look like:

[split_files_m]: ./split_files.m

````matlab
% SPLIT_FILES - Split BATMAN's large .mff files into single-block files
%
% This is the first stage of the BATMAN processing chain. The input to this
% stage are the raw .mff files. The produced output is a set of
% single-block .pset/pseth files (meegpipe's own data format). By
% single-block we mean a single condition block (Baseline, PVT, RS, RSQ)
% within a given experimental manipulation.


% Start in a completely clean state
close all;
clear all;
clear classes;
restoredefaultpath;

% Add meegpipe to your path, and initialize it
addpath(genpath('/data1/toolbox/meegpipe_v0.0.8'));
meegpipe.initialize;

% The output directory where we want to store the splitted data files
OUTPUT_DIR = '/data1/projects/meegpipe/batman_tut/gherrero/split_files_output';

% Some (optional) parameters that you may want to play with when experimenting
% with your processing pipeline
PARALELLIZE = true; % Should each file be processed in parallel?
DO_REPORT   = true; % Should full HTML reports be generated?

% Create an instance of your data splitting pipeline
myPipe = batman.split_files_pipeline(...
    'GenerateReport', DO_REPORT, ...
    'Parallelize',    PARALELLIZE);

% Note that we have not yet written function splitting_pipeline!

% Generate links to the relevant data files into the output directory. This
% step is equivalent to copying the relevant data files into the output
% directory but has the advantage of saving valuable disk space. The
% command below will only work at somerengrid.
files = somsds.link2rec('batman', 'subject', [1 2], 'folder', OUTPUT_DIR);

% files should now be a cell array containing the full paths to the files
% that are to be splitted (or, rather, the full paths to the symbolic links
% that point to those files).

% This is kind of obvious...
run(myPipe, files{:});
````


## The splitting pipeline


Our `split_files.m` above made used of certain `batman.split_files_pipeline`
that took care of creating an instance of the data processing pipeline. It is
now time to write that script.


### Node 1: Importing the `.mff` files

Obviously, the first step in our processing pipeline needs to be converting the
raw `.mff` files into [physioset][physioset] objects, which is the data
structure that _meegpipe_'s processing nodes understand.

[physioset]: ../../+physioset/@physioset/README.md

Importing data from various disk file formats into a _physioset_ object is
always performed with a [physioset_import][physioset_import_node] node. Below
you have a schematic diagram of such a node:

[physioset_import_node]: ../../+meegpipe/+node/+physioset_import/README.md

![physioset_import node](./img/physioset_import_node.png "physioset_import node")

__NOTE:__ In the diagram above I have depicted also the _data selection_ steps
that take place before and after the actual node processing. These steps are
common to all classes of processing nodes and, in the following, they will not
be shown in the node diagrams. See the documentation of the generic [node][node]
interface for more information.

[node]: ../../+meegpipe/+node/README.md

Nodes of class `physioset_import` admit only one configuration option,
_Importer_. The user needs to set it to an object of one of the physioset
importer classes that are available in package [physioset.import][physioset_import_pkg].
For instance, you may use a `physioset.import.edfplus` object for
importing [EDF+][edfplus_format] files. In our case, the raw data files are in
`.mff` format and thus we need a `physioset.import.mff` importer.

Realize that, although the `physioset_import` node has just one configuration
option (_Importer_), the actual data importer object has several properties that
allow you to specify various importing options. Let's build a _default_ mff
importer object to find out what properties it has:

[physioset_import_pkg]: ../../+physioset/+import/README.md
[edfplus_format]: http://www.edfplus.info/

````matlab
>> physioset.import.mff

ans =

mff
Package: physioset.import


      ReadDataValues : true
           Precision : double
            Writable : true
           Temporary : true
           ChunkSize : 500000000
   AutoDestroyMemMap : false
          ReadEvents : true
            FileName :
          FileNaming : inherit
             Sensors : []
        EventMapping : [1x1 mjava.hash]
          MetaMapper : @(data)regexp(get_name(data),'(?<subject_id>0+\d+)_.+','names')
         EventMapper : []
           StartTime :
````

One property that you may often want to override is the `Precision` property,
which determines the numeric precision that is used to store the values
contained in the generated `physioset` object. The code below will create a
`physioset_importer` node that will convert an `.mff` data file into
a `physioset` object of `single` precision:

````matlab
import meegpipe.node.*; % Just to avoid writing fully qualified node class names

% Let's build our data importer with a custom Precision value
myImporter = physioset.import.mff('Precision', 'single');

% Now, let's build a physioset_import node that uses the importer object above
myNode = physioset_import.new('Importer', myImporter);
````

And that's it for the first node in our file splitting pipeline. Let's move on
to the next node.


### Node 2: Splitting the imported physioset

Nodes of class [split][split_node] allow you to split a physioset object into
several (possible overlapping) data subsets. Below you have the node schema:

[split_node]: ../../+meegpipe/+node/+split/README.md


![split node](./img/split_node.png "split node")


Nodes of class `split` determine the time range of each data split by selecting
a subset of events among the events present in the input physioset. The subset
of relevant events is selected using a user-defined `EventSelector`. To
understand how this works, let's take a look at the default value of
`EventSelector`:

````matlab
defEvSel = get_config(split.new, 'EventSelector')

ans =

class_selector
Package: physioset.event


               Class : { 'split_begin'}
                Type : {}
             Negated : false
                Name :
````

So, by default, the `split` node looks for `split_begin` events in the input
physioset and splits the dataset according to the `Offset` and `Duration` of
such events. The example below illustrates this behavior with a random toy
physioset:

````matlab
% Create a toy physioset object that contains random values
myImporter = physioset.import.matrix('FileName', 'mytoydata');
myPhysioset = import(myImporter, rand(3,1000));

% In your current directory there should now be a file called `mytoydata.pset`
% which contains the physioset data values.

% Add two split_begin events at samples 100 and 500, with a duration of 500 samples
myEvArray = physioset.event.std.split_begin([100 500], 'Duration', 500);
add_event(myPhysioset, myEvArray);

% Split it!
myNode = split.new;
run(myNode, myPhysioset);
````

The messages displayed in the MATLAB command window should provide you all the
details regarding the two data splits that were just generated. The same
information together with instructions on how to retrieve the produced data
splits is found in the HTML report generated by the `split` node. Open the
report and click on the _Split mytoydata_split1_ link, under the _Produced data
splits_ section. At the bottom of the page you will find the code snippet that
you need to retrieve the first split physioset. In my case the code snippet
reads:

````matlab
data = pset.load(['/Volumes/DATA/work/' ...
    'mytoydata.meegpipe/split-3fa41e_gherrero_MACI64-R2012a/' ...
    'mytoydata_split1.pseth'])
````

Let's check whether the first split contains the expected values:

````matlab
% The first split should contain these data values:
myManualSplit1 = myPhysioset(:, 100:599);

% This should not throw an error if everything went OK
assert(all(data(:)==myManualSplit1(:)));
````




#### Splitting strategy


![splitting strategy](./img/batman_protocol_subblock.png "splitting strategy")