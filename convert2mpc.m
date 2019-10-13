function [mpc, contingencies] = convert2mpc(...
                   rawfile_name, ropfile_name, inlfile_name, confile_name)
%CONVERT2MPC Converts the GOCompetition dataset to a MatPower case
%	MPC, CONTINGENCIES = CONVERT2MPC()
%	MPC, CONTINGENCIES = CONVERT2MPC(RAWFILE_NAME, ROPFILE_NAME,
%                     INLFILE_NAME,CONFILE_NAME)
%	Converts the GOCompetition dataset to a MatPower case. This function
%	needs the CASE.RAW, CASE.ROP, CASE.INL and CASE.CON files (PSSE 
%   rev 33.7 raw file), as downloaded from the GOCompetition website.
%
%	If the directory where the downloaded dataset files is not the current
%	directory nor in the Matlab Path, or the filenames have been changed
%	from their default values, specify location and filename
%	in the arguments.
%	Else, you may use the function without arguments.
%
%
%   Current version: v11
%
%   About this version:
%   Added field MPC.ISXFMR to tell apart xfmrs and non-xfmr branches.
%   Added field CONTINGENCIES.OFFLINE to save the data about elements
%   offline that appear in the list of CONS. These elements are no longer
%   included in CONTINGENCIES.BRANCH nor CONTINGENCIES.GEN
%
%   About version 10
%   Fixed bug detected when parsing some INL files.
%   Fixed generators are now converted to (possibly negative) loads. Fixed
%   generators are generators for which PMIN=PMAX and QMIN=QMAX. The
%   information about these generators is located in field fixedGen.
%   
%   About version 9:
%   Some datasets were issued without the 'END OF ... DATA, BEGIN OF
%   GENERATOR DATA', which was used to define section limits. This function
%   was adapted accordingly.
%   About version 8:
%   Added contingency limits for voltages
%
%   Adapted to new problem formulation and datasets, using RAW, CON, INL
%   and ROP files.
%
%   Added FIXEDGEN field to MPC struct. FIXEDGEN is an NG-by-1 vector,
%   where NG is the number of generators in the system (counting offline
%   generators). The i-th element in the vector is 1 if the generator was
%   converted to a load, 0 otherwise.
%   
%   Added INDEXMAP field to MPC struct.
%   indexMap
%           .branch
%                  .psse2mpc
%                  .mpc2psse
%           .gen
%                  .psse2mpc
%                  .mpc2psse
%
%   These are maps from PSSE generator and branch keys to
%   matpower case indices and vice-versa.
%   The maps are of the Map class (similar to python dictionaries)
%   The map's keys are PSS/E's unique generator keys, built as
%   the concatenation of bus number and generator ID
%   e.g. '1-1' for a generator;
%   and branch unique keys, built as the concatenation of
%   busFromNumber, busToNumber and circuit
%   e.g. '1-2-BL' for a branch
%   The map's value is the generator or branch index in the 
%   corresponding MPC matrix
%   The map in the opposite direction is identical but with key and values
%   inverted.
%
%   Added field COMP, a matrix where each row represents a bus in which
%   there is a fixed or switched reactive shunt, and with the following
%   columns:
%   1: SW_BUS: Bus number
%   2: FS: Fixed shunt (MVAR @ 1.0 p.u)
%   3: SSMAX: Switched shunt max (MVAR @ 1.0 p.u)
%   4: SSMIN: Switched shunt min (MVAR @ 1.0 p.u)
%   5: SW_STATUS: Switched shunt status (1/0)
%
%   Added field ISXFMR
%   ISXFMR is a NB-by-1 column boolean vector, where NB is the number of
%   branches. The i-th element is true if branch corresponding to row i of
%   branch matrix is a xfmr, false otherwise.
%
%
%   Copyright (c) 2018, Gers USA
%   by Tomas Valencia tvalencia@gersusa.com

%% Define filenames default values
if nargin<4
    confile_name = 'case.con';
    if nargin<3
        inlfile_name = 'case.inl';
        if nargin<2
            ropfile_name = 'case.rop';
            if nargin==0
                rawfile_name = 'case.raw';
            end
        end
    end
end

%% Define constants
define_constants();

%% Main function

% Get system topology from PSSE raw file
mpc = psse2mpc_mod(rawfile_name,0,33.7);

% Get generator droop participation factors
genDroops = readParseINLFile(inlfile_name);

% Get generator cost data
[genCosts, maxNpairs] = readParseRopFile(ropfile_name);

% Write read droop and cost data into MatPower case
indexMap = buildIndexMap(rawfile_name);
mpc.indexMap = indexMap;
genIndexMap = indexMap.gen.psse2mpc;
genKeys = keys(genDroops);
gencost = zeros(length(genKeys),2*maxNpairs+4);
for i=1:length(genKeys)
    genIndex = genIndexMap(genKeys{i});
    mpc.gen(genIndex,21) = genDroops(genKeys{i});
    gencost(genIndex,1) = 1;
    genCostPCW = genCosts(genKeys{i});
    [nPairs,genCostPCW_corr] = getNPairs(genCostPCW);
    gencost(genIndex,4) = nPairs;
    gencost(genIndex,5:4+2*nPairs) = genCostPCW_corr;
end
mpc.gencost = gencost;

% Convert fixed generators to loads
fG = mpc.gen(:,GEN_STATUS)~=0 & mpc.gen(:,PMIN)==mpc.gen(:,PMAX) & ...
    mpc.gen(:,QMIN)==mpc.gen(:,QMAX);
% Using a for instead of vector assigning because there is no guarantee
% that values are ordered
fGI = find(fG);
for b=1:sum(fG)
    fGB = mpc.bus(:,BUS_I)==mpc.gen(fGI(b),GEN_BUS);
    mpc.bus(fGB,PD) = mpc.bus(fGB,PD) - mpc.gen(fGI(b),PMAX);
    mpc.bus(fGB,QD) = mpc.bus(fGB,QD) - mpc.gen(fGI(b),QMAX);
end
% Turn off fixed gens and save data
mpc.gen(fG,GEN_STATUS)=0;
mpc.fixedGen = fG;

% Get contingency data
[branOutMap, genOutMap] = readParseCONFile(confile_name);

% Check if any of them are offline elements, save them and delete them
offline.branch = [];
offline.gen = [];
if ~isempty(branOutMap)
    cBrI = cell2mat(values(indexMap.branch.psse2mpc,values(branOutMap)));
    offB = mpc.branch(cBrI,BR_STATUS)==0;
    brk = keys(branOutMap);
    brv = values(branOutMap);
    if any(offB)
        offline.branch = containers.Map(brk(offB),brv(offB));
        remove(branOutMap,brk(offB));
    end
end
if ~isempty(genOutMap)
    cGI = cell2mat(values(indexMap.gen.psse2mpc,values(genOutMap)));
    offG = mpc.gen(cGI,GEN_STATUS)==0;
    gk = keys(genOutMap);
    gv = values(genOutMap);
    if any(offG)
        offline.gen = containers.Map(gk(offG),gv(offG));
        remove(genOutMap,gk(offG));
    end
end
if numel(unique(mpc.bus(:,BUS_AREA)))>1
    onGen(:,1)=mpc.gen(cGI,GEN_BUS);
    for i=1:size(onGen,1)
        genAreas(i,1)=onGen(i,1);
        genAreas(i,2)=mpc.bus(find(onGen(i,1)==mpc.bus(:,BUS_I),1),BUS_AREA);
    end
    areasGen=arrayfun(@(x) sum(genAreas(:,2)==x), unique(mpc.bus(:,BUS_AREA)));
    if any(areasGen==1)
        areaOneGen=find(areasGen==1);
        oneGen=genAreas(ismember(genAreas(:,2),areaOneGen),1);
        oneGen=find(ismember(mpc.gen(:,GEN_BUS),oneGen));
        oneGenCI = arrayfun(@(x) find(strcmp(mpc.indexMap.gen.mpc2psse(x),values(genOutMap))),oneGen);
        offline.oneGenInArea = containers.Map(gk(oneGenCI),gv(oneGenCI));
        remove(genOutMap,gk(oneGenCI));
%     else
%         offline.oneGenInArea = [];
    end
end

contingencies.branch = branOutMap;
contingencies.gen = genOutMap;
contingencies.offline = offline;
end
%% Nested Functions
%=====================================================================
%=====================================================================
%=====================================================================

function [ indexMap ] = buildIndexMap( rawfile_name )
%BUILDINDEXMAP 
%   This function builds maps from PSSE generator and branch keys to
%   matpower case indices and vice-versa.
%   The maps are of the Map class (similar to python dictionaries)
%   The map's keys are PSS/E's unique generator keys, built as
%   the concatenation of bus number and generator ID
%   e.g. '1-1' for a generator;
%   and branch unique keys, built as the concatenation of
%   busFromNumber, busToNumber and circuit
%   e.g. '1-2-BL' for a branch
%   The map's value is the generator or branch index in the 
%   corresponding MPC matrix
%   The map in the opposite direction is identical but with key and values
%   inverted.
%   The result is returned in a structure:
%
%   indexMap
%           .branch
%                  .psse2mpc
%                  .mpc2psse
%           .gen
%                  .psse2mpc
%                  .mpc2psse
%
%
%   genIndexMap, branchIndexMap  = BUILDINDEXMAP( rawfile_name )
[~, ~, ext] = fileparts(rawfile_name);
if nargin<1
    error('Not enough arguments')
elseif ~strcmpi(ext,'.RAW')
    disp(ext)
    error('Not a valid RAW file')
end

% Read RAW file
rawFile = fileread(rawfile_name);

% Get indices of beginning of sections which are the ones starting with '0
% / [text]'
sectionIndices = regexp(rawFile,'\s+0 */ *.*?(?=\n)');



%=====================================================
% Create generator index map
%=====================================================

% Get generator information section, which is: 
ind1 = sectionIndices(3); % Generation dispatch data is section 4
ind2 = sectionIndices(4);
genSection = rawFile(ind1:ind2-1);
genSection = strsplit(strtrim(genSection),'\n');

genSection(1)=[]; % Discard 'BEGIN GENERATOR DATA' text

% Build the map  
genKeys = cell(length(genSection),1);
for i=1:length(genSection)
    genRow = genSection{i};
    items = strsplit(genRow,',');
    genKey = strcat(strtrim(items{1}),...
                    '-',strrep(strtrim(items{2}),'''',''));
	genKeys{i} = genKey;
end

gen.psse2mpc = containers.Map(genKeys,1:length(genSection));
gen.mpc2psse = containers.Map(1:length(genSection),genKeys);
indexMap.gen = gen;

%=====================================================
% Create branch index map
%=====================================================

% Get branch information section, which is 
ind1 = sectionIndices(4); % Non-xfmr branch data is section 5
ind2 = sectionIndices(5);
branchSection = rawFile(ind1:ind2-1);
branchSection = strsplit(strtrim(branchSection),'\n');

branchSection(1)=[]; % Discard 'BEGIN BRANCH DATA' text

% Get xfmr information section, which is 
ind1 = sectionIndices(5); % Xfmr data is section 6
ind2 = sectionIndices(6);
xfmrSection = strtrim(rawFile(ind1:ind2-1));
% Discard 'BEGIN XFMR DATA' text
indNL = regexp(xfmrSection,'\n','once');
xfmrSection = strtrim(xfmrSection(indNL:end));

xfmrSection = strrep(xfmrSection,char(10),'');
xfmrSection = strrep(strtrim(xfmrSection),char(13),',');
xfmrSection = strsplit(strtrim(xfmrSection),',');

% For xfmr records, cannot use '\n' as record separator,
% got to read all fields and divide by 43 (total records qty per xfmr)
xfmrSection = reshape(xfmrSection,43,[])';

% Build the map  
branchKeys = cell(length(branchSection)+size(xfmrSection,1),1);
% Add non-xfmr branches first
for i=1:length(branchSection)
    branchRow = branchSection{i};
    items = strsplit(branchRow,',');
    branchFrom = num2str(min(str2double(items{1}),...
                                              str2double(items{2})));
    branchTo = num2str(max(str2double(items{1}),...
                                              str2double(items{2})));
    branchKey = strcat(branchFrom,...
                       '-',branchTo,...
                       '-',strrep(strtrim(items{3}),'''',''));
	branchKeys{i} = branchKey;
end
% Add 2-winding xfmrs
j=1;
for i=1:size(xfmrSection,1)
    if ~strcmpi(strtrim(xfmrSection{i,3}),'0')
        % Ignore 3-winding xfmrs (there's not supposed to be any)
        continue
    end
    branchFrom = num2str(min(str2double(xfmrSection{i,1}),...
                                           str2double(xfmrSection{i,2})));
    branchTo = num2str(max(str2double(xfmrSection{i,1}),...
                                           str2double(xfmrSection{i,2})));
    xfmrKey = strcat(branchFrom,...
                       '-',branchTo,...
                       '-',strrep(strtrim(xfmrSection{i,4}),'''',''));
	branchKeys{length(branchSection)+j} = xfmrKey;
    j=j+1;
end
% If there were any 3-winding xfmrs, there are empty rows (the ones we
% skipped over). Delete them
branchKeys(length(branchSection)+j:end)=[];
% Create maps
branch.psse2mpc = containers.Map(branchKeys,1:(length(branchSection)+j-1));
branch.mpc2psse = containers.Map(1:(length(branchSection)+j-1),branchKeys);
indexMap.branch = branch;
%branchIndexMap = containers.Map(branchKeys,1:(length(branchSection)+j-1));
end

%=====================================================================
%=====================================================================
%=====================================================================

function [ branchOutMap, genOutMap ] = readParseCONFile( confile_name )
%READPARSECONFILE Reads and parses data in valid CON file
%   This function reads the data contained in a valid CON file and returns
%   it as maps in which the keys are the contingency label and the values
%   are the unique key of the branch or generator that is set out of
%   service in the contingency.
%
%[ branchOutMap, genOutMap ] = READPARSECONFILE( confile_name )

[~, ~, ext] = fileparts(confile_name);
if nargin<1
    error('Not enough arguments')
elseif ~strcmpi(ext,'.CON')
    disp(ext)
    error('Not a valid CON file')
end

% Read CON file
conFile = fileread(confile_name);
% Parse file into list of contingencies, defined as everything
% encompassed by 'CONTINGENCY' and 'END'
contingList = regexpi(conFile,'(?<=CONTINGENCY).*?[\n\r ]+END[\n\r]+','match');
branchOutKeys = {};
genOutKeys = {};
branchOutVals = {};
genOutVals = {};
for i=1:length(contingList)
    % For each contingency in the list
    contingency = strsplit(contingList{i},'\n');
    contingKey = strtrim(contingency{1});
    contingDef = strtrim(contingency{2});
    % Check if it's a branch-out contingency
    branchOut = regexpi(contingDef,...
        '(DISCONNECT|OPEN|TRIP) +(BRANCH|LINE) +FROM +BUS +(\d+) +TO +BUS +(\d+) +(CIRCUIT|CKT) +(''?[\w]{1,2}''?).*?','tokens');
    if ~isempty(branchOut)
        branchOut = branchOut{1};
        branchOutFrom = num2str(min(str2double(branchOut{3}),...
                                              str2double(branchOut{4})));
        branchOutTo = num2str(max(str2double(branchOut{3}),...
                                              str2double(branchOut{4})));
        branchOutVal = strcat(branchOutFrom,'-',branchOutTo,'-',...
            strrep(branchOut{6},'''',''));
        branchOutKeys{end+1} = contingKey;
        branchOutVals{end+1} = branchOutVal;
    end
    % Check if it's a gen-out contingency
    genOut = regexpi(contingDef,...
        '(REMOVE|TRIP) +(MACHINE|UNIT) +(''?[\w]{1,2}''?) +FROM +BUS +(\d+).*?','tokens');
    if ~isempty(genOut)
        genOut = genOut{1};
        genOutVal = strcat(genOut{4},'-',strrep(genOut{3},'''',''));
        genOutKeys{end+1} = contingKey;
        genOutVals{end+1} = genOutVal;
    end
end
if ~isempty(branchOutKeys)
    branchOutMap = containers.Map(branchOutKeys,branchOutVals);
else
    branchOutMap = [];
end
if ~isempty(genOutKeys)
    genOutMap = containers.Map(genOutKeys,genOutVals);
else
    genOutMap = [];
end
end

%=====================================================================
%=====================================================================
%=====================================================================

function genDroops = readParseINLFile(inlfile_name)
%READPARSEINLFILE
%   This function reads an INL file and retrieves the droop
%   participation factors for each generator. The results are 
%   returned as a map.
%   The map is of the Map class (similar to python dictionaries)
%   The map's keys are PSS/E's unique generator keys, built as
%   the concatenation of bus number and generator ID
%   e.g. '1-1'
%   The map's value is the generator droop participation factor.
%
%   genDroops  = READPARSEINLFILE( rawfile_name )
[~, ~, ext] = fileparts(inlfile_name);
if nargin<1
    error('Not enough arguments')
elseif ~strcmpi(ext,'.INL')
    disp(ext)
    error('Not a valid INL file')
end
    fid = fopen(inlfile_name,'r');
    genKeys = {};
    droops = [];
    while true
        tline = fgetl(fid);
        if ~ischar(tline)
            break
        end
        % Last line of file contains a '0'
        % Do not try to parse it
        if isempty(regexpi(tline,'[0Q] */*.*?','once')) ||...
                regexpi(tline,'[0Q] */*.*?','once')~=1
            items = strsplit(tline,',');
            genKey = strcat(strtrim(items{1}),...
                        '-',strrep(strtrim(items{2}),'''',''));
            genKeys{end+1} = genKey;
            droops = [droops,str2num(items{6})];
        end
    end
    fclose(fid);
    genDroops = containers.Map(genKeys,droops);
end

%=====================================================================
%=====================================================================
%=====================================================================

function [ genCostMap, maxNpairs ] = readParseRopFile( ropfile_name )
%READPARSEROPFILE Reads and parses ROP file to retrieve generator cost data
%   This function reads an ROP file and retrieves the generator cost data
%   The data is returned as a Map (similar to python's dictionaries)
%   with unique generator id's as keys and 1-by-2*N matrix containing
%   the coordinate pairs (P_i,C_i) for each generator as values.
%   The maximum number of pairs defining the piecewise linear cost function
%   can also be returned.
%   genCostMap  = READPARSEROPFILE( ropfile_name )
[~, ~, ext] = fileparts(ropfile_name);
if nargin<1
    error('Not enough arguments')
elseif ~strcmpi(ext,'.ROP')
    disp(ext)
    error('Not a valid ROP file')
end

% Read ROP file
ropFile = fileread(ropfile_name);

% Three dictionaries are read from the ROP file
% 1. First dictionaries contains active dispatch table for each generator
% 2. Second dict. contains piecewise cost table for each active dispatch
% table
% 3. Third dict. contains definition of each piecewise cost table, i.e. the
% pairs (P_i,C_i)
% See GOP formulation Appendix A.3 for further information

%=========================================================
% Read active dispatch table for each generator
%=========================================================
% Get indices of beginning of sections which are the ones starting with '0
% / [text]'
sectionIndices = regexp(ropFile,'0 */ *.*?\n');
ind1 = sectionIndices(5); % Generation dispatch data is section 5
ind2 = sectionIndices(6);
genDispSect = ropFile(ind1:ind2-1);
genDispSect = strsplit(strtrim(genDispSect),'\n');

genDispSect(1)=[]; % Discard 'begin Generator Dispatch data'
% Read data, parse it and make Map (dictionary)
genKeys = cell(length(genDispSect),1);
dispVals = zeros(length(genDispSect),1);
for i=1:length(genDispSect)
    genRow = genDispSect{i};
    items = strsplit(genRow,',');
    genKey = strcat(strtrim(items{1}),...
                    '-',strrep(strtrim(items{2}),'''',''));
	genKeys{i} = genKey;
    dispVals(i) = str2double(items{4});
end
genDispMap = containers.Map(genKeys,dispVals);


%=========================================================
% Read cost table for each active dispatch table
%=========================================================
% Create map where keys are active dispatch table numbers and
% values are cost tables according to Active Power Dispatch Tables

ind1 = sectionIndices(6); % Active power tables data is section 6
ind2 = sectionIndices(7);
actPowDispSect = ropFile(ind1:ind2-1);
actPowDispSect = strsplit(strtrim(actPowDispSect),'\n');

actPowDispSect(1)=[]; % Discard 'begin Active Power Dispatch ...'

% Create map where keys are table unique keys and
% values are piecewise linear cost table number 
actPowKeys = cell(length(actPowDispSect),1);
actDispVals = zeros(length(actPowDispSect),1);
for i=1:length(actPowDispSect)
    actDispRow = actPowDispSect{i};
    items = strsplit(actDispRow,',');
    actKey = str2double(items{1});
	actPowKeys{i} = actKey;
    actDispVals(i) = str2double(items{7});
end
actPowDispMap = containers.Map(actPowKeys,dispVals);

%=========================================================
% Read each cost table definition (data pairs)
%=========================================================
% Create map where keys are cost table numbers and
% values are a 1-by-2*N matrix with cost pairs (P_i,C_i)
ind1 = sectionIndices(10); % Cost curve tables is section 10
ind2 = sectionIndices(11);
pcwLinSect = ropFile(ind1:ind2-1);
pcwLinSect = strsplit(strtrim(pcwLinSect),'\n');

pcwLinSect(1)=[]; % Discard 'begin piece-wise linear...'
pcwKey = {};
pcwVal = {};

% For each table, read how many pairs of points it has (N)
nextRow = 1;
maxNpairs = 0;
while nextRow<=length(pcwLinSect)
    items = strsplit(pcwLinSect{nextRow},',');
    dispVal = str2double(items{1});
    pcwKey{end+1} = dispVal;
    nPoints = str2double(items{3});
    if nPoints > maxNpairs
        maxNpairs = nPoints;
    end
    costTable = zeros(1,2*nPoints);
    % For each pair of points, read each coordinate and 
    % set it along a one-axis array
    for i=1:nPoints
        costRow = strsplit(pcwLinSect{nextRow+i},',');
        costTable(2*(i-1)+1) = str2double(costRow{1});
        costTable(2*i) = str2double(costRow{2});
    end
    pcwVal{end+1} = costTable;
    nextRow = nextRow + nPoints +1;
end
pcwLinMap = containers.Map(pcwKey,pcwVal);

% Create Map with pcwLinCost data for each generator
clearvars -except *Map maxNpairs
genKeys = keys(genDispMap);
pcwLinVals = cell(length(genKeys),1);
for i=1:length(genKeys)
    pcwLinVals{i} = pcwLinMap(actPowDispMap(genDispMap(genKeys{i})));
end
genCostMap = containers.Map(genKeys,pcwLinVals);
end

function [NPairs,genCostPCW_corr] = getNPairs(genCostPCW)
%GETNPAIRS Get actual number of distinct pairs of points in piecewise
%linear generator cost definition. This functions gets a 1-by-2*N vector,
%where N is the assumed number of pairs, and checks if all X values are
%distinct. If they are not, they must be all equal to the last unique value.
%If so, the number of distinct pairs is returned in NPAIRS and the vector
%GENCOSTPCW_CORR is returned containing only the distinct pairs.
gcostXY = reshape(genCostPCW,2,[])';
[~,IA] = unique(gcostXY(:,1));
if length(IA)~=size(gcostXY,1)
    xEnd = gcostXY(IA(end)+1:end,1);
    yEnd = gcostXY(IA(end)+1:end,2);
    if ~all(xEnd==xEnd(1) & yEnd==yEnd(1))
        error('Error in definition of generator costs in ROP file.')
    end
end
NPairs = length(IA);
genCostPCW_corr = genCostPCW(1:2*NPairs);
end

