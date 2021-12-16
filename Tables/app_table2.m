close all
% clear all % ensure we're loading new data
addpath('../Functions',...
        '../Functions/packages/CircStat',...
        '../Data/')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

global  idx_x idx_y idx_flag...
        idx_speed idx_theta...
        idx_localSpeed idx_localStd idx_localMinMax...
        idx_state

%%% Options %%%
assembleData = 1; %set to 0 to load existing table data
printStats = 1; %set to 0 to surpress printing motion stats as you go
%%%

tableDataFile = 'table_data.mat';


if assembleData

if ~exist('recording','var')
    disp('Now loading data...')
    tic;
    load('data_recording.mat')
    fprintf('That took %f seconds \n', toc)
end
if ~exist('recording2','var')
    disp('Now loading data...')
    tic;
    load('data_recording2.mat')
    fprintf('That took %f seconds \n', toc)
end
% the loaded file includes the following variables:
%   1) recording = a struct array
%   2-10) idx_* = the global feature indices for data_final

%%% Set up the Tables %%%
%Construct a new table, using metadata from recording
    nRecs = size(recording,2);
    
    %Process the row names
    rowNames = {'All'};
    for j = 1:nRecs
        recName = strrep(recording(j).file, '_', '\_');
        rowNames{end+1} = recName;
    end
    for j = 1:nRecs
        clips = recording(j).data(:,1);
        clips = cellfun(@(C) strrep( C , '_', '\_'), clips, 'UniformOutput',false);
        rowNames(end+1:end+numel(clips)) = clips;
    end
    rowNames = rowNames';
    
    %Generate the column names
    varBand = 'Band';
    varNSpdAll = '$N_sp$ (all)';
    varPerStop = '\% Stop';
    varPerWalk = '\% Walk';
    varPerHop = '\% Hop';
    varSpeedAll = 'Speed All';
    varSpeedAllStd = 'Speed All Std';
    varSpeedStop = 'Speed Stop';
    varSpeedStopStd = 'Speed Stop Std';
    varSpeedWalk = 'Speed Walk';
    varSpeedWalkStd = 'Speed Walk Std';
    varSpeedHop = 'Speed Hop';
    varSpeedHopStd = 'Speed Hop Std';
    varNames = {varBand, varNSpdAll, varPerStop, varPerWalk, varPerHop,...
                        varSpeedAll, varSpeedAllStd, varSpeedStop, varSpeedStopStd,...
                        varSpeedWalk, varSpeedWalkStd, varSpeedHop, varSpeedHopStd};
    varTypes = { 'string','double','double','double','double',...
                          'double','double','double','double',...
                          'double','double','double','double'};
    
% initialize the table 
sz = [numel(rowNames) numel(varNames)];

motionTable = table('Size',sz,...
                'VariableTypes',varTypes,...
                'RowNames',rowNames,...
                'VariableNames',varNames);
            
% Fill in the 'Band' column with our metadata
% band6 is denoted band4 in the paper
idx_133 = 2+nRecs:2+nRecs+5;
idx_098 = idx_133(end)+1:idx_133(end)+7;
idx_096 = idx_098(end)+1:idx_098(end)+7;
idx_146 = idx_096(end)+1:idx_096(end)+7;
bands = {'all bands', 'all band 1', 'all band 3', 'all band 2', 'all band 6', 'band 1', 'band 3', 'band 2', 'band 6'};
motionTable.('Band')(1:5) = bands(1:5);
motionTable.('Band')(idx_133) = bands(6);
motionTable.('Band')(idx_098) = bands(7);
motionTable.('Band')(idx_096) = bands(8);
motionTable.('Band')(idx_146) = bands(9);


nClips = 0;
for file_idx = 1:4
    nClips = nClips + size(recording(file_idx).data, 1);
end
totalSpeeds = cell(nClips, 4);
totalAngles = cell(nClips, 4);
totalStates = cell(nClips, 4);
clipNum = 1;

for file_idx = 1:4

if (1 <= file_idx) && (file_idx <= 2)
    data_all = recording(file_idx).data;
elseif (3 <= file_idx) && (file_idx <= 4)
    data_all = recording2(file_idx).data;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Data to be Plotted %%%
matNums = 1:numel(data_all(:,1));
%matNums = 1; % for testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Up %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize data variable
dataSpeeds = cell(length(matNums),4);
dataAngles = cell(length(matNums),4);
dataStates = cell(length(matNums),4);

for m = matNums
    
    disp(['Now assembling data from vid' data_all{m,1}]);
    
    tic
    %%Here we assemble the data
    data_struct = data_all{m,2};

    [data_final, neighbors] = struct2data(data_struct);

    % assemble the individual data for this data set
    speeds = indivData(data_final,idx_speed);
    angles = indivData(data_final,idx_theta);
    states = indivData(data_final,idx_state);
    
    dataSpeeds(m,:) = speeds;
    dataAngles(m,:) = angles;
    dataStates(m,:) = states;
    
    [Nstates, ~] = computeStats(states,'state', printStats);
    motionTable.(varNSpdAll)(1+nRecs+clipNum) = Nstates{1};
    motionTable.('\% Stop')(1+nRecs+clipNum) = 100*Nstates{2}/Nstates{1};
    motionTable.('\% Walk')(1+nRecs+clipNum) = 100*Nstates{3}/Nstates{1};
    motionTable.('\% Hop')(1+nRecs+clipNum) =100*Nstates{4}/Nstates{1};

    [Nspeeds, stats_speed] = computeStats(speeds,'speed', printStats);
    motionTable.('Speed All')(1+nRecs+clipNum) = stats_speed{1,1}(1);
    motionTable.('Speed All Std')(1+nRecs+clipNum) = stats_speed{1,1}(3);
    motionTable.('Speed Stop')(1+nRecs+clipNum) = stats_speed{1,2}(1);
    motionTable.('Speed Stop Std')(1+nRecs+clipNum) = stats_speed{1,2}(3);
    motionTable.('Speed Walk')(1+nRecs+clipNum) = stats_speed{1,3}(1);
    motionTable.('Speed Walk Std')(1+nRecs+clipNum) = stats_speed{1,3}(3); 
    motionTable.('Speed Hop')(1+nRecs+clipNum) = stats_speed{1,4}(1);
    motionTable.('Speed Hop Std')(1+nRecs+clipNum) = stats_speed{1,4}(3);
    
    % to shift the angles
    meanA = circ_mean(angles{1});
    
    totalSpeeds(clipNum,:) = speeds;
    totalAngles(clipNum,:) = cellfun(@(M) M-meanA+pi/2, angles, 'UniformOutput', false);
    totalStates(clipNum,:) = states;
    clipNum = clipNum + 1;
    
    %%%  NEED TO CONFIRM THIS IS WORKING AS INTENDED %%%
    %%Here we calculate the intervals over which a locust is in each one of
    %%these states. 
    %stateDurations = StateDrt(data_final);
    %%%  NEED TO CONFIRM THIS IS WORKING AS INTENDED %%%

    fprintf('Assembling that data took %f seconds \n', toc)
end

% Compute stats for this recording
[Nstates,~] = computeStats(dataStates,'state', printStats);
motionTable.(varNSpdAll)(1+file_idx) = Nstates{1};
motionTable.('\% Stop')(1+file_idx) = 100*Nstates{2}/Nstates{1};
motionTable.('\% Walk')(1+file_idx) = 100*Nstates{3}/Nstates{1};
motionTable.('\% Hop')(1+file_idx) =100*Nstates{4}/Nstates{1};

[Nspeeds, stats_speed] = computeStats(dataSpeeds,'speed', printStats);
motionTable.('Speed All')(1+file_idx) = stats_speed{1,1}(1);
motionTable.('Speed All Std')(1+file_idx) = stats_speed{1,1}(3);
motionTable.('Speed Stop')(1+file_idx) = stats_speed{1,2}(1);
motionTable.('Speed Stop Std')(1+file_idx) = stats_speed{1,2}(3);
motionTable.('Speed Walk')(1+file_idx) = stats_speed{1,3}(1);
motionTable.('Speed Walk Std')(1+file_idx) = stats_speed{1,3}(3); 
motionTable.('Speed Hop')(1+file_idx) = stats_speed{1,4}(1);
motionTable.('Speed Hop Std')(1+file_idx) = stats_speed{1,4}(3);
    
end

% Compute stats for this recording
[Nstates,~] = computeStats(totalStates,'state', printStats);
motionTable.(varNSpdAll)(1) = Nstates{1};
motionTable.('\% Stop')(1) = 100*Nstates{2}/Nstates{1};
motionTable.('\% Walk')(1) = 100*Nstates{3}/Nstates{1};
motionTable.('\% Hop')(1) = 100*Nstates{4}/Nstates{1};

[Nspeeds, stats_speed] = computeStats(totalSpeeds,'speed', printStats);
motionTable.('Speed All')(1) = stats_speed{1,1}(1);
motionTable.('Speed All Std')(1) = stats_speed{1,1}(3);
motionTable.('Speed Stop')(1) = stats_speed{1,2}(1);
motionTable.('Speed Stop Std')(1) = stats_speed{1,2}(3);
motionTable.('Speed Walk')(1) = stats_speed{1,3}(1);
motionTable.('Speed Walk Std')(1) = stats_speed{1,3}(3); 
motionTable.('Speed Hop')(1) = stats_speed{1,4}(1);
motionTable.('Speed Hop Std')(1) = stats_speed{1,4}(3);

save(tableDataFile,'motionTable',...
                   '-append')
else
   load(tableDataFile, 'motionTable') 
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Associated Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataByState = indivData(M,idx_feature)
    %%% inputs
    %feature = nlocusts by 1 by ntimesteps matrix
    %%% outputs
    % dataByState = a 1 by 4 cell array of data separated by motion state: 
    %               {all, stop, Walk, hop}
    
    global idx_state
    
    newData = M(:,idx_feature,2:end); %leave off first frame since no motion state is not assigned

    states = M(:,idx_state,2:end);
    
    dataByState = { newData, newData( states==0 ),...
                    newData( states==1 ), newData( states==2 )};
    for column = 1:size(dataByState,2)
        data = dataByState{column};
        data = rmmissing( data(:) );
        dataByState{column} = data;
    end
end 

function wideData = shapeData(thinData, idx_data, s)
    % recreates old shape
    newData = NaN(length(idx_data),1);
    newData(idx_data) = thinData;
    wideData = reshape(newData, s(1), 1, s(3)-1);
end

function [N, stats] = computeStats(data, type, printFlag)
%returns nSpeeds, meanSpeeds, stdSpeeds, nAngles, circm, stdAngles, circPolar, stdPolar
% Note: because we are shifting all angles, circm = mean angle = pi/2.

stats = cell(1,4); %initialize variable
N = cell(1,4);

motion = {'all', 'stop', 'walk', 'hop'};

if strcmp(type,'speed')
    for column = 1:size(data,2)
        speeds = cell2mat(data(:,column));
        nSpeeds = size(speeds, 1);
        meanSpeed = mean(speeds);
        medSpeed = median(speeds);
        stdSpeed = std(speeds);
        stats{column} = [meanSpeed medSpeed stdSpeed];
        if printFlag
            fprintf([
             'Number of ' motion{column} ' Speeds: %d \n'...
             '     Mean: %f | Median: %f | STD: %f \n'...
             ],...
             nSpeeds, meanSpeed, medSpeed, stdSpeed)
        end
        N{column} = nSpeeds;
    end
elseif strcmp(type,'angle')
    for column = 1:size(data,2)
        angles = cell2mat(data(:,column));
        % Computing some circular statistics
        nAngles = size(angles,1);
        meanAngle = circ_mean(angles);
        bootmean = bootstrp(100,@circ_mean,angles);
        stdAngle = std(bootmean);
        %%% IS STD RIGHT?!!!!
        stats{column} = [meanAngle stdAngle];
        if printFlag
            fprintf([
             'Number of ' motion{column} ' Angles: %d \n'...
             '     Mean: %f | STD: %f \n'...
             ],...
             nAngles, meanAngle, stdAngle)
        end
        N{column} = nAngles;
    end
elseif strcmp(type,'state')
        nLoc = size( cell2mat(data(:,1)), 1);
        nStop = size( cell2mat(data(:,2)), 1);
        nCrawl = size( cell2mat(data(:,3)), 1);
        nHop = size( cell2mat(data(:,4)), 1);
        N = {nLoc, nStop, nCrawl, nHop};
        if printFlag
            fprintf([
             'Number of ' motion{1} ' States: %d \n'...
             '     Stop: %d | Crawl %d | Hop: %d \n'...
             ],...
             nLoc, nStop, nCrawl, nHop)
        end
        stats = {'no stats for movement state'};
end

end

function stateDurations = StateDrt(M)
    stateDurations = NaN(size(M,1)*size(M,3),3);
    states_here = M(:,idx_state,:);
    states_here = permute(states_here,[3,2,1]);
    states_here = transpose(states_here(:));
    for st = 0:2
        start1 = strfind([0,states_here==st],[0 1]);
        end1 = strfind([states_here==st,0],[1 0]);
        d = transpose(end1 - start1 + 1);
        d = d.*(1/25);
        stateDurations(1:size(d,1),st+1) = d;
    end
end