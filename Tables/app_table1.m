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
%%%

tableDataFile = 'table_data.mat';

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

if assembleData

nClips = 0;
nRecs = size(recording,2);
for file_idx = 1:nRecs
    nClips = nClips + size(recording(file_idx).data, 1);
end
totalData = cell(nClips, 6);
clipNum = 1;

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
    varN = '$N$';
    varDen = 'Density';
    varPolar = 'Polarization';
    varEnt = 'Entropy';
    varDenStd = 'STD Density';
    varPolarStd = 'STD Polarization';
    varEntStd = 'STD Entropy';
    varNames = { varBand, varN, varDen, varDenStd, varPolar, varPolarStd, varEnt, varEntStd};
    varTypes = { 'string','double','double','double','double','double','double','double'};
    
% initialize the table 
sz = [numel(rowNames) numel(varNames)];

collectiveTable = table('Size',sz,...
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
collectiveTable.('Band')(1:5) = bands(1:5);
collectiveTable.('Band')(idx_133) = bands(6);
collectiveTable.('Band')(idx_098) = bands(7);
collectiveTable.('Band')(idx_096) = bands(8);
collectiveTable.('Band')(idx_146) = bands(9);

for file_idx = 1:nRecs

if (1 <= file_idx) && (file_idx <= 2)
    data_all = recording(file_idx).data;
elseif (3 <= file_idx) && (file_idx <= 4)
    data_all = recording2(file_idx).data;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Data to be Plotted %%%
matNums = 1:numel(data_all(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(recording(file_idx).corners)
    scale = recording(file_idx).scale; % pi/cm
    trans = [];
    fieldDims = recording(file_idx).fieldDims / scale;
    area = fieldDims(2)*fieldDims(4)*0.01^2;
else
    cornersPix = recording(file_idx).corners;
    fieldDimsPix = recording(file_idx).fieldDims;
    [trans, scale, fieldDims, newR_A] = projTrans(cornersPix, fieldDimsPix); % in cm
    area = newR_A*0.01^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Up %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize data variable
allData = cell(length(matNums),6);

for m = matNums
    
    disp(['Now assembling data from vid' data_all{m,1}]);
    
    tic
    %%Here we assemble the data
    data_struct = data_all{m,2};

    [data_final, neighbors] = struct2data(data_struct);

    rng default %for randi used in function below
    % assemble the collective qualities from this data set
    [totalLocs, density, aveDir, polarization, entropy, vidStats]...
        = collectiveData(data_final,area);

    allData{m,1} = density;
    allData{m,2} = aveDir;
    allData{m,3} = polarization;
    allData{m,4} = entropy;
    allData{m,5} = vidStats;
    allData{m,6} = totalLocs;
    
    collectiveTable(nRecs+1+clipNum,...
        {varN, varDen, varDenStd, varPolar, varPolarStd, varEnt, varEntStd}) = num2cell(vidStats(1:end-2));
    totalData(clipNum,:) = {density, aveDir, polarization, entropy, vidStats, totalLocs};
    clipNum = clipNum + 1;
    
    fprintf('Assembling that data took %f seconds \n', toc)
end
    bandLocs = cell2mat(allData(:,6));
    bandMat = cell2mat(allData(:,[1 3 4]));
    bandMean = mean(bandMat,1,'omitnan');
    bandStd = std(bandMat,[],1,'omitnan');

    collectiveTable(1+file_idx,...
        {varN, varDen, varDenStd, varPolar, varPolarStd, varEnt, varEntStd}) = ...
        {sum(bandLocs), bandMean(1), bandStd(1), bandMean(2), bandStd(2), bandMean(3), bandStd(3)};
end

totalDataLocs = cell2mat(totalData(:,6));
totalDataMat = cell2mat(totalData(:,[1 3 4]));
totalDataMean = mean(totalDataMat,1, 'omitnan');
totalDataStd = std(totalDataMat,[],1, 'omitnan');

collectiveTable(1,...
    {varN, varDen, varDenStd, varPolar, varPolarStd, varEnt, varEntStd}) = ...
    {sum(totalDataLocs), totalDataMean(1), totalDataStd(1),...
    totalDataMean(2), totalDataStd(2), totalDataMean(3), totalDataStd(3)};

%%% Assign data from all the clips %%%
clipsData = cell2mat( totalData(:,5) );


save(tableDataFile,'collectiveTable',...
                   '-append')
else
   load(tableDataFile, 'collectiveTable') 
end

fprintf(['Total locusts points: %d | Mean Density: %f | \n'...
         '                               Mean Polarization: %f | Mean Entropy: %f \n'],...
         collectiveTable.(varN)(1), collectiveTable.(varDen)(1), collectiveTable.(varPolar)(1), collectiveTable.(varEnt)(1) )

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Associated Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [totalLocs, density, aveDir, polarization, entropy, vidStats] = collectiveData(data,area)

global  idx_x idx_y idx_flag...
        idx_speed idx_theta...
        idx_localSpeed idx_localStd idx_localMinMax...
        idx_state

[Nlocs, Nfeats, Ntimes] = size(data);

%%% Initialize Variables %%%
density = nan(Ntimes,1);
aveDir = nan(Ntimes,1);
polarization = nan(Ntimes,1);
entropy = nan(Ntimes,1);
totalLocs = 0;

for t = 1:Ntimes
    
    nLocusts = size(rmmissing(   data(:,idx_x,t)   ),1);
    totalLocs = totalLocs + nLocusts;
    
    %%% Density %%%
    density(t) = nLocusts/area;
    
    %%% Polarization & Average Direction %%%
    cVects = rmmissing( data(:,idx_speed,t).*exp(1i*data(:,idx_theta,t)) );
    cVects_norm = cVects./abs(cVects); % creates NaNs when abs(cVects) = 0
    % If we set cVects_norm = 0 for each cVects = 0, we would have angle(0)
    % included in our mean. Matlab says angle(0) = 0, while I say angle(0)
    % should be undefined, and therefore be excluded.
    % So we simply omit NaNs from out computation of the mean.
    aveDirVect = mean(cVects_norm, 'omitnan'); 
    aveDir(t) = angle(aveDirVect);
    polarization(t) = abs(aveDirVect);
    
    %%% Entropy %%%
    % based on Buhl 2011, Baldasarre 2009
    angles = rmmissing(   data(:,idx_theta,t)   );
    
    C = 72; % number of classes, each 5 degrees
    edges = linspace(0,2*pi,C);
    
    nAngles = size(angles,1)-1; %subtract one because we're about to exclude one
    
    if nAngles >= 0 % sometimes there are no angles!
        rng default %set outside function
        idx_ex = randi(nAngles+1); %choose a random to determine where theta = 0 is
        
        degreeShift = 360/(C*2);
        thetaRotate = angles(idx_ex)-degreeShift/180*pi;
        angles(idx_ex) = []; %exclude that one
        
        angles = angles - thetaRotate; %rotate everyone else
    else
        nAngles = 0;
    end
    
    
    
    [binCounts, edges] = histcounts(mod(angles,2*pi), edges);
    binCounts( binCounts == 0 ) = []; % throw out the zero binCounts
    % for debugging
%     if sum(binCounts) ~= nAngles
%         fprintf([   'Sum of BinCounts = %d != %d = nAngles \n'...
%                     'Locust removed: %d \n'...
%                     '--- \n'], sum(binCounts),nAngles,idx_ex)
%     end

    % for the normalization factor
    I = floor(nAngles/C);
    R = mod(nAngles,C);
        
    if nAngles>170 % to approximate log(factorial(n));
        logFact_nAngles = logFact(nAngles); 
        logFact_binCounts = logFact(binCounts);
        
        log_wmax = logFact_nAngles - sum( (C-R)*logFact(I) + R*logFact(I+1) );
        k = 1/log_wmax;
    else % or to compute directly
        logFact_nAngles = log(factorial(nAngles));
        logFact_binCounts = log(factorial(binCounts));
        
        wmax = factorial(nAngles)/(   factorial(I+1)^(R)*factorial(I)^(C-R)   );
        k = 1/log(  wmax  );
    end
    % for debugging
%     if isnan(logFact_nAngles)
%         error('logFact_nAngles was NaN')
%     end
%     if isnan( sum( logFact_binCounts) ) 
%         error('the sum of logFact_binCounts was NaN')
%     end
    
        % for debugging
%     if isnan(wmax)
%         error('wmax was NaN')
%     end
%     if isnan(k) || k == 0 || k == Inf
%         error('k was NaN')
%     end
    
    entropy(t) = k*( logFact_nAngles - sum( logFact_binCounts ) );
    % for debugging
%     if isnan( entropy(t) ) || entropy(t) == Inf
%         error('entropy was NaN')
%     end
end

aveDir = unwrap( aveDir, [], 1 ); %removes discontinuities from -pi to pi
% was in vidValues mean(aveDir,'omitNaN')
vidStats = [totalLocs, mean(density,'omitNaN'), std(density,'omitNaN'),...
                        mean(polarization,'omitNaN'), std(polarization,'omitNaN'),...
                        mean(entropy,'omitNaN'), std(entropy,'omitNaN'),...
                        circ_mean( rmmissing(aveDir) ), circ_std( rmmissing(aveDir) )];

% output = density, aveDir, polarization, entropy, vidValues

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lnofnFact = logFact(n)
        lnofnFact = (n + 1/2).*log(n) - n + log(sqrt(2*pi));
end