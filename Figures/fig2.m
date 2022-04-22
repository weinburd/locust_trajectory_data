%plots density, polarization, and entropy as in Figure 2

close all
% clear all % ensure we're loading new data
addpath('../Data/',...
    '../Functions',...
    '../Functions/packages/CircStat')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

figDataFile = 'fig_data.mat';

%%% Options %%%
saveFigs = 1;
assembleData = 0; %and save it to figDataFile
    %else loads data from figDataFile
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 1 %%%
number = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

nClips = 0;
nRecs = size(recording,2);
for file_idx = 1:nRecs
    nClips = nClips + size(recording(file_idx).data, 1);
end
totalData = cell(nClips, 6);
clipNum = 1;

%file_idx = 1; %vid133
%file_idx = 2; %vid098
%file_idx = 3; %vid096
%file_idx = 4; %vid146

%%% Set up the Table %%%
% tableName = 'collectiveTable';
% collectiveTable = initTable(tableName,[],0); %[] loads table, recording creates a new one
varN = '$N$';
varDen = 'Density';
varPolar = 'Polarization';
varEnt = 'Entropy';
varDenStd = 'STD Density';
varPolarStd = 'STD Polarization';
varEntStd = 'STD Entropy';

%collectiveTable = collectiveTable(:,{varN, varDen, varPolar, varEnt});
varTypes = {'double','double','double','double','double','double','double'};
varNames = {varN, varDen, varDenStd, varPolar, varPolarStd, varEnt, varEntStd};
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
    sz = [numel(rowNames) numel(varNames)];

collectiveTable = table('Size',sz,...
                    'VariableTypes',varTypes,...
                    'RowNames',rowNames,...
                    'VariableNames',varNames);
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
    %collectiveTable(clipNum+1+nRecs,1:4) = num2cell(vidValues);
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

save(figDataFile,   'totalData', 'collectiveTable',...
                    '-append')

else
    load(figDataFile, 'totalData', 'collectiveTable')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yLabel = 'Polarization';
h_polar = plotColPoints(number+4, totalData(:,[1 3]), collectiveTable.Band(6:end),yLabel);
yLabel = 'Entropy Index';
h_entrop = plotColPoints(number+5, totalData(:,[1 4]), collectiveTable.Band(6:end),yLabel);

fprintf(['Total locusts points: %d | Mean Density: %f | \n'...
         '                               Mean Polarization: %f | Mean Entropy: %f \n'],...
         collectiveTable.(varN)(1), collectiveTable.(varDen)(1), collectiveTable.(varPolar)(1), collectiveTable.(varEnt)(1) )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFigs
    
    figPath = 'figs/';
    
    filename = 'fig2_polar';
    print(h_polar,[figPath filename, '.eps'],'-depsc')
    
    filename = 'fig2_entrop';
    print(h_entrop,[figPath filename, '.eps'],'-depsc')
    
end

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


end

function h = plotColPoints(figNum, totalData, bands, yLabel)

data = struct();

band_list = unique( bands, 'stable');
nBands = numel(band_list);

for i = 1:nBands
    idx = find( cellfun(@(c) strcmp(c,band_list{i}), bands) );

    data(i).x = cell2mat(totalData(idx,1));
    data(i).y = cell2mat(totalData(idx,2));
end

factor = 2;

h = figure(figNum);
set(h,'Units','Inches');
pos = get(h, 'Position');
set(h,'PaperPositionMode','Manual') % Setting this to 'manual' unlinks 'Position' (what gets displayed) and 'PaperPosition' (what gets printed)
set(h,'PaperPosition',[ 0 0 4 3]*factor);
set(h,'Position',[ 0 0 4 3]*factor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure Options %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
titlesize = 12*factor;
axislabelsize = 10*factor;
fontSize = axislabelsize;
ticklabelsize = 8*factor;
subfiglabelsize = 12*factor;

%linewidth
lwidth = 1.5*factor;
mrkSize = 5*factor;

color = lines(6); % default matlab colors...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scatter Plot! %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alph = .4;
sz = mrkSize;
for i = 1:size(data,2)
    g = scatter(data(i).x(1:5:end), data(i).y(1:5:end),sz,color(i,:), 'o','filled',...
        'MarkerEdgeAlpha',alph, 'MarkerFaceAlpha',alph);
    hold on
end
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Labels %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca, 'FontSize', ticklabelsize)
xlabel('Density (locusts/meter$^2$)','FontSize',axislabelsize)
ylabel(yLabel,'FontSize',axislabelsize)

meanX = mean(cell2mat(totalData(:,1)),"omitnan");
meanY = mean(cell2mat(totalData(:,2)),"omitnan");

if strcmp(yLabel,'Polarization')
    ymin = 0.45;
    %break_axis('bottom_reset','0','position',0.47)

    hold on
    p1 = plot([meanX meanX],[ymin 1],'k--','Linewidth',lwidth);
    % second p2 added below after the legend
    hold off

    leg = legend( '', 'Location', 'southeast' );
elseif strcmp(yLabel,'Entropy Index')
    ymin = 0.56;    
    
    hold on
    p1 = plot([meanX meanX],[ymin 1],'k--','Linewidth',lwidth);
    % second p2 added below after the legend
    hold off

    leg = legend( '', 'Location', 'northeast' );
end
% change legend colors to take out transparency
band_list{4} = 'band 4';
for i = [1 3 2 4]
    hold on
    scatter([NaN], [NaN], sz, color(i,:), 'o','filled',...
            'DisplayName',band_list{i});
end
p2 = plot([0 425],[meanY meanY], 'k--','Linewidth',lwidth,'DisplayName','mean values');
uistack(p1,'bottom')
uistack(p2,'bottom')
hold off


axis([0 425 ymin 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Touch Up %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %remove whitespace
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    %%% Add a y-axis break %%%
    pos = get(ax, 'position');

    ytics = ax.YTick;
    yticLabels = ax.YTickLabel;

    position = (ytics(1)+ymin)/2;
    lims = get(gca, 'ylim');
    cent = pos(2) + pos(4) * (position - lims(1))/(diff(lims));

    % define the extents of the bar
    X1 = [pos(1) - 0.005; 
          pos(1) + 0.005];
    X2 = [pos(1) - 0.015; 
          pos(1) + 0.015];
    Y = [cent - 0.005; 
         cent + 0.005];
     
    % Make the bar and ticks
    h1 = annotation(h, 'rectangle', [X1(1), Y(1), diff(X1), diff(Y)]);
    annotation(h, 'line', [X2(1) X2(2)], [Y(1) Y(1)]);
    annotation(h, 'line', [X2(1) X2(2)], [Y(2) Y(2)]);
    % reset the color to white
    h1.FaceColor = 'w';
    h1.Color = 'w';

    ax.YTick = [ymin ytics];
    ax.YTickLabel = vertcat('0', yticLabels);
end

function lnofnFact = logFact(n)
        lnofnFact = (n + 1/2).*log(n) - n + log(sqrt(2*pi));
end