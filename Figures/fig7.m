%bar graphs of the proportion of neighbors to the side of the focal locust

close all
%clear all % ensure we're loading new data

addpath('../Data/',...
        '../Functions',...
            '../Functions/packages/CircStat')

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

global  idx_x idx_y idx_flag...
        idx_speed idx_theta...
        idx_localSpeed idx_localStd idx_localMinMax...
        idx_state...
        numNeighbors numFocal

figDataFile = 'fig_data_reshape.mat';

% To restore original...
% edit reshapeFactor, dx, dy, numNeighbors, and uncomment caxis([0 1.6])

%%% Options %%%
saveFigs = 1;
assembleData = 0; %and save it to figDataFile
    %else loads data from figDataFile

% the max number of neighbors around each focal individual
% can be an integer or 'all', but 'all' cannot distinguish focalState
numNeighbors = 1; % 100 gets all of them % 'all';
d = 4; % 7 cm max radius for angles
% for comparing rescaled must set this at maxmum to 14/3 because the
% original data goes out to a radius of 14 cm and the rescaling factor is 1/3

% for RESHAPING Data
% suppose a locust is approximatly 5mm wide and 15mm long
reshapeData = {'none', [0 0]}; 
%reshapeData = {'rescale', [.5 1.5]}; 
%reshapeData = {'worst', [0 1.5]}; 
% 'rescale', [.5 1.5] will rescale by: 0.5*[1/.5 1/1.5]
% or
% 'subtract', [.5 1.5] will subtract the distance within ellipse
% with width = .5 and length = 1.5cm, at the angle of each relative
% neighbor
% 'subtractVert', [.5 1.5] will subtract the vertical distance within
% ellipse (same ellipse as 'subtract')
% or
% 'worst', [0, 1.5] will subtract a vertical distance determined by
% modeling the locusts as linesegments with length 1.5cm and finding the
% minimum distance between the focal's head and the neighbor.
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 20s %%%
number = 20;
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

totalData = collectData(recording, recording2, d, reshapeData);

rescaleData = collectData(recording, recording2, d, {'rescale', [.5 1.5]});

% totalData = 27 x 16 cell array
%           = clips x [relative neighbor positions, 
%                      relative neighbor angles, 
%                      proportion of side neighbors,
%                      mean density] 
%                      for each of all, stop, crawl, hop

clear box

vmc = 4; % number of features in totalData for each motion state, defined in collectData()
% convert to a matrix array
allData = cell2mat(totalData(:,3));
stopData = cell2mat(totalData(:,3+vmc));
crawlData = cell2mat(totalData(:,3+2*vmc));
hopData = cell2mat(totalData(:,3+3*vmc));
box.prop = [allData.prop stopData.prop crawlData.prop hopData.prop]';
box.motion(1:27) = "all";
box.motion(28:54) = "stop";
box.motion(55:81) = "walk";
box.motion(82:108) = "hop";
box.type(1:108) = "original";

allRescale = cell2mat(rescaleData(:,3));
stopRescale = cell2mat(rescaleData(:,3+vmc));
crawlRescale = cell2mat(rescaleData(:,3+2*vmc));
hopRescale = cell2mat(rescaleData(:,3+3*vmc));
box.prop = [box.prop; [allRescale.prop stopRescale.prop crawlRescale.prop hopRescale.prop]'];
box.motion = [box.motion box.motion];
box.type(109:216) = "rescaled";

scatter.dens = cell2mat(totalData(:,4));
scatter.prop_orig = [allData.prop];
scatter.prop_rescale = [allRescale.prop];

bodyShapedData = matfile(figDataFile, 'Writable',true);

bodyShapedData.("boxplot") = box;
bodyShapedData.("scatter") = scatter;
bodyShapedData.("numNei") = numNeighbors;
bodyShapedData.("minDist") = d;

% save(figDataFile,   'binN', 'binCtrs',...
%                     'binSizes','-append')

else
    bodyShapedData = matfile(figDataFile);
    box = bodyShapedData.boxplot;
    scatter = bodyShapedData.scatter;
    numNeighbors = bodyShapedData.numNei;
    d = bodyShapedData.minDist;
    
%     load(figDataFile,   'binN', 'binCtrs',...
%                         'binSizes')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

mytitle = sprintf("Proportion Side Neighbors ($K = %d$, $d = %d$ cm)", numNeighbors, d );
mytitle = sprintf("Proportion of Side Neighbors");

figbox = figure(19);
% some categorical stuff here
catmotion = categorical(box.motion, ["all", "space", "stop","space1",  "walk","space2",  "hop"]);
b = boxchart( catmotion, box.prop, "GroupByColor", box.type, "boxwidth", 0.8 );
hold on
ax = gca(figbox);
plot(ax.XLim, [0.5 0.5], '--k')
hold off
ax.XTickLabel = ["all", " ","stop"," ", "walk"," ", "hop"];
ax.XAxis.FontSize = 20;
ylabel("Proportion of Side Neighbors", "FontSize", 20)
legend(["original", "rescaled"], "FontSize", 14)
title( mytitle, "FontSize", 24 )

color = lines(6);

bands(1:6,1) = "Band 1";
bands(7:13,1) = "Band 3";
bands(14:20,1) = "Band 2";
bands(21:27,1) = "Band 4";

figscat = figure(57);
gdot = gscatter( scatter.dens, scatter.prop_orig, bands, color(1:4,:), '.', 36);
hold on
gdothide = gscatter( [0], [0], {'dot'}, 'k', '.', 36);
gcrosshide = gscatter([0], [0], {'cross'}, 'k', 'x', 12);
gcrosshide(1).LineWidth = 3;
gcross = gscatter( scatter.dens, scatter.prop_rescale, bands, color(1:4,:), 'x', 12);
% add a black dotted line at 0.5
ax = gca(figscat);
plot(ax.XLim, [0.5 0.5], '--k')
hold off
for i = 1:numel(gcross)
    gcross(i).LineWidth = 3;
end
xlim([0 300])
ylim([0.475 0.65])
xlabel("Mean Density", "FontSize", 20)
ylabel("Proportion of Side Neighbors", "FontSize", 20)
legend([gdothide; gcrosshide; gdot],["original", "rescaled", "Band 1", "Band 3", "Band 2", "Band 4"], "FontSize", 14)
title( mytitle, "FontSize", 24 )

% TODO: 
%   NO Make y-axis symmetric about 0.50
%   DONE add a dotted line at proportion = 0.5
%   DONE do reshape data as well as original data
%   DONE experiment with "d" and ensure we're not including too many isometric
%       neighbors far out
%   DONE modify the save statement to include this stuff
%   DONE in scatterplot, color dots by band#
%   test reshaping on simulated (uniform) data?

fprintf('Plotting all the figures took %f seconds \n', toc)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Figure, Table, LaTeX %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFigs
    
    figPath = 'figs/';
    
    filename = sprintf('fig7_box_K%d_d%d', numNeighbors,d);
    print(figbox,[figPath filename, '.eps'],'-depsc')
    
    filename = sprintf('fig7_scatter_K%d_d%d', numNeighbors,d);
    print(figscat,[figPath filename, '.eps'],'-depsc')
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Associated Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Moved to separate function file neighPositions.m
% function[thisData, varargout] = neighPositions(data, neighbors, varargin)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [totalLocs, density vidStats] = collectiveData(data,area)
%modified from Table 1 (I think?)
%[totalLocs, density, aveDir, polarization, entropy, vidStats] = collectiveData(data,area)

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
    
%     %%% Polarization & Average Direction %%%
%     cVects = rmmissing( data(:,idx_speed,t).*exp(1i*data(:,idx_theta,t)) );
%     cVects_norm = cVects./abs(cVects); % creates NaNs when abs(cVects) = 0
%     % If we set cVects_norm = 0 for each cVects = 0, we would have angle(0)
%     % included in our mean. Matlab says angle(0) = 0, while I say angle(0)
%     % should be undefined, and therefore be excluded.
%     % So we simply omit NaNs from out computation of the mean.
%     aveDirVect = mean(cVects_norm, 'omitnan'); 
%     aveDir(t) = angle(aveDirVect);
%     polarization(t) = abs(aveDirVect);
    
%     %%% Entropy %%%
%     % based on Buhl 2011, Baldasarre 2009
%     angles = rmmissing(   data(:,idx_theta,t)   );
%     
%     C = 72; % number of classes, each 5 degrees
%     edges = linspace(0,2*pi,C);
%     
%     nAngles = size(angles,1)-1; %subtract one because we're about to exclude one
%     
%     if nAngles >= 0 % sometimes there are no angles!
%         rng default %set outside function
%         idx_ex = randi(nAngles+1); %choose a random to determine where theta = 0 is
%         
%         degreeShift = 360/(C*2);
%         thetaRotate = angles(idx_ex)-degreeShift/180*pi;
%         angles(idx_ex) = []; %exclude that one
%         
%         angles = angles - thetaRotate; %rotate everyone else
%     else
%         nAngles = 0;
%     end
%     
%     
%     
%     [binCounts, edges] = histcounts(mod(angles,2*pi), edges);
%     binCounts( binCounts == 0 ) = []; % throw out the zero binCounts
%     % for debugging
% %     if sum(binCounts) ~= nAngles
% %         fprintf([   'Sum of BinCounts = %d != %d = nAngles \n'...
% %                     'Locust removed: %d \n'...
% %                     '--- \n'], sum(binCounts),nAngles,idx_ex)
% %     end
% 
%     % for the normalization factor
%     I = floor(nAngles/C);
%     R = mod(nAngles,C);
%         
%     if nAngles>170 % to approximate log(factorial(n));
%         logFact_nAngles = logFact(nAngles); 
%         logFact_binCounts = logFact(binCounts);
%         
%         log_wmax = logFact_nAngles - sum( (C-R)*logFact(I) + R*logFact(I+1) );
%         k = 1/log_wmax;
%     else % or to compute directly
%         logFact_nAngles = log(factorial(nAngles));
%         logFact_binCounts = log(factorial(binCounts));
%         
%         wmax = factorial(nAngles)/(   factorial(I+1)^(R)*factorial(I)^(C-R)   );
%         k = 1/log(  wmax  );
%     end
%     
%     
%     entropy(t) = k*( logFact_nAngles - sum( logFact_binCounts ) );
% 
end

% aveDir = unwrap( aveDir, [], 1 ); %removes discontinuities from -pi to pi
% was in vidValues mean(aveDir,'omitNaN')
vidStats = [totalLocs, mean(density,'omitNaN'), std(density,'omitNaN'),...
                          ];
%                         mean(polarization,'omitNaN'), std(polarization,'omitNaN'),...
%                         mean(entropy,'omitNaN'), std(entropy,'omitNaN'),...
%                         circ_mean( rmmissing(aveDir) ), circ_std( rmmissing(aveDir) )];

% output = density, aveDir, polarization, entropy, vidValues

end

function totalData = collectData(recording, recording2, d, reshapeData)

global  idx_x idx_y idx_flag...
        idx_speed idx_theta...
        idx_localSpeed idx_localStd idx_localMinMax...
        idx_state...
        numNeighbors numFocal

nClips = 0;
nRecs = size(recording,2);
for file_idx = 1:nRecs
    nClips = nClips + size(recording(file_idx).data, 1);
end
totalData = cell(nClips, 16);
clipNum = 1;

%file_idx = 1; %vid133
%file_idx = 2; %vid098
%file_idx = 3; %vid096
%file_idx = 4; %vid146

for file_idx = 1:nRecs

    if (1 <= file_idx) && (file_idx <= 2)
        data_all = recording(file_idx).data;
    elseif (3 <= file_idx) && (file_idx <= 4)
        data_all = recording2(file_idx).data;
    end

    tic

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Data to be Plotted %%%
    matNums = 1:numel(data_all(:,1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(recording(file_idx).corners)
        scale = recording(file_idx).scale; % pi/cm
        trans = [];
        fieldDims = recording(file_idx).fieldDims / scale;
        w = fieldDims(2)-fieldDims(1);
        h = fieldDims(4)-fieldDims(3);
        area = w*h*0.01^2;
    else
        cornersPix = recording(file_idx).corners;
        fieldDimsPix = recording(file_idx).fieldDims;
        [trans, scale, fieldDims, newR_A] = projTrans(cornersPix, fieldDimsPix); % in cm
        area = newR_A*0.01^2;
    end

    vmc = 4; 
    allData = cell(length(matNums),vmc);
    stopData = cell(length(matNums),vmc);
    crawlData = cell(length(matNums),vmc);
    hopData = cell(length(matNums),vmc);
    numFocal = 0;

    for m = matNums

        disp(['Now assembling data from vid' data_all{m,1}]);

        tic
        % where the magic happens
        if isa(numNeighbors,'double')
            %%Here we assemble the data
            data_struct = data_all{m,2};

            [data_final, neighbors] = struct2data(data_struct);
            % extract position, orientation, and motion state
            data = [data_final(:,idx_x:idx_y,:) data_final(:,idx_theta,:) data_final(:,idx_state,:)];

            %%% for debugging... %%%
            inDebug.count = 0;
            inDebug.countnan = 0;
            inDebug.countempty = 0;
            inDebug.countnoneighs = 0;

            % assemble a new neighbor list considering only the prescribed number of neighbors
            [thisData, outDebug]...
                = neighPositions(data, neighbors, inDebug, reshapeData);

            %thisData = nNeighs by 3 matrix = [xpos ypos focalState]
            neighborAngles = atan2(thisData(:,2),thisData(:,1));
            %d = 7; % 7 cm max radius for angles -- set above
            idx_closeNeigh = ( vecnorm(thisData(:,1:2),2,2) < d );

            [totalLocs, density, vidStats] = collectiveData(data,area);

            allData{m,1} = thisData(:,1:2);
            allData{m,2} = neighborAngles( idx_closeNeigh );
            % proportion of side neighbors in this clip
            allData{m,3} = countNei( allData{m,2} );
            allData{m,4} = vidStats(2);


            idx_stop = (thisData(:,3) == 0);
            stopData{m,1} = thisData(idx_stop,1:2);
            stopData{m,2} = neighborAngles(idx_stop & idx_closeNeigh);
            % proportion of side neighbors in this clip
            stopData{m,3} = countNei( stopData{m,2} );

            idx_crawl = (thisData(:,3) == 1);
            crawlData{m,1} = thisData(idx_crawl,1:2);
            crawlData{m,2} = neighborAngles(idx_crawl & idx_closeNeigh);
            % proportion of side neighbors in this clip
            crawlData{m,3} = countNei( crawlData{m,2} );

            idx_hop = (thisData(:,3) == 2);
            hopData{m,1} = thisData(idx_hop,1:2);
            hopData{m,2} = neighborAngles(idx_hop & idx_closeNeigh);
            % proportion of side neighbors in this clip
            hopData{m,3} = countNei( hopData{m,2} );

            totalData(clipNum,:) = { allData{m,:}, stopData{m,:}, crawlData{m,:}, hopData{m,:} };
            clipNum = clipNum + 1;

        elseif strcmp(numNeighbors,'all') %%% THIS IS NOT WORKING CURRENTLY
            %this should also work, to simply use the giant neighbor cloud
            allData{m} = data_all{m,3}; 
            totalData{clipNum,1} = data_all{m,3};
            clipNum = clipNum + 1;
        end

        fprintf('Assembling that data took %f seconds \n', toc)
    end

end

end

function quad_nei = countNei(data)

quad_nei.side = sum( abs(data) > 3/4*pi | abs(data) < pi/4 );
quad_nei.all = numel(data);
quad_nei.nan = sum( isnan(data) );
quad_nei.prop = quad_nei.side / quad_nei.all;

end