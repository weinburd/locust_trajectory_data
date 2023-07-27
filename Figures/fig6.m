%plots relative neighbor density as appears in Figure 4 and Figure 1 (right).

close all
clear all % ensure we're loading new data

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
assembleData = 1; %and save it to figDataFile
    %else loads data from figDataFile

% the max number of neighbors around each focal individual
% can be an integer or 'all', but 'all' cannot distinguish focalState
numNeighbors = 100; % 100 gets all of them % 'all';
d = 7; % 7 cm max radius for angles

% for RESHAPING Data
% suppose a locust is approximatly 5mm wide and 15mm long
reshapeData = {'none', [0 0]}; 
reshapeData = {'rescale', [.5 1.5]}; 
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
        
        allData{m,1} = thisData(:,1:2);
        allData{m,2} = neighborAngles( idx_closeNeigh );
        
        idx_stop = (thisData(:,3) == 0);
        stopData{m,1} = thisData(idx_stop,1:2);
        stopData{m,2} = neighborAngles(idx_stop & idx_closeNeigh);
        
        idx_crawl = (thisData(:,3) == 1);
        crawlData{m,1} = thisData(idx_crawl,1:2);
        crawlData{m,2} = neighborAngles(idx_crawl & idx_closeNeigh);
        
        idx_hop = (thisData(:,3) == 2);
        hopData{m,1} = thisData(idx_hop,1:2);
        hopData{m,2} = neighborAngles(idx_hop & idx_closeNeigh);
        
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

% convert to a matrix array
allData = cell2mat(totalData(:,1));
stopData = cell2mat(totalData(:,1+vmc));
crawlData = cell2mat(totalData(:,1+2*vmc));
hopData = cell2mat(totalData(:,1+3*vmc));

% histogram options
dx = 0.5;
dy = 0.5;

binSizes = [dx dy];

%%% Histogram ALL %%%
plotrad = 14;
[binN.all, binCtrs.all] = RelNeiBins(allData, binSizes, plotrad);

plotrad = 7;
%%% Histogram STOP %%%
[binN.stop, binCtrs.stop] = RelNeiBins(stopData, binSizes, plotrad);

%%% Histogram CRAWL %%%
[binN.crawl, binCtrs.crawl] = RelNeiBins(crawlData, binSizes, plotrad);

%%% Histogram HOP %%%
[binN.hop, binCtrs.hop] = RelNeiBins(hopData, binSizes, plotrad);


save(figDataFile,   'binN', 'binCtrs',...
                    'binSizes','-append')

else
    load(figDataFile,   'binN', 'binCtrs',...
                        'binSizes')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% histogram options
% binSizes = loaded from the file

%%% FIGURE ALL %%%
plotrad = 14;
figTitle = sprintf("Relative Neighbor Density");
h_all = plotRelNeiDen(number, binN.all, binCtrs.all, binSizes, plotrad, figTitle, reshapeData);

plotrad = 7;
%binSizes = [0.25 0.25];
%%% FIGURE STOP %%%
figTitle = sprintf("Around Stationary Locusts");
h_stationary = plotRelNeiDen(number+1, binN.stop, binCtrs.stop, binSizes, plotrad, figTitle, reshapeData);

%%% FIGURE CRAWL %%%
figTitle = sprintf("Around Walking Locusts");
h_walk = plotRelNeiDen(number+2, binN.crawl, binCtrs.crawl, binSizes, plotrad, figTitle, reshapeData);

%%% FIGURE HOP %%%
figTitle = sprintf("Around Hopping Locusts");
h_hop = plotRelNeiDen(number+3, binN.hop, binCtrs.hop, binSizes, plotrad, figTitle, reshapeData);

fprintf('Plotting all the figures took %f seconds \n', toc)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Figure, Table, LaTeX %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFigs
    
    figPath = 'figs/';
    
    filename = sprintf('fig6_density_%s',reshapeData{1});
    print(h_all,[figPath filename, '.eps'],'-depsc')
    
    filename = sprintf('fig6_stationary_%s',reshapeData{1});
    print(h_stationary,[figPath filename, '.eps'],'-depsc')
    
    filename = sprintf('fig6_walk_%s',reshapeData{1});
    print(h_walk,[figPath filename, '.eps'],'-depsc')
    
    filename = sprintf('fig6_hop_%s',reshapeData{1});
    print(h_hop,[figPath filename, '.eps'],'-depsc')
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Associated Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Moved to separate function file neighPositions.m
% function[thisData, varargout] = neighPositions(data, neighbors, varargin)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [binN, binCtrs] = RelNeiBins(neighbors, binSizes, plotrad)

% histogram options
dx = binSizes(1);
dy = binSizes(2);
%plotrad = 14;
xedges = -plotrad:dx:plotrad;
yedges = -plotrad:dy:plotrad;
edges = {xedges, yedges};
% histogram counts
[N, ctrs] = hist3(neighbors,'Edges',edges,'CdataMode','auto','Normalization','pdf');
binN = N;
binCtrs = ctrs;

end

function h = plotRelNeiDen(figNum, binN, binCtrs, binSizes, plotrad, figTitle,varargin)

reshapeData = varargin{1};

factor = 2; % choose 1 or 2

wid = 3.75*factor;
hei = 3*factor;
h = figure(figNum);
set(h,'Units','Inches');
pos = get(h, 'Position');
set(h,'PaperPositionMode','Manual') % Setting this to 'manual' unlinks 'Position' (what gets displayed) and 'PaperPosition' (what gets printed)
set(h,'PaperPosition',[ 0 0 wid hei]);
set(h,'Position',[ 0 0 wid hei]);

% histogram options
dx = binSizes(1);
dy = binSizes(2);
%plotrad = 14;
xedges = -plotrad:dx:plotrad;
yedges = -plotrad:dy:plotrad;
edges = {xedges, yedges};
% histogram counts
N = binN;
ctrs = binCtrs;
%[N, ctrs] = hist3(neighbors,'Edges',edges,'CdataMode','auto','Normalization','pdf');

% convert to "relative density", according to Buhl et al. 2012

% When computing relative density, we need to account for 
% 1) area of circle or square, 
% 2) total number of neighbors in the plot
n = sum(N,'all');
if plotrad == 14
    N = N*1/(dx*dy)*pi*plotrad^2/n;
    focalSize = 3*factor;
    focalLine = [1/factor 2]; %length width
elseif plotrad == 7
    N = N*1/(dx*dy)*(2*plotrad)^2/n;
    focalSize = 8*factor;
    focalLine = [.75 2*factor];
end

xctrs = ctrs{1}; yctrs = ctrs{2}; % these are the bin centers
xctrs = xctrs(1:end-1); yctrs = yctrs(1:end-1); %cut off last empty bin
N = N(1:end-1,1:end-1);
% plot as a surface
surf(xctrs,yctrs,N', 'EdgeColor','None','FaceColor','interp')
colormap('turbo') % also try 'jet' and default 'parula'
xlim([xctrs(1) xctrs(end)])
ylim([yctrs(1) yctrs(end)])
daspect([1 1 1]) %pbaspect([1 1 1])
if plotrad == 7
    %caxis([0 1.6])
end
view(2)

%%% Add focal locust %%%
if ~(strcmp(reshapeData{1},'subtract') || strcmp(reshapeData{1},'subtractVert'))
hold on
quiver3(0,0,1.1*max(max(N)),0,focalLine(1),0,...
        'o','MarkerSize',focalSize,'LineWidth',focalLine(2),...
        'MarkerFaceColor','w','Color','w','AutoScale','off')
hold off
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Labels %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = gca;
axpos = ax.Position;
ax.PositionConstraint = 'innerposition';

figSubtitle = sprintf('With Body-Shape Correction: %s', reshapeData{1});

set(ax, 'FontSize', ticklabelsize)
xlabel('$\Delta x$ (cm)','FontSize',axislabelsize)
ylabel('$\Delta y$ (cm)','FontSize',axislabelsize)
title(figTitle,'FontSize',titlesize)
subtitle(figSubtitle,'FontSize',axislabelsize)

% colorbar labels
if strcmp(figTitle,'Relative Neighbor Density') || strcmp(figTitle,'Around Hopping Locusts')
    cbar = 1;
else
    cbar = 0;
    cwidth = 0;
end

if factor == 2
    cLabelAdjust = 1.3;
    cPosAdjust = [ 0 0 0 0 ];
elseif factor == 1
    cLabelAdjust = 1.25;
    cPosAdjust = [0.07 0 -.02 0];
end

if cbar
    c = colorbar(ax);
    cbpos = c.Position;
    c.Position = cbpos + cPosAdjust;
    c.Label.String = 'Relative Density';
    c.Label.FontSize = fontSize;
    c.Label.Rotation = 270;
    c.Label.Position(1) = c.Label.Position(1) + cLabelAdjust;
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    cwidth = c.Position(3)*5*factor;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Touch Up %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %remove whitespace
    ax.Units = 'Inches';
    %ax = gca;
    outerpos = ax.OuterPosition;
    axpos = ax.Position;
    ti = ax.TightInset; 
    %colorbarX = c.Position(3);
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
    left =  ti(1); %axpos(1) +
    bottom = axpos(2);% + ti(2);
    ax_width = axpos(3);% - ti(1) - ti(3);
    ax_height = axpos(4);% - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    ax.Units = 'Inches';
    %paperPos = [0 0 wid wid*outerpos(4)/outerpos(3)];
    wid = ax.Position(3)*(1 + cwidth);% + ax.TightInset(1) + ax.TightInset(3));
    paperPos = [0 0 wid hei];
    set(h,'PaperPosition',paperPos);
    set(h,'Position',paperPos);
end
