close all
% clear all % ensure we're loading new data

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

figDataFile = 'fig_data.mat';

%%% Options %%%
saveFigs = 1;
assembleData = 0; %and save it to figDataFile
    %else loads data from figDataFile

% the max number of neighbors around each focal individual
% can be an integer or 'all', but 'all' cannot distinguish focalState
numNeighbors = 100; % 100 gets all of them % 'all';
d = 7; % 7 cm max radius for angles
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
        % extract position and orientation information and motion state
        data = [data_final(:,idx_x:idx_y,:) data_final(:,idx_theta,:) data_final(:,idx_state,:)];

        %%% debugging... %%%
        count = 0;
        countnan = 0;
        countempty = 0;
        countnoneighs = 0;
    
        % assemble a new neighbor list considering only the prescribed number of neighbors
        [thisData, count, countnan, countempty, countnoneighs]...
            = neighPositions(data, neighbors, count, countnan, countempty, countnoneighs);
        
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
[binN_all, binCtrs_all] = RelNeiBins(allData, binSizes, plotrad);

plotrad = 7;
%%% Histogram STOP %%%
[binN_stop, binCtrs_stop] = RelNeiBins(stopData, binSizes, plotrad);

%%% Histogram CRAWL %%%
[binN_crawl, binCtrs_crawl] = RelNeiBins(crawlData, binSizes, plotrad);

%%% Histogram HOP %%%
[binN_hop, binCtrs_hop] = RelNeiBins(hopData, binSizes, plotrad);

save(figDataFile,   'binN_all', 'binCtrs_all',...
                    'binN_stop', 'binCtrs_stop',...
                    'binN_crawl', 'binCtrs_crawl',...
                    'binN_hop', 'binCtrs_hop', '-append')

else
    load(figDataFile,   'binN_all', 'binCtrs_all',...
                        'binN_stop', 'binCtrs_stop',...
                        'binN_crawl', 'binCtrs_crawl',...
                        'binN_hop', 'binCtrs_hop')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% histogram options
dx = 0.5;
dy = 0.5;

binSizes = [dx dy];

%%% FIGURE ALL %%%
plotrad = 14;
figTitle = sprintf("Relative Neighbor Density");
h_all = plotRelNeiDen(number, binN_all, binCtrs_all, binSizes, plotrad, figTitle);

plotrad = 7;
%binSizes = [0.25 0.25];
%%% FIGURE STOP %%%
figTitle = sprintf("Around Stationary Locusts");
h_stationary = plotRelNeiDen(number+1, binN_stop, binCtrs_stop, binSizes, plotrad, figTitle);

%%% FIGURE CRAWL %%%
figTitle = sprintf("Around Walking Locusts");
h_walk = plotRelNeiDen(number+2, binN_crawl, binCtrs_crawl, binSizes, plotrad, figTitle);

%%% FIGURE HOP %%%
figTitle = sprintf("Around Hopping Locusts");
h_hop = plotRelNeiDen(number+3, binN_hop, binCtrs_hop, binSizes, plotrad, figTitle);

fprintf('Plotting all the figures took %f seconds \n', toc)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Figure, Table, LaTeX %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFigs
    
    figPath = 'figs/';
    
    filename = 'fig1_density';
    print(h_all,[figPath filename, '.eps'],'-depsc')
    
    filename = 'fig4_stationary';
    print(h_stationary,[figPath filename, '.eps'],'-depsc')
    
    filename = 'fig4_walk';
    print(h_walk,[figPath filename, '.eps'],'-depsc')
    
    filename = 'fig4_hop';
    print(h_hop,[figPath filename, '.eps'],'-depsc')
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Associated Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [thisData, varargout] = neighPositions(data, neighbors, varargin)

global numNeighbors numFocal

count = varargin{1};
countnan = varargin{2};
countempty = varargin{3};
countnoneighs = varargin{4};

Ntimesteps = size(neighbors,2);
thisData = cell(Ntimesteps,1);

    for t = 1:Ntimesteps
        % find the nonempty entries for this frame, i.e. the indices of locusts that have neighbors
        locusts = find( ~cellfun( @isempty, neighbors(:,t) ) );
        nowNeighbors = cell(size(locusts,1),1);
        for locust = locusts' % locusts must be a ROW vector
            neigh_idx = neighbors{locust,t};
            if isnan(neigh_idx) % focal locust is either not in the frame or is in an edge box
                % This never happens anymore because I took out NaNs for locusts out of the frame or on the edge.
                %allData = cat(1,allData,[NaN, NaN]); %unnesseccary
                countnoneighs = countnoneighs + 1;
            else
                numFocal = numFocal + 1;
                rel_posn = data(neigh_idx,1:2,t)-data(locust,1:2,t);
        
                % convert to complex numbers
                % y-component gets a negative sign due to ImageJ pixel coordinate
                crel_posn = rel_posn(:,1)-1i*rel_posn(:,2); 
                
                if isa(numNeighbors,'double')
                    % consider only the prescribed number of neighbors
                    [~,I] = mink(abs(crel_posn), numNeighbors);
                    crel_posn = crel_posn(I);
                    last_idx = min(numNeighbors,length(crel_posn));
                elseif isa(numNeighbors,'char')
                    % consider all neighbors
                    last_idx = length(crel_posn);
                end
                
                % rotate by focal locust's orientation data(locust,3,t)
                % and by a factor of pi/2 so that focal locust is facing up
                crel_posn = abs(crel_posn)...
                            .*exp(1i*angle(crel_posn) - 1i*data(locust,3,t) + 1i*pi/2);
                        
                %debugging
                %global THETAS
                %THETAS = [THETAS data(locust,3,t)];
                
                % convert back to (x,y) positions
                rel_posn = [real(crel_posn) imag(crel_posn)];
                
                %%% debugging... %%%
                if isnan(rel_posn(1:last_idx,:)) %(data(locust,3,t)) %using data catches some places where the heading is undefined but there are no neighbors
                % focal locust was in the frame, but heading = NaN because it either just entered or just left the frame
                    countnan = countnan + size(rel_posn(1:last_idx,:),1);
                elseif sum(sum(isnan(rel_posn(1:last_idx,:))))>0 %if there are ANY NaN values
                    error('There is a NaN hiding in rel_posn, but it is not all NaNs')
                end
                if isempty(rel_posn(1:last_idx,:)) %(neigh_idx)
                % focal locust was in the frame and not in an edge box, but has no neighbor within min(dx,dy)
                    countempty = countempty + 1;
                end
                count = count + max(1,size(rel_posn(1:last_idx,:),1));
                % should have: count == countnan+countempty+size(rmmissing(allData),1)
                %        also: count == countempty+size(allData,1)
                
                theseNeighbors = rel_posn(1:last_idx,:);
            end
            nNeighs = size(theseNeighbors,1);
            nowNeighbors{locusts==locust} = [ theseNeighbors data(locust,4,t)*ones(nNeighs,1) ];
        end
        thisData{t} = cell2mat(nowNeighbors);
    end

thisData = cell2mat(thisData);

varargout{1} = count;
varargout{2} = countnan;
varargout{3} = countempty;
varargout{4} = countnoneighs;

end

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

function h = plotRelNeiDen(figNum, binN, binCtrs, binSizes, plotrad, figTitle)

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
    caxis([0 1.6])
end
view(2)

%%% Add focal locust %%%
hold on
quiver3(0,0,1.1*max(max(N)),0,focalLine(1),0,...
        'o','MarkerSize',focalSize,'LineWidth',focalLine(2),...
        'MarkerFaceColor','w','Color','w','AutoScale','off')
hold off

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

set(ax, 'FontSize', ticklabelsize)
xlabel('$\Delta x$ (cm)','FontSize',axislabelsize)
ylabel('$\Delta y$ (cm)','FontSize',axislabelsize)
title(figTitle,'FontSize',titlesize)

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
