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

rMax = 14;
%%% Radius vs. Anisotropy -- STOP %%%
Data = stopData;
aniData = Data( vecnorm( Data, 2, 2 ) < rMax, : );
[Ms_stop, binCs_stop] = computeAnisotropy(aniData, rMax);

%%% Radius vs. Anisotropy -- CRAWL %%%
Data = crawlData;
aniData = Data( vecnorm( Data, 2, 2 ) < rMax, : );
[Ms_crawl, binCs_crawl] = computeAnisotropy(aniData, rMax);

%%% Radius vs. Anisotropy -- HOP %%%
Data = hopData;
aniData = Data( vecnorm( Data, 2, 2 ) < rMax, : );
[Ms_hop, binCs_hop] = computeAnisotropy(aniData, rMax);

save(figDataFile,   'Ms_stop', 'binCs_stop',...
                    'Ms_crawl', 'binCs_crawl',...
                    'Ms_hop', 'binCs_hop', '-append') %

else
    load(figDataFile,   'Ms_stop', 'binCs_stop',...
                        'Ms_crawl', 'binCs_crawl',...
                        'Ms_hop', 'binCs_hop')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% histogram options
dx = 0.5;
dy = 0.5;

binSizes = [dx dy];

rMax = 14;
%%% FIGURE Radius vs. Anisotropy -- STOP %%%
figTitle = sprintf("Around Stationary Locusts");
h_stationary = plotAnisotropy(number + 8, Ms_stop, binCs_stop, figTitle);

%%% FIGURE Radius vs. Anisotropy -- CRAWL %%%
figTitle = sprintf("Around Walking Locusts");
h_walk = plotAnisotropy(number + 9, Ms_crawl, binCs_crawl, figTitle);

%%% FIGURE Radius vs. Anisotropy -- HOP %%%
figTitle = sprintf("Around Hopping Locusts");
h_hop = plotAnisotropy(number + 10, Ms_hop, binCs_hop, figTitle);

fprintf('Plotting all the figures took %f seconds \n', toc)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Figure, Table, LaTeX %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFigs
    
    figPath = 'figs/';
    
    filename = 'fig5_stationary';
    print(h_stationary,[figPath filename, '.eps'],'-depsc')
    
    filename = 'fig5_walk';
    print(h_walk,[figPath filename, '.eps'],'-depsc')
    
    filename = 'fig5_hop';
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
function [MsCell, binCtrs] = computeAnisotropy(Data, rMax)

disp('Now computing the anisotropy by radius...')
tic;

Data = rmmissing(Data);

dists = vecnorm(Data(:,1:2),2,2);
angles = atan2( Data(:,2), Data(:,1));
NAngles = size(angles,1);

% we want overlapping annuli, each of width 1/2 and spaced by a distance of 1/4
binSize = 1;
% to create overlapping bins, we use histcounts on finer bins first
% binsP stands for binsPrime, the smaller/finer bins
% binsPinBins is the number of smaller/finer bins in one of my big bins
binsPinBins = 4;
binPsize = binSize/binsPinBins; 
nBinsP = round(rMax/binPsize);

[Np, edgesP, binP_idx] = histcounts(dists, nBinsP);

nBins = nBinsP-(binSize/binPsize-1);
binCtrs = linspace((binsPinBins/2)*binPsize,rMax-(binsPinBins/2)*binPsize,nBins);

% Initialize variables for looping
% Ms = [ M^c_1 -M^s_1 M^c_2 -M^s_2 ]; for each bin
Ms = zeros(nBins,4);
Mstds = zeros(nBins,4);
phis = zeros(nBins,2);
binNAngle = zeros(nBins,1);
% and for the uniform distribution, 
% it's same for all modes and both components
MsUni = zeros(nBins,1);
uniStds = zeros(nBins,1);


for bin = 1:nBins
    bin_idx = (binP_idx == bin);
    for numBins = 1:binsPinBins-1
        bin_idx = bin_idx | (binP_idx == bin+numBins); %union the other bins in
    end
    
    thisBinAngles = angles(bin_idx);
    thisBinN = size(thisBinAngles,1);
    binNAngle(bin) = thisBinN;
    
    n = 25;
    bootMs = bootstrp(n, @(samples) [   circ_moment(samples,[],1),...
                                        circ_moment(samples,[],2)...
                                     ], thisBinAngles);
    % scale the Ms by:
    %   density in annulus
    %   1/(density in disc with rad = 14)
    AofAnn = pi*(binCtrs(bin)+binSize).^2-pi*(binCtrs(bin)-binSize).^2;
    AofDisc = pi*14^2;
    fact = (thisBinN/AofAnn)/(NAngles/AofDisc);

    bootMs = bootMs*fact;

    % computed analytically, which is 0
    %MsUni = 0; set above
    % and one standard deviation from it's mean,
    uniStds(bin) = 1/sqrt(2*thisBinN) + 1i/sqrt(2*thisBinN);
    % since both components are the same, no need for both
    uniStds(bin) = 1/sqrt(2*thisBinN);
    uniStds(bin) = uniStds(bin)*fact;
        
    %%% Standard Trig Moments
    % Ms = [ M^c_1 -M^s_1 M^c_2 -M^s_2 ]; for each bin
    RIbootMs = [real(bootMs(:,1)) -imag(bootMs(:,1)) real(bootMs(:,2)) -imag(bootMs(:,2)) ];
    meanM = mean( RIbootMs, 1 );
    stdM = std( RIbootMs, 1);

    bootPhis = angle( bootMs );
    meanPhi = circ_mean(bootPhis,[],1);
    stdPhi = circ_std(bootPhis,[],[],1);
    
    Ms(bin,:) = meanM;
    Mstds(bin,:) = stdM;
    phis(bin,:) = meanPhi;
end

MsCell = {Ms, Mstds, uniStds};

fprintf('That took %f seconds \n', toc)
end

function h = plotAnisotropy(figNum, MsCell, binCtrs, figTitle)

Ms = MsCell{1};
Mstds = MsCell{2};
uniStds = MsCell{3};

disp('Now making the anisotropy plot...')
tic;

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

nBins = length(binCtrs);
MsUni = zeros(nBins,1);

% for confidence intervals
alph = 0.01; % this gives us the 100(1-alph)% = 99% CI
% tstar = tinv(1-alph,n-1);
% zstar = icdf('normal',1-alph/2,0,1);

% Ms = [ M^c_1 -M^s_1 M^c_2 -M^s_2 ]; for each bin
p = plot( binCtrs, Ms(:,2),'-.', binCtrs, Ms(:,3),'-', 'LineWidth', lwidth );
uistack(p(1),'top')
hold on
for mode = [2 3]
    y = Ms(:,mode);
    x = binCtrs;
    std_dev = Mstds(:,mode);
    curve1 = y + 5*std_dev; % tstar*std_dev/sqrt(n); %
    curve2 = y - 5*std_dev; % tstar*std_dev/sqrt(n);
    x2 = [x, fliplr(x)];
    inBetween = [curve1', fliplr(curve2')];
    f(mode) = fill(x2, inBetween, [.8 .8 .8], 'FaceAlpha', 0.5, 'LineStyle', ':');
    uistack(f(mode),'bottom')
end

% now for the unidistribution
y = MsUni;
x = binCtrs;
std_dev = uniStds(:,1);
curve1 = y + 5*std_dev; %zstar*std_dev/sqrt(n);
curve2 = y - 5*std_dev; %zstar*std_dev/sqrt(n);
x2 = [x, fliplr(x)];
inBetween = [curve1', fliplr(curve2')];
pUni = plot( x, y, 'LineWidth', lwidth, 'Color', 'k');
fUni = fill(x2, inBetween, [.8 .8 .8], 'FaceAlpha', 0.5, 'LineStyle', ':');
uistack(pUni, 'bottom')
uistack(fUni,'bottom')

hold off

% figure(figNum+1)
% plot( binCtrs, phis(:,1) )

ax = gca;

% also uncomment this line, below:
% axes(ax2) % brings it to the front 

insetFlag = 0;
if insetFlag
    % NEED TO FIX POSITION OF LEGEND IF insetFlag = true
    ax2 = axes('Position',[.5 .6 .4 .3]*factor);
    box on

    % 
    % %%% Max Anisotropy Angle Distributions
    % for mode = 1:2
    %     [M, maxBin] = max(Ms(:,mode));
    %     %maxBin = find( binCtrs==2.25 );
    %     bin_idx = (binP_idx == maxBin);
    %     for numBins = 1:binsPinBins-1
    %         bin_idx = bin_idx | (binP_idx == maxBin+numBins); %union the other bins in
    %     end
    %     
    %     maxBinAngles = angles(bin_idx);
    %     maxBinN = size(thisBinAngles,1);
    %     
    %     C = 36; % 2*pi/(pi/6);
    %     [counts,edges] = histcounts(maxBinAngles,linspace(-pi,pi,C)-pi/(2*C));
    % 
    %     idx_cut = round(C*1/4);
    %     first = 1:idx_cut;
    %     second = idx_cut+1:C-1;
    %     dphi = edges(2)-edges(1);
    %     phis = edges(1:end-1);
    %     phis = [phis(second) phis(first)+2*pi]-dphi/2;
    %     plot(phis,[counts(second) counts(first)]/dphi)
    %     hold on
    % end
    % hold off
    %%% Perpendicular anisotropy measures
    % Ms = [ M^c_1 -M^s_1 M^c_2 -M^s_2 ]; for each bin
    pmini = plot( binCtrs, Ms(:,1), binCtrs, Ms(:,4), 'LineWidth', lwidth );
    uistack(pmini(1),'top')
    hold on
    for mode = [1 4]
        y = Ms(:,mode);
        x = binCtrs;
        std_dev = Mstds(:,mode);
        curve1 = y + 5*std_dev; % tstar*std_dev/sqrt(n); %
        curve2 = y - 5*std_dev; % tstar*std_dev/sqrt(n);
        x2 = [x, fliplr(x)];
        inBetween = [curve1', fliplr(curve2')];
        fmini(mode) = fill(x2, inBetween, [.8 .8 .8], 'FaceAlpha', 0.5, 'LineStyle', ':');
        uistack(fmini(mode),'bottom')
    end
    
    % now for the unidistribution
    y = MsUni;
    x = binCtrs;
    std_dev = uniStds(:,1);
    curve1 = y + 5*std_dev; %zstar*std_dev/sqrt(n);
    curve2 = y - 5*std_dev; %zstar*std_dev/sqrt(n);
    x2 = [x, fliplr(x)];
    inBetween = [curve1', fliplr(curve2')];
    pUniMini = plot( x, y, 'LineWidth', lwidth, 'Color', 'k');
    fUniMini = fill(x2, inBetween, [.8 .8 .8], 'FaceAlpha', 0.5, 'LineStyle', ':');
    uistack(pUniMini, 'bottom')
    uistack(fUniMini,'bottom')
    
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Labels %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ax)
set(ax, 'FontSize', ticklabelsize)
xlabel('Distance from Focal (cm)','FontSize',axislabelsize)
ylabel('Fraction of Anisotropy','FontSize',axislabelsize)
title(figTitle,'FontSize',titlesize)
%ylim([-6*10^(-5), 16*10^-5]) % when using an old scaling...
ylimits = [-.03 0.11];
ylim(ylimits)
if strcmp(figTitle,'Around Hopping Locusts')
    legend([p(1) p(2) pUni f(2)],...
            {'front-back asymmetry', 'four-fold anisotropy',...
            'anisotropy of uniform dist.', '$\pm 5$ standard deviations'},...
            'Location','northeast')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Touch Up %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %remove whitespace
    %ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];

if insetFlag
    axes(ax2) % brings it to the front
    ylim(ylimits)
end

fprintf('That took %f seconds \n', toc)
end