%plots histograms of instantaneous speed as appear in Figure 3 and Figure 1 (left).

close all
% clear all % ensure we're loading new data
addpath('../Data/')

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

figDataFile = 'fig_data.mat';

%%% Options %%%
saveFigs = 1;
assembleData = 0; %and save it to figDataFile
    %else loads data from figDataFile
    
    printStats = 1;
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 10s %%%
number = 10;
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

nRecs = size(recording,2);

nClips = 0;
for file_idx = 1:4
    nClips = nClips + size(recording(file_idx).data, 1);
end
totalSpeeds = cell(nClips, 4);

clipNum = 1;

%file_idx = 1; %vid133
%file_idx = 2; %vid098
%file_idx = 3; %vid096
%file_idx = 4; %vid146

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

    for m = matNums

        disp(['Now assembling data from vid' data_all{m,1}]);

        tic
        %%Here we assemble the data
        data_struct = data_all{m,2};

        [data_final, neighbors] = struct2data(data_struct);

        % assemble the individual data for this data set
        speeds = indivData(data_final,idx_speed);

        dataSpeeds(m,:) = speeds;

        totalSpeeds(clipNum,:) = speeds;
        clipNum = clipNum + 1;

        fprintf('Assembling that data took %f seconds \n', toc)
    end

end

save(figDataFile, 'totalSpeeds', '-append')

else
    load(figDataFile, 'totalSpeeds')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

%%% FIGURE Speeds ALL %%%
figTitle = sprintf('Histogram of Instantaneous Speeds');
h_all = plotIndividual(number, totalSpeeds(:,1), 'speed', figTitle);

%%% FIGURE Speeds STOP %%%
figTitle = sprintf('Stationary Locusts');
h_stationary = plotIndividual(number+2, totalSpeeds(:,2), 'speed', figTitle);

%%% FIGURE Speeds CRAWL %%%
figTitle = sprintf('Walking Locusts');
h_walk = plotIndividual(number+4, totalSpeeds(:,3), 'speed', figTitle);

%%% FIGURE Speeds HOP %%%
figTitle = sprintf('Hopping Locusts');
h_hop = plotIndividual(number+6, totalSpeeds(:,4), 'speed', figTitle);

fprintf('Plotting all the figures took %f seconds \n', toc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Figure, Table, LaTeX %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFigs
    figPath = 'figs/';
    
    filename = 'fig1_speeds';
    print(h_all,[figPath filename, '.eps'],'-depsc')
    
    filename = 'fig3_stationary';
    print(h_stationary,[figPath filename, '.eps'],'-depsc')
    
    filename = 'fig3_walk';
    print(h_walk,[figPath filename, '.eps'],'-depsc')
    
    filename = 'fig3_hop';
    print(h_hop,[figPath filename, '.eps'],'-depsc')
    
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


function h = plotIndividual(figNum, allData, dataFlag, figTitle)

Data = cell2mat(allData);

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

if strcmp(dataFlag, 'speed')
    
    topedge = ceil(max(Data));
    edges = 0:0.5:topedge;
    histogram(Data,edges);
    set(gca,'YScale','log')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Labels %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(gca, 'FontSize', ticklabelsize)
    xlabel('Speed (cm/sec)','FontSize',axislabelsize)
    % ylabel('log(counts)','FontSize',axislabelsize)
    ylabel(sprintf('log(counts)'),'FontSize',axislabelsize)
    title(figTitle,'FontSize',titlesize)
    xlim([0, 32])
%     if strcmp(figTitle,'Histogram of Instantaneous Speeds')
        ymax = 1.25e+6;
%     else
%         ymax = 1e+6;
%     end
    ylim([10^3, ymax])
elseif strcmp(dataFlag,'angle')
    C = 72;
    edges = linspace(0,2*pi,C)-2*pi/(2*72); % bins centered around zero
    hst = polarhistogram(Data(:),edges); % polarhistogram() applies mod 2*pi
    hct = hst.BinCounts;
    hbe = hst.BinEdges;
    %alternative bin edges that don't include -pi = pi
    hb2 = hbe(1:end-1) + diff(hbe) / 2;
    
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    title(figTitle,'FontSize',titlesize)
end

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
    
end