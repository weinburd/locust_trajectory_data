%%% Create a movie from data_final %%%
%function mov_locusts(recording, file_idx, clipNum, vidFlag)

% Requires the directory '../Data/'

close all
addpath('../Functions',...
        '../Data/examples')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
% global  idx_x idx_y idx_flag...
%         idx_speed idx_theta...
%         idx_localSpeed idx_localStd idx_localMinMax...
%         idx_state...
%         numNeighbors numFocal

idx_x = 1;
idx_y = 2;
idx_theta = 5;
idx_state = 9;

load( 'data_recording_examples.mat', 'recording_examples' )
recording = recording_examples;
file_idx = 1;
clipNum = 7; %can change to 9

%%% Options %%%
saveFigs = 1;
% assembleData = 0; %and save it to figDataFile
    %else loads data from figDataFile
%%%
vidFlag = 1; %whether or not to load data from the video clip

% Unpack the data and metadata
vidName = recording(file_idx).data{clipNum,1};
clipName = char(extractBetween(vidName,'tracks_','.xml'));

data_struct = recording(file_idx).data{clipNum,2};
[data_final, ~] = struct2data(data_struct);

if vidFlag
    v = ['preprocessed_' clipName(1:end-4) '.avi']; %seem to need to remove the final "_720"
    % The video file is incompatible with MacOS
    % See an explanation here: 
    % https://www.mathworks.com/matlabcentral/answers/395560-why-does-videoreader-not-load-a-video-in-r2017a-or-later-that-used-to-load-in-r2016b-on-mac
    video = VideoReader(v);
end

if isempty(recording(file_idx).corners)
    scale = recording(file_idx).scale; % pi/cm
    trans = [];
    fieldDims = recording(file_idx).fieldDims / scale;
else
    cornersPix = recording(file_idx).corners;
    fieldDimsPix = recording(file_idx).fieldDims;
    [trans, scale, fieldDims, newR_A] = projTrans(cornersPix, fieldDimsPix); % in cm
end
wid = fieldDims(2)-fieldDims(1);
hei = fieldDims(4)-fieldDims(3);

% chose some size parameters
figWidth = 16;
titleSize = 20; %fontsize
labelSize = 24;
ticklabelSize = 18;
trail_width = 5; % sample trajectory widths
alpha = .4; % sample trajectory transparency

%% Set up a figure for plotting
figNum = 1011;
h = figure(figNum);
set(h,'Units','Inches');
pos = get(h, 'Position');
set(h,'PaperPositionMode','Manual') 
% Setting this to 'manual' unlinks 'Position' (what gets displayed) and 'PaperPosition' (what gets printed)
set(h,'PaperPosition',[ 0 0 figWidth figWidth/wid*hei]);
set(h,'Position',[ 0 0 figWidth figWidth/wid*hei]);
movegui(h,'north')


pbaspect([wid,hei,1]);
% colors
orange = [230, 159, 0]/255;
skyblue = [86, 180, 233]/255;
blueishgreen = [0, 158, 115]/255;
yellow = [240,228, 66]/255;
blue = [0, 114, 178]/255;
redpurple = [204, 121, 167]/255;
vermillion = [213, 94, 0]/255;
colors = lines(7);

%% Extract some locust data
[Nlocs, Nfeats, Ntimes] = size(data_final);

% rng(42)
% nii = 352; % Capture trajectories for 3 locusts
% loc_idx = randsample(352, nii)
% loc_filt = false(Nlocs, 1);
% loc_filt(loc_idx) = true;
% loc_filt(:) = true;

x_line = squeeze( data_final(:, idx_x, :) );
y_line = squeeze( data_final(:, idx_y, :) );
state_line = squeeze(data_final(:, idx_state, :) );
    idx_stop = state_line == 0;
    idx_crawl = state_line == 1;
    idx_hop = state_line == 2;
theta_line = squeeze( data_final(:, idx_theta, :) );

x_line = x_line(:,:)*scale;
y_line = fieldDims(4)*scale-y_line(:,:)*scale;

%% For finding good locusts
% % Try to idenity the (horizontally) longest trajectorys
% [r, c] = find(~isnan(x_line));
% idx = accumarray(r, c, [], @(x) {[min(x) max(x)]} );
% idx = vertcat(idx{:});
% 
% x_first = x_line( sub2ind(size(x_line), 1:nii, idx(:,1)') );
% x_last = x_line( sub2ind(size(x_line), 1:nii, idx(:,2)') );
% x_diff = x_first-x_last;
% 
% [M, I] = maxk(x_diff, 20);
% 
% I;

%% Set up video frame
T_end = 183;

    if vidFlag
        vidFrame = read(video,T_end);
        image(flipud(vidFrame))
        hold on
    end
    
%% With locusts chosen, plot lines
loc_idx = [7 92 106]
loc_filt = false(Nlocs, 1);
loc_filt(loc_idx) = true;

hold on

for ii = loc_idx
    stop_steps = find(idx_stop(ii,1:T_end));
    crawl_steps = find(idx_crawl(ii,1:T_end));
    hop_steps = find(idx_hop(ii,1:T_end));
    
    for jj = stop_steps
        step = jj-1:jj;
        plot(x_line(ii, step), y_line(ii, step), "color", [colors(7,:) alpha], 'LineWidth', trail_width )
    end
    for jj = crawl_steps
        step = jj-1:jj;
        plot(x_line(ii, step), y_line(ii, step), "color", [yellow alpha], 'LineWidth', trail_width)
    end
    for jj = hop_steps
        step = jj-1:jj;
        plot(x_line(ii, step), y_line(ii, step), "color", [blueishgreen alpha], 'LineWidth', trail_width)
    end
end

%myQuiv([true;false;true], blueishgreen, x_line(loc_filt,T_end), y_line(loc_filt,T_end), theta_line(loc_filt,T_end), scalefactor)
%myQuiv([false;true;false], yellow, x_line(loc_filt,T_end), y_line(loc_filt,T_end), theta_line(loc_filt,T_end), scalefactor)

% a non-existent scatter plot to create the legend.
% try as I might I could not make the little circles bigger!
q_legend = { scatter(nan, nan, 'MarkerEdgeColor', colors(7,:), "LineWidth", 1.5), ...
             scatter(nan, nan, 'MarkerEdgeColor', yellow, "LineWidth", 1.5), ...
             scatter(nan, nan, 'MarkerEdgeColor', blueishgreen, "LineWidth", 1.5)};

qs = myQuiv(idx_stop(:,T_end), colors(7,:), x_line(:,T_end), y_line(:,T_end), theta_line(:,T_end), scalefactor);
qc = myQuiv(idx_crawl(:,T_end), yellow, x_line(:,T_end), y_line(:,T_end), theta_line(:,T_end), scalefactor);
qh = myQuiv(idx_hop(:,T_end), blueishgreen, x_line(:,T_end), y_line(:,T_end), theta_line(:,T_end), scalefactor);

%% add a legend?
%TODO
l = legend([q_legend{:}], ...
       {"Stationary Locust", "Walking Locust", "Hopping Locust"},...
       "Location", "SouthEast", ...
       "FontSize", ticklabelSize+2);

%% Axes and stuff

scalefactor = sqrt(wid*hei)*0.02*scale;
imageAxis = [fieldDims(1:2) fieldDims(3) fieldDims(4)]*scale;
%[fieldDims(1:2) -fieldDims(4) -fieldDims(3)]*scale;
axis(imageAxis)
ax = gca;

ax.FontSize = ticklabelSize;
xtlab = linspace(0,100,11);
ytlab = linspace(0,50,6);
xticks(xtlab*scale)
yticks(ytlab*scale)
xticklabels(cellfun(@num2str,num2cell(xtlab),'un',0))
yticklabels(cellfun(@num2str,num2cell(ytlab),'un',0))

%title(['Time ',num2str(0,'%3.1f'), ' (sec)'],'FontSize',titleSize);
title('Single Video Frame with Tracked Positions and Sample Trajectories', 'FontSize',labelSize+2)
xlabel('$x$ (centimeters)','FontSize',ticklabelSize+2)
ylabel('$y$ (centimeters)','FontSize',ticklabelSize+2)

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
    ax_height = outerpos(4) - ti(2) - ti(4) - .015;
    ax.Position = [left bottom ax_width ax_height];


%% save figure and data
%add assemble data functionality?

if saveFigs
    figPath = 'figs/';
    
    filename = 'fig1_traj';
    print(h,[figPath filename, '.eps'],'-depsc')
end


%% Associated Functions %%

    function q = myQuiv(idx,color, xposn, yposn, theta, scalefactor)
        q = quiver(xposn(idx),yposn(idx),scalefactor*cos(theta(idx)),scalefactor*sin(theta(idx)),...
               'o','MarkerSize', 14,'Color',color,'AutoScale','off', 'LineWidth', 1.5);
        q.ShowArrowHead = 'off';
        q.MaxHeadSize = 0;
    end

%end