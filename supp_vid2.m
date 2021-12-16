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
vidFlag = 1;

% Unpack the data and metadata
vidName = recording(file_idx).data{clipNum,1};
clipName = char(extractBetween(vidName,'tracks_','.xml'));

data_struct = recording(file_idx).data{clipNum,2};
[data_final, ~] = struct2data(data_struct);

if vidFlag
    v = ['preprocessed_' clipName '.avi'];
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


% Create a movie writer 'Figures/' OR 'C:\Users\jazzy\Desktop\'
vid = VideoWriter(clipName,'MPEG-4');
vid.Quality = 100;
%vid = VideoWriter(['Figures/' clipName],'Uncompressed AVI');
vid.FrameRate = 25;
open(vid);

% chose some size parameters
figWidth = 16;
titleSize = 20; %fontsize
labelSize = 14;
ticklabelSize = 12;

% Set up a figure for plotting
figNum = 101;
h = figure(figNum);
set(h,'Units','Inches');
pos = get(h, 'Position');
set(h,'PaperPositionMode','Manual') 
% Setting this to 'manual' unlinks 'Position' (what gets displayed) and 'PaperPosition' (what gets printed)
set(h,'PaperPosition',[ 0 0 figWidth figWidth/wid*hei]);
set(h,'Position',[ 0 0 figWidth figWidth/wid*hei]);
movegui(h,'north')

title(['Time ',num2str(0,'%3.1f'), ' (sec)'],'FontSize',titleSize);
xlabel('centimeters','FontSize',labelSize)
ylabel('centimeters','FontSize',labelSize)

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

% Time loop through the data
[Nlocs, Nfeats, Ntimes] = size(data_final);


for t=1:Ntimes
    set(gca,'nextplot','replacechildren');

    if vidFlag
        vidFrame = read(video,t);
        image(flipud(vidFrame))
        hold on
    end

    xposn = data_final(:,idx_x,t);
    yposn = data_final(:,idx_y,t);

    idx = ~isnan(xposn);
    
    state = data_final(idx,idx_state,t);
    idx_stop = state == 0;
    idx_crawl = state == 1;
    idx_hop = state == 2;
    theta = data_final(idx,idx_theta,t);
    xposn = xposn(idx)*scale;
    yposn = fieldDims(4)*scale-yposn(idx)*scale;

    myQuiv(idx_stop,colors(7,:))
    hold on;
    myQuiv(idx_crawl,yellow)
    myQuiv(idx_hop,blueishgreen)
    %axis([0-0.02*Lx 1.02*Lx 0-0.02*Ly 1.02*Ly]); 
    title(['Time ',num2str(t*0.04,'%3.1f'), ' (sec)'],'FontSize',titleSize);
    drawnow;
    box on


    %%% Write this frame to video
    frame = getframe(figNum);
    writeVideo(vid,frame);
    
    hold off;
end

% save video at the end
close(vid)

%% Associated Functions %%

    function myQuiv(idx,color)
        quiver(xposn(idx),yposn(idx),scalefactor*cos(theta(idx)),scalefactor*sin(theta(idx)),...
               'o','MarkerSize', 14,'Color',color,'AutoScale','off', 'LineWidth', 1.5);
    end

%end