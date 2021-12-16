function recording = initData(file_idx,varargin)

nOpt = numel(varargin);
if nOpt > 0
    acc = varargin{1};
end

if (1 <= file_idx) && (file_idx <= 2)
    if nOpt>0
        dataFile = ['Data/examples/data_recording_' acc '.mat'];
    else
        dataFile = 'Data/data_recording.mat';
    end
    disp(['Now loading your data from ' dataFile]);
    tic
    load( dataFile )    
    fprintf(['That took %f seconds', newline],toc)
elseif (3 <= file_idx) && (file_idx <= 4)
    dataFile = 'Data/data_recording2.mat';
    disp(['Now loading your data from ' dataFile]);
    tic
    load( dataFile )
    recording = recording2;
    fprintf(['That took %f seconds', newline],toc)
end


if nOpt > 0
    recording = eval(['recording_' acc]);
end
%recording = struct();%comment out so as not to overwrite

% the name of the file
recording(1).file = 'VS00133noteLightVariation.MTS';
% the dimensions in pixels of the the video
recording(1).fieldDims = [0 1824 0 1026];
% either the scale factor or the corners of the board outside the video field of view 
recording(1).scale = 18.34113771; % pix/cm
% location of corners of the board after cropping, possibly outside the video frame
recording(1).corners = []; %cannot find them for this video
% all the data
if nOpt == 0
    names = { '133_1min_225',...
              '133_1min_650',...
              '133_1min_750',...
              '133_1min_850',...
              '133_1min_950',...
              '133_1min_1115'};
    if ~isfield(recording(1),'data')
        recording(1).data = cell(length(names),4);
    end
    recording(1).data(:,1) = names;
end

% the name of the file
recording(2).file = 'VS00098.MTS';
% the dimensions in pixels of the the video
recording(2).fieldDims = [0 1894 0 976];
% either the scale factor or the corners of the board outside the video field of view 
recording(2).scale = 17.26023729; % pix/cm
% location of corners of the board after cropping, possibly outside the video frame
recording(2).corners = []; %cannot find them for this video
% all the data
names = {'098_1min_430',...
         '098_1min_530',...
         '098_1min_630',...
         '098_1min_730',...
         '098_1min_5830',...
         '098_1min_5930',...
         '098_1min_10030'};
if ~isfield(recording(2),'data')
    recording(2).data = cell(length(names),4);
end
recording(2).data(:,1) = names;

% the name of the file
recording(3).file = '00096.MTS';
% the dimensions in pixels of the the video
recording(3).fieldDims = [0 1496 0 802];
% either the scale factor or the corners of the board outside the video field of view 
recording(3).scale = []; % pix/cm
% location of corners of the board after cropping, possibly outside the video frame
% [ top right; top left; bottom left; bottom right]
recording(3).corners = [1509.250  -28.125;  -145.188  -37.625;  -184.625  813.625;  1568.833  803.667]; 
% all the data
names = {   '096_1min_500',...
            '096_1min_600',...
            '096_1min_700',...
            '096_1min_800',...
            '096_1min_900',...
            '096_1min_1000',...
            '096_1min_1100'     };
if ~isfield(recording(3),'data')
    recording(3).data = cell(length(names),4);
end
recording(3).data(:,1) = names;

% the name of the file
recording(4).file = '00146_Low_density.MTS';
% the dimensions in pixels of the the video
recording(4).fieldDims = [0 1650 0 844];
% either the scale factor or the corners of the board outside the video field of view 
recording(4).scale = []; % pix/cm
% location of corners of the board after cropping, possibly outside the video frame
% [ top right; top left; bottom left; bottom right]
recording(4).corners = [1651.333  -18.083;  -165.417  -3.333;  -220.0093079  924.3918079;  1717.392111  927.8515579]; 
% all the data
names = {   '146_1min_515',...
            '146_1min_615',...
            '146_1min_715',...
            '146_1min_815',...
            '146_1min_915',...
            '146_1min_1015',...
            '146_1min_1115'     };
if ~isfield(recording(4),'data')
    recording(4).data = cell(length(names),4);
end
recording(4).data(:,1) = names;


end