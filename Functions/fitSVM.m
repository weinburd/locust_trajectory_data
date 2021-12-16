function model = fitSVM(varargin)
%%% Fit an SVM to assign motion state %%%
%%% Depends on data_accuracy.mat
%%% Input:
% M = an Nlocust by 9 by Ntimestep matrix
%       each locust has 9 features  at each timestep
%       the features are [x y flag theta speed MagAveVelocity localStdSpeed min(fwdMax,bkwdMax) state]
%       the x,y spatial units are centimeters
%       flag specifies whether or not this data was interpolated
%       columns 4 to 8 are inferred from position
%           theta = heading direction (radians) computed from velocity?
%           speed (cm/sec) is computed as forward/backward difference? with/without smoothing of positions depending on smoothFlag
%           MagAveVelocity = average is taken over a centered window of 15 timesteps
%           localStdSpeed = standard deviation of the speeds in a centered window of 15 timesteps
%           min(fwdMaxSpd,bkwdMaxSpd): 
%               fwdMaxSpd = maximum speed in the future 15 timesteps
%               bkwdMaxSpd = maximum speed in the past 15 timesteps
%
%%% Output:
% model = a fitcecoc() model that can predict the motion state of locusts as
%           NaN = locust not in view or only in view for one frame (first timestep is all NaNs)
%           0 = stationary
%           1 = crawling
%           2 = hopping


global  idx_x idx_y idx_flag...
        idx_speed idx_theta...
        idx_localSpeed idx_localStd idx_localMinMax...
        idx_state

% idx_x = 1; idx_y = 2; idx_flag = 3;
idx_speed = 4; %idx_theta = 5;
idx_localSpeed = 6; idx_localStd = 7; idx_localMinMax = 8;
idx_state = 9;

if nargin > 0
    kernel = varargin{1};
elseif nargin == 0
    kernel = 'gaussian';
end

disp(['Now training your ' kernel ' SVM... ']);

datafile = 'Data/examples/data_recording_examples.mat';
load(datafile, 'recording_examples');
manual_data = struct2data(recording_examples(1).data{8,2});
%manual_data = data_accuracy{2,2}; % Jacob's manual tracking data
featCols = [idx_speed, idx_localSpeed, idx_localStd, idx_localMinMax]; 
%[4 6 7 8]; %[speed localSpeed localStd min(maxFwd, maxBkwd)]

%%% Prepare Training Data %%%
Mtrain = manual_data;

dataTrain = [];
for idx = [featCols idx_state]
    [thinData,~] = unshapeData(Mtrain,idx);
    dataTrain = [dataTrain thinData];
end

statesTrain = dataTrain(:,end);

%%% Remove nonlocust states %%%
idx_black = (statesTrain == -5);
dataTrain(idx_black,:) = [];
statesTrain = dataTrain(:,end); %remove them from the actual labels/states
%%% Change rotating (purple) to stationary (red) %%%
statesTrain( statesTrain == -1 ) = 0; % 

%%% Fit an SVM %%%
X = dataTrain(:,1:end-1);
labels = statesTrain; %dataTrain(:,end);
% a linear classifier
tLin = templateLinear('Learner','svm','Lambda',0.257);

% an optimized SVM classifier
% optimization done in SVM.m
tSVM = templateSVM( 'KernelFunction','gaussian',...
                    'Standardize', false, 'BoxConstraint', 1.0515, 'KernelScale', 3.3579 );
coding = 'onevsall';

tic
if strcmp(kernel,'gaussian')
    %%% state of the art %%%
    model = fitcecoc(X,labels,'Learners',tSVM,'Coding',coding);
elseif strcmp(kernel,'linear')
    %%% linear version %%%
    model = fitcecoc(X',labels','Learners',tLin,'Coding','onevsone','ObservationsIn','columns');
end
fprintf(['Fitting your %s SVM took %f seconds.', newline],kernel, toc)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Associated Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [thinData, idx_data] = unshapeData(M,idx_feature)
    %%% inputs
    %feature = nlocusts by 1 by ntimesteps matrix
    %%% outputs
    %a one dimensional array with length nlocusts*(ntimesteps) - (#ofNaNvalues)

    newData = M(:,idx_feature,2:end); %leave off first frame since no motion state can be assigned

    % to throw out data points that don't have well-defined speeds
    % these are the points at the start and end of tracks (I think?)
    idx_speed = 4;
    speed = M(:,idx_speed,2:end);
    flag = isnan(speed);
    newData(flag) = NaN;
    
    thinData = newData(:);
    idx_data = ~isnan(thinData);
    thinData = thinData(idx_data);
end 

function wideData = shapeData(thinData, idx_data, s)
    % recreates old shape
    newData = NaN(length(idx_data),1);
    newData(idx_data) = thinData;
    wideData = reshape(newData, s(1), 1, s(3)-1);
end

end