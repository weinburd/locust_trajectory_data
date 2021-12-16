function data_final = dataForm(data_clean,model,smoothFlag, interpThetaFlag, varargin)
%%% Input:
% data_clean = an Nlocusts by 3 by Ntimesteps matrix
%               each locust has [x y flag] at each timestep
%               the spatial units are centimeters
%               the flag values are:
%                   0 = absent data because locust is not in frame
%                   1 = data from when a locust is present
%                   2 = missing data to be interpolated (track had a missing spot)
% smoothFlag = logical 
%               True - use smoothed position data when computing velocity
% optional inputs:
% manualFlag = is a flag to import the data including motion state data
% plotFlag = a flag whether or not to plot a scatter plot colored by state
%%% Output:
% data_final = an Nlocust by 9 by Ntimestep matrix
%       each locust has 9 features  at each timestep
%       the features are:
%   [x y flag speed theta MagAveVelocity localStdSpeed min(fwdMax,bkwdMax) motionState]
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
%       column 9 is motionState, either inferred from SVM or manually assigned (for manually tracked data sets)
disp('Now computing inferred data');
tic
nOptInputs = numel(varargin);
if nOptInputs > 1
    manualFlag = varargin{2};
    plotFlag = varargin{1};
elseif nOptInputs > 0
    manualFlag = false;
    plotFlag = varargin{1};
else
    manualFlag = false;
    plotFlag = false;
end

Mfinal3D = data_clean;
M = zeros(size(Mfinal3D));

global  idx_x idx_y idx_flag...
        idx_speed idx_theta...
        idx_localSpeed idx_localStd idx_localMinMax...
        idx_state%...
%        idx_speed_NS
% idx_speed_NS = 10;

idx_x = 1; idx_y = 2; idx_flag = 3;
idx_speed = 4; idx_theta = 5;
idx_localSpeed = 6; idx_localStd = 7; idx_localMinMax = 8;
idx_state = 9;

if manualFlag
    M(:,idx_x:idx_y,:) = Mfinal3D(:,1:2,:);
    M(:,idx_flag,:) = Mfinal3D(:,4,:); % This is the column of flags
    states = Mfinal3D(:,3,:);
else
    M = Mfinal3D; %%set M equal to your imported and cleaned 3D matrix
end
s = size(M);
    
%%separating xpos, ypos, and creating xvel yvel and total velocity
%Add a speed and direction (rad) column to M
if smoothFlag
    KeepNaN = double(~isnan(M(:,1,:)));
    KeepNaN(KeepNaN==0) = NaN;
    xsmooth = smoothdata(M(:,1,:),3,"gaussian",8,'omitnan');
    ysmooth = smoothdata(M(:,2,:),3,"gaussian",8,'omitnan');
    xsmooth = xsmooth.*KeepNaN;
    ysmooth = ysmooth.*KeepNaN;
    xpos = xsmooth;
    ypos = ysmooth*-1;%just to flip the axis so up is positive

    %flag = M(:,3, :);
    [~,~,xvel] = gradient(xpos);
    [~,~,yvel] = gradient(ypos);
    xvel = xvel*25; %convert to cm/s
    yvel = yvel*25; %convert to cm/s

    speed = (xvel.^2 + yvel.^2).^.5; %speed of locust in cm/frm
    direct = atan2(yvel,xvel); %direction locust is moving (rad)
else
    xpos_NS = M(:,1,:);
    ypos_NS = M(:,2,:);
    ypos_NS = ypos_NS*-1;

    %flag = M(:,3, :);
    [~,~,xvel_NS] = gradient(xpos_NS);
    [~,~,yvel_NS] = gradient(ypos_NS);
    xvel_NS = xvel_NS*25; %convert to cm/s
    yvel_NS = yvel_NS*25; %convert to cm/s

    speed = (xvel_NS.^2 + yvel_NS.^2).^.5; %speed of locust in cm/frm
    direct = atan2(yvel_NS,xvel_NS); %direction locust is moving (rad)
                                     %atan2 outputs values in [-pi,pi]
end

vvec = [cos(direct),sin(direct)];

idx_zerospeed = speed == 0;
direct(idx_zerospeed) = NaN; %get rid of 0,0 velocities
z = (speed>.01); %logical matrix that keeps represents only non-zero velocities
trackVel = mean(speed,3,'omitnan');

M(:,idx_theta,:) = direct;
M(:,idx_speed,:) = speed;

% First Spike Ends: mean = 4.5, 5.2; median = 2, 4; std = 4.7, 4.9
% Second Spike Ends: mean = 16.0, 17.8; median = 12, 14; std = 14.0, 15.4
% Could choose: Correlation time = median + 1 std = 7 to 9
% only want to associate data points with others within this correlation time
% so set l = 9; k = 9; j = 2*9? Shouldn't these be reversed? i.e. l = k = 2*9; j=9?
l = 15; %9; %17
k = 15; %9; %19
j = 15; %18; %20
smoothedPos = 0;
[localSpeed, localStd, minMax]  = localStats(M, l, k, j, smoothedPos);
M(:,idx_localSpeed,:) = localSpeed;
M(:,idx_localStd,:) = localStd;
M(:,idx_localMinMax,:) = minMax;

% %%Using locally-averaged instantaneous velocities to determine if a locust is stopped,
% %%walking, or hopping. Outputs a num locusts x 2 x num frames matrix, where
% %%the first column is 0 = stopped, 1 = walking, or 2 = hopping and the
% %%second column is the time
% go = .0533;%% locusts with low standard deviations (less than dev) and locally-averaged velocities greater than this (in cm/s) are considered walking
% dev = 6;%%locally-averaged standard deviations greater than this (in cm/s) are considered hopping
% moving = zeros(s(1),1,s(3));
% moving(localStd>dev) = 2;%set the hoppers
% moving(localSpeed>go & localStd<=dev) = 1;%set the walkers
% moving(isnan(speed))= NaN; %set the NaNs
% moving(moving==10)=NaN; %what is this doing???
% placeHold = NaN(size(moving,1),1,size(moving,3));
% moving = cat(2,moving, placeHold);
% for tm = 0:size(moving,3)-1
%     appendable = tm*ones(s(1),1);
%     moving(:,2,tm+1) = appendable;
% end
% ph = moving(:,2,:);
% ph(flag) = NaN;
% moving(:,2,:)= ph;
% 

% sch = moving(:,1,:); %this is a num locusts x 1 x num frames matrix. Each value corresponds to stop (0), walk(1), or hop(2) for a locust at a given frame
% M = cat(2,M,sch);%appending sch to M

%%% Assign Motion States %%%
if manualFlag
    states( states == -65536 ) = 0;         % stationary (red)
    states( states == -65281 ) = -1; %0     % rotating (purple)
    states( states == -16776961 ) = 1;      % crawling (blue)
    states( states == -16711936 ) = 2;      % hopping (green)
    states( states == -16777216 ) = -5; %-1 % not a locust (black)
    
    M(:,idx_state,:) = states;
    plotFlag = 0; % plot function here needs reshaped data as below
else
    featCols = [idx_speed idx_localSpeed idx_localStd idx_localMinMax];
    % unshape and stick data together
    dataTest = [];
    for idx = featCols
        [thinData, idx_data] = unshapeData(M,idx);
        dataTest = [dataTest thinData];
    end
    %%% Predict States for Test Data %%%
    SVM_states = predict(model,dataTest','ObservationsIn','columns');
    states = shapeData(SVM_states,idx_data,s);
    states = cat(3, NaN(s(1),1,1), states);
    M(:,idx_state,:) = states;
end
    
%%% Interpolate Heading Direction %%%
%%% (Based on motion state) %%%

if interpThetaFlag
%starts as 0 tag, go until non zero tag
s = size(M);

%%% An ALTERNATIVE Heading Interpolation method that I never completed
%%% I realized it was unnecessary to redo it.
% stationary = reshape( M(:,idx_state,:), s(1), s(3) );
% [locusts, times] = find( stationary == 0 );
% pairs = [locusts times];
% pairs = sortrows(pairs, 1);
% Npairs = size( pairs, 1);
% 
% clear idx
% j = 1;
% for idx = 1:Npairs
%     locust = pairs(idx,1);
%     time = pairs(idx,2);
%     
%     lastLocust = pairs( max(1,idx-1),1 );
%     lastTime = pairs( max(1,idx-1),2 );
%     
%     if locust == pairs(max(idx-1,1),1) && time-1 == pairs(max(idx-1,1),2)
%         j = j+1;
%     else
%         %interpolate
%         tPreStat = max(lastTime-j,1);
%         tPostStat = min(lastTime+1,s(3));
%         if tPostStat-tPreStat==j
%             theta1 = M(lastLocust,idx_theta,tPreStat);
%             theta2 = M(lastLocust,idx_theta,tPostStat);
% 
%             
%             M(lastLocust,idx_theta,tPreStat+1:tPostStat-1) = interpThetas(theta1,theta2,j);
%         end
%         
%         j = 1; % reset j
%     end
%     
% end
% count = 0 % for debugging

for locust = 1:s(1)
    for timestep = 1:s(3)
        % if the locust is stationary, then interpolate a new position for it
        if M(locust, idx_state, timestep) == 0 
            ahead = sum(M(locust, idx_state, timestep:end), 'all','omitnan')>0;
            behind = sum(M(locust, idx_state, 1:timestep), 'all','omitnan')>0;
            %if ~ahead && ~behind %if never moving, eliminate
                %M(locust,:,:) = [];
            if ahead && ~behind %if only future position data, set initial angles to that future angle
                k = find(M(locust,idx_state,timestep+1:end),1,'first');
                M(locust, idx_theta, timestep) = M(locust,idx_theta,timestep+k);
            elseif ~ahead && behind %if only previous data, set angle to past angle
                M(locust, idx_theta, timestep) = M(locust, idx_theta, timestep-1);
            elseif ahead && behind %if previous and future data, interpolate 
                %disp("alert")
                k = find(M(locust,idx_state,timestep+1:end),1,'first');
                b = find(M(locust,idx_state,1:timestep),1,'last');
                %disp(timestep)
                %disp(k)
                %disp(b)
                j = k+timestep-b;
                th_last = M(locust,idx_theta,timestep-1);
                th_start = M(locust,idx_theta,b);
                th_end = M(locust,idx_theta,timestep+k);

                newTheta = interpTheta(th_start,th_end,th_last,j);
                
                
%                 %for debugging
%                 oldTheta = M(locust,idx_theta,timestep-1)+(M(locust,idx_theta,timestep+k)-M(locust,idx_theta,b))/(k+timestep-b);
%                 if abs(newTheta - oldTheta) > 1e-10
%                     th_start-th_end
%                     newTheta
%                     oldTheta
%                 end

                M(locust, idx_theta, timestep) = newTheta; 
            end
        end
      
        % there's some question as to where we should interpolate thetas
        % option 1: all stationary locusts (see Jacob's interp above)
        % option 2: all locusts who do have a speed, but not a well-defined by velocity (implemented below)
        % option 3: could do both (would need to modify so that we use the correct angles for when locusts begin moving again)

        % if the locust is detected AND the speed is well-defined AND the theta value is NaN, 
        % then interpolate a position for it
        if logical( M(locust, idx_flag, timestep) )...   
           && ~isnan( M(locust, idx_speed, timestep) )...
           && isnan( M(locust, idx_theta, timestep) )
           % order here matters for performace because Matlab stops evaluating as soon as it hits a "false"

            % ahead = 0 if there are only NaNs ahead, 1 otherwise. behind is similar
            ahead = sum( ~isnan(M(locust, idx_theta, timestep:end)), 'all','omitnan')>0;
            behind = sum( ~isnan(M(locust, idx_theta, 1:timestep)), 'all','omitnan')>0;
            %if ~ahead && ~behind %if never moving, eliminate
                %M(locust,:,:) = [];
            if ahead && ~behind %if only future position data, set initial angles to that future angle
                k = find( ~isnan(M(locust,idx_theta,timestep+1:end)),1,'first');
                M(locust, idx_theta, timestep) = M(locust,idx_theta,timestep+k);
            elseif ~ahead && behind %if only previous data, set angle to past angle
                M(locust, idx_theta, timestep) = M(locust, idx_theta, timestep-1);
            elseif ahead && behind %if previous and future data, interpolate 
                %disp("alert")
                k = find( ~isnan(M(locust,idx_theta,timestep+1:end)),1,'first');
                b = find( ~isnan(M(locust,idx_theta,1:timestep)),1,'last');
                    %disp(timestep)
                    %disp(k)
                    %disp(b)
                % for debugging
                if M(locust,idx_theta,timestep-1) ~= M(locust,idx_theta,b)
                    disp("Last theta is not the theta given by find(..., 'last')")
                end
                j = k+timestep-b;
                th_last = M(locust,idx_theta,timestep-1);
                th_start = M(locust,idx_theta,b);
                th_end = M(locust,idx_theta,timestep+k);

                newTheta = interpTheta(th_start,th_end,th_last,j);
                
                
%                 %for debugging
%                 oldTheta = M(locust,idx_theta,timestep-1)+(M(locust,idx_theta,timestep+k)-M(locust,idx_theta,b))/(k+timestep-b);
%                 if abs(newTheta - oldTheta) > 1e-10
%                     th_start-th_end
%                     newTheta
%                     oldTheta
%                 end

                M(locust, idx_theta, timestep) = newTheta;
            end
            % for debugging
            %if isnan(M(locust, idx_theta, timestep)) && (ahead || behind)
            %    count = count + 1
            %end
        end

    end %for timesteps
end %for locusts
end %if interpThetaFlag

if plotFlag
    %%% Scatter plot %%%
    figNum = 31;
    h_SVM = scatterSVM(dataTest, SVM_states, figNum);
end

data_final = M;
fprintf(['That took %f seconds', newline],toc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Associated Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [localSpeed, localStd, minMax]  = localStats(M, l, k, j, smoothFlag)
    %%% inputs:
    %M = nlocusts by nfeatures by ntimesteps matrix
    %l = 15; %this is the window size for local velocity (the total window size is symmetrical about and includes the point of focus)
    %k = 15; %this is window size for standard deviation
    %j = 15; %this is window size for min(fwdMax, bkwdMax)
    
    %%% outputs:
    %three local statistics, each has size: nlocusts by 1 by ntimesteps
    
    %global  idx_x idx_y
    
    if smoothFlag
        %using unsmoothed data 
        KeepNaN = double(~isnan(M(:,idx_x,:)));
        KeepNaN(KeepNaN==0) = NaN;
        xsmooth = smoothdata(M(:,idx_x,:),3,"gaussian",8,'omitnan');
        ysmooth = smoothdata(M(:,idx_y,:),3,"gaussian",8,'omitnan');
        xsmooth = xsmooth.*KeepNaN;
        ysmooth = ysmooth.*KeepNaN;
        xpos = xsmooth;
        ypos = ysmooth*-1;%just to flip the axis so up is positive
        
    else
        %using unsmoothed data 
        xpos_NS = M(:,idx_x,:);
        ypos_NS = M(:,idx_y,:);
        ypos_NS = ypos_NS*-1;

        xpos = xpos_NS;
        ypos = ypos_NS;
    end

    [~,~,xvelns] = gradient(xpos);%creating non-smoothed velocity data
    [~,~,yvelns] = gradient(ypos);
    xvelns = xvelns*25;%convert to cm/s
    yvelns = yvelns*25;
    nsv = (xvelns.^2 + yvelns.^2).^.5; %non-smoothed speed
    localxVelocity = movmean(xvelns,l,3,'omitNaN');
    localyVelocity = movmean(yvelns,l,3,'omitNaN');

    localSpeed = (localxVelocity.^2 + localyVelocity.^2).^.5;%local mean non-smoothed velocity
    localStd = movstd(nsv,k,0,3,'omitNaN');
    fwdmax = movmax(nsv,[0 j],3,'omitNaN');
    bkmax = movmax(nsv,[j 0],3,'omitNaN');
    minMax = min(fwdmax,bkmax,'omitNaN');

    flag = isnan(nsv); %flag to make sure we don't assign time coordinates to points where we have no actual data
    localSpeed(flag) = NaN; %use flag to make sure we don't get average velocities where there is no data
    localStd(flag) = NaN;
    minMax(flag) = NaN;
end

function [thinData, idx_data] = unshapeData(M,idx_feature)
    %%% inputs
    %feature = nlocusts by 1 by ntimesteps matrix
    %%% outputs
    %a one dimensional array with length nlocusts*(ntimesteps) - (#ofNaNvalues)

    newData = M(:,idx_feature,2:end); %leave off first frame since no motion state can be assigned

    % to throw out data points that don't have well-defined speeds
    % these are the points at the start and end of tracks (I think?)
    this_speed = M(:,idx_speed,2:end);
    flag = isnan(this_speed);
    newData(flag) = NaN;
    
    thinData = newData(:);
    idx_data = ~isnan(thinData);
    thinData = thinData(idx_data);
    %length(thinData)
end 

function wideData = shapeData(thinData, idx_data, s)
    % recreates old shape
    newData = NaN(length(idx_data),1);
    newData(idx_data) = thinData;
    wideData = reshape(newData, s(1), 1, s(3)-1);
end

% function newThetas = interpThetas(theta1,theta2,j)
%     %returns a vector of j theta values, linearly spaced between theta1 and theta2
%     js = 1:j;
%     newThetas = theta1 + js/(j+1)*(theta2-theta1);
%     newThetas = reshape( newThetas, 1,1,j);
% end

function newTheta = interpTheta(theta1,theta2,theta0,j)
    % returns a new theta value
    % newTheta is computed by taking a step of size 1/j from theta0
    % along the line from theta1 to theta2
    
    if theta1-theta2 > pi
        theta1 = theta1-2*pi;
    elseif theta2-theta1 > pi
        theta2 = theta2-2*pi;
    end

    newTheta = theta0 + 1/j*(theta2-theta1);
    newTheta = angle(exp(1i*newTheta)); %
end

function h_SVM = scatterSVM(dataPred, SVM_states, figNum)

%     states_color = zeros( length(man_states), 3 );
%     stopped = man_states==0; crawling = man_states==1; hopping = man_states==2;
%     states_color( stopped , 1 ) = 255;   % red = [255 0 0]
%     states_color( crawling , 3 ) = 255; % blue = [0 0 255]
%     states_color( hopping, 2 ) = 255;    % green = [0 255 0]
    
    SVMstates_color = zeros( length(SVM_states), 3 );
    stopped = SVM_states==0; crawling = SVM_states==1; hopping = SVM_states==2;
    SVMstates_color( stopped , 1 ) = 255;   % red = [255 0 0]
    SVMstates_color( crawling , 3 ) = 255; % blue = [0 0 255]
    SVMstates_color( hopping, 2 ) = 255;    % green = [0 255 0]
    
    sizes = 10;
    alph=0.3;
    
%     figure(figNum)
%     
%     h_true = scatter3( dataPred(:,2), dataPred(:,3), dataPred(:,4), sizes, states_color,...
%                     'MarkerEdgeAlpha', alph, 'MarkerFaceAlpha', alph...
%                     );
%     %axis([0 90 0 40 0 70])
%     %axis([0 30 0 40 0 70])
%     title(['Motion Statistics Colored by Manually Assigned State (', student, ' Video)'])
%     %xlabel('Instantaneous Speeds')
%     xlabel('Mag. of Ave. Velocity')
%     ylabel('Standard Deviation of Speed')
%     zlabel('Min( fwdMax, bkwdMax )')
    
    figure(figNum)
    h_SVM = scatter3( dataPred(:,2), dataPred(:,3), dataPred(:,4), sizes, SVMstates_color,...
                    'MarkerEdgeAlpha', alph, 'MarkerFaceAlpha', alph...
                    );
    %axis([0 90 0 40 0 70])
    %axis([0 30 0 40 0 70])
    title(['Motion Statistics Colored by SVM Assigned State'])
    %xlabel('Instantaneous Speeds')
    xlabel('Mag. of Ave. Velocity')
    ylabel('Standard Deviation of Speed')
    zlabel('Min( fwdMax, bkwdMax )')

end

end