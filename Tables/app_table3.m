close all
% clear all % ensure we're loading new data
addpath('../Functions',...
        '../Functions/packages/CircStat',...
        '../Data/')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

global  idx_x idx_y idx_flag...
        idx_speed idx_theta...
        idx_localSpeed idx_localStd idx_localMinMax...
        idx_state...
        numNeighbors numFocal

%%% Options %%%
assembleData = 1; %set to 0 to load existing table data
printStats = 1; %set to 0 to surpress printing motion stats as you go
%%%

tableDataFile = 'table_data.mat';

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
% the loaded file includes the following variables:
%   1) recording = a struct array
%   2-10) idx_* = the global feature indices for data_final

%%% Set up the Tables %%%
%Construct a new table, using metadata from recording
    nRecs = size(recording,2);
    
    %Process the row names
    rowNames = {'All'};
    for j = 1:nRecs
        recName = strrep(recording(j).file, '_', '\_');
        rowNames{end+1} = recName;
    end
    for j = 1:nRecs
        clips = recording(j).data(:,1);
        clips = cellfun(@(C) strrep( C , '_', '\_'), clips, 'UniformOutput',false);
        rowNames(end+1:end+numel(clips)) = clips;
    end
    rowNames = rowNames';
    
    %Generate the column names
    varBand = 'Band';
    
    nall = '$n$ (All)';
    M1all = '$\vert \mathbf{M}_1 \vert$ (All)';
    phi1all = '$\phi_1$ (All)';
    Ms1all = '$-M^s_1$ (All)';
    Mc1all = '$M^c_1$ (All)';
    M2all = '$\vert \mathbf{M}_2 \vert$ (All)';
    phi2all = '$\phi_2$ (All)';
    Ms2all = '$-M^s_2$ (All)';
    Mc2all = '$M^c_2$ (All)';
    all_vars = {nall, M1all, phi1all, Ms1all, Mc1all, M2all, phi2all, Ms2all, Mc2all};
    
    nstop = '$n$ (Stop)';
    M1stop = '$\vert \mathbf{M}_1 \vert$ (Stop)';
    phi1stop = '$\phi_1$ (Stop)';
    Ms1stop = '$-M^s_1$ (Stop)';
    Mc1stop = '$M^c_1$ (Stop)';
    M2stop = '$\vert \mathbf{M}_2 \vert$ (Stop)';
    phi2stop = '$\phi_2$ (Stop)';
    Ms2stop = '$-M^s_2$ (Stop)';
    Mc2stop = '$M^c_2$ (Stop)';
    stop_vars = {nstop, M1stop, phi1stop, Ms1stop, Mc1stop, M2stop, phi2stop, Ms2stop, Mc2stop};
    
    ncrawl = '$n$ (Walk)';
    M1crawl = '$\vert \mathbf{M}_1 \vert$ (Walk)';
    phi1crawl = '$\phi_1$ (Walk)';
    Ms1crawl = '$-M^s_1$ (Walk)';
    Mc1crawl = '$M^c_1$ (Walk)';
    M2crawl = '$\vert \mathbf{M}_2 \vert$ (Walk)';
    phi2crawl = '$\phi_2$ (Walk)';
    Ms2crawl = '$-M^s_2$ (Walk)';
    Mc2crawl = '$M^c_2$ (Walk)';
    crawl_vars = {ncrawl, M1crawl, phi1crawl, Ms1crawl, Mc1crawl, M2crawl, phi2crawl, Ms2crawl, Mc2crawl};

    nhop = '$n$ (Hop)';
    M1hop = '$\vert \mathbf{M}_1 \vert$ (Hop)';
    phi1hop = '$\phi_1$ (Hop)';
    Ms1hop = '$-M^s_1$ (Hop)';
    Mc1hop = '$M^c_1$ (Hop)';
    M2hop = '$\vert \mathbf{M}_2 \vert$ (Hop)';
    phi2hop = '$\phi_2$ (Hop)';
    Ms2hop = '$-M^s_2$ (Hop)';
    Mc2hop = '$M^c_2$ (Hop)';
    hop_vars = {nhop, M1hop, phi1hop, Ms1hop, Mc1hop, M2hop, phi2hop, Ms2hop, Mc2hop};


    varNames = {varBand, all_vars{:}, stop_vars{:}, crawl_vars{:}, hop_vars{:} };
    %define variable types
    varTypes = cell(size(varNames));
    varTypes(:) = {'double'};
    varTypes(1) = {'string'};
    
% initialize the table 
sz = [numel(rowNames) numel(varNames)];

neighborTable = table('Size',sz,...
                'VariableTypes',varTypes,...
                'RowNames',rowNames,...
                'VariableNames',varNames);


%%% Extract indices for each state...
all_idx = cell2mat( cellfun(@(c) find( strcmp(varNames,c) ), all_vars, 'uni', false) );

stop_idx = cell2mat( cellfun(@(c) find( strcmp(varNames,c) ), stop_vars, 'uni', false) );

crawl_idx = cell2mat( cellfun(@(c) find( strcmp(varNames,c) ), crawl_vars, 'uni', false) );

hop_idx = cell2mat( cellfun(@(c) find( strcmp(varNames,c) ), hop_vars, 'uni', false) );



nClips = 0;
nRecs = size(recording,2);
for file_idx = 1:nRecs
    nClips = nClips + size(recording(file_idx).data, 1);
end
totalData = cell(nClips, 16);
clipNum = 1;


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
%matNums(5) = [];% 133 remove b/c low density
% matNums([5 7]) = [];% 098 remove b/c low density, 
% 4 5 6 7 are all a bit odd, but 5 and 7 are the worst
%matNums = 1; % for testing 
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

vmc = 4; %vidMotionCells
%{ [posx posy], angles, FC1stats, FC2stats }
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
        %d = 7; % 7 cm max radius for angles
        idx_closeNeigh = ( vecnorm(thisData(:,1:2),2,2) < d );
        % filtered out via & indexing with idx_motion
        %   FCstats = 1 by 2 cell = { nAngles,  {(mode1) magM, phi, -Ms, Mc},
        %                                       {(mode2) magM, phi, -Ms, Mc}  }
        
        allData{m,1} = thisData(:,1:2);
        allData{m,2} = neighborAngles( idx_closeNeigh );
        FCstats = computeFCstats(allData{m,2});
        %FC2stats?
        allData(m,3:4) = FCstats(2:3);

        
        idx_stop = (thisData(:,3) == 0);
        stopData{m,1} = thisData(idx_stop,1:2);
        stopData{m,2} = neighborAngles(idx_stop & idx_closeNeigh);
        FCstats = computeFCstats(stopData{m,2});
        %FC2stats?
        stopData(m,3:4) = FCstats(2:3);
        
        idx_crawl = (thisData(:,3) == 1);
        crawlData{m,1} = thisData(idx_crawl,1:2);
        crawlData{m,2} = neighborAngles(idx_crawl & idx_closeNeigh);
        FCstats = computeFCstats(crawlData{m,2});
        %FC2stats?
        crawlData(m,3:4) = FCstats(2:3);
        
        idx_hop = (thisData(:,3) == 2);
        hopData{m,1} = thisData(idx_hop,1:2);
        hopData{m,2} = neighborAngles(idx_hop & idx_closeNeigh);
        FCstats = computeFCstats(hopData{m,2});
        %FC2stats?
        hopData(m,3:4) = FCstats(2:3);
        
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

% compute FC stats for this video
FCstats = computeFCstats(allData(:,2));
neighborTable(1+file_idx, all_idx) = horzcat(FCstats(1), FCstats{2}, FCstats{3}) ;

FCstats = computeFCstats(stopData(:,2));
neighborTable(1+file_idx,stop_idx) = horzcat(FCstats(1), FCstats{2}, FCstats{3}) ;

FCstats = computeFCstats(crawlData(:,2));
neighborTable(1+file_idx,crawl_idx) = horzcat(FCstats(1), FCstats{2}, FCstats{3}) ;

FCstats = computeFCstats(hopData(:,2));
neighborTable(1+file_idx,hop_idx) = horzcat(FCstats(1), FCstats{2}, FCstats{3}) ;

end

% assign FC stats to table
idxes = [all_idx; stop_idx; crawl_idx; hop_idx];

for i = 1:size(idxes,1)
    FC1stats = vertcat( totalData{:,3+vmc*(i-1)} );
    FC2stats = vertcat( totalData{:,4+vmc*(i-1)} );

    idx = 1+nRecs+1:nRecs+1+size(FC1stats,1);
    nAngles = num2cell(cellfun(@(C) length(rmmissing(C)), totalData(:,2+vmc*(i-1))));
    neighborTable(idx,idxes(i,:)) = horzcat(nAngles, FC1stats, FC2stats );
end

%%% Longform... %%%
% 
% FC1stats = vertcat( totalData{:,3} );
% FC2stats = vertcat( totalData{:,4} );
% nAngles_all = num2cell(cellfun(@(C) length(rmmissing(C)), totalData(:,2)));
% neighborTable(idx,all_idx) = horzcat(nAngles_all, FC1stats, FC2stats );
% 
% FC1stats = vertcat( totalData{:,3+vmc} );
% FC2stats = vertcat( totalData{:,4+vmc} );
% nAngles_stop = num2cell(cellfun(@(C) length(rmmissing(C)), totalData(:,2+vmc)));
% neighborTable(idx,stop_idx) = horzcat(nAngles_stop, FC1stats, FC2stats );
% 
% FC1stats = vertcat( totalData{:,3+2*vmc} );
% FC2stats = vertcat( totalData{:,4+2*vmc} );
% nAngles_crawl = num2cell(cellfun(@(C) length(rmmissing(C)), totalData(:,2+2*vmc)));
% neighborTable(idx,crawl_idx) = horzcat(nAngles_crawl, FC1stats, FC2stats );
% 
% FC1stats = vertcat( totalData{:,3+3*vmc} );
% FC2stats = vertcat( totalData{:,4+3*vmc} );
% nAngles_hop = num2cell(cellfun(@(C) length(rmmissing(C)), totalData(:,2+3*vmc)));
% neighborTable(idx,hop_idx) = horzcat(nAngles_hop, FC1stats, FC2stats );

% compute FC stats for all videos!
for i = 1:size(idxes,1)
    FCstats = computeFCstats(totalData(:,2+vmc*(i-1)));
    neighborTable(1,idxes(i,:)) = horzcat(FCstats(1), FCstats{2}, FCstats{3});
end

%%% Outdated longform... %%%
% %all motion
% FCstats = computeFCstats(totalData(:,2));
% neighborTable(1,all_idx) = {FCstats{1}{:} FCstats{2}{2:end}};
% %stop
% FCstats = computeFCstats(totalData(:,2+vmc));
% neighborTable(1,stop_idx) = {FCstats{1}{:} FCstats{2}{2:end}};
% %crawl
% FCstats = computeFCstats(totalData(:,2+2*vmc));
% neighborTable(1,crawl_idx) = {FCstats{1}{:} FCstats{2}{2:end}};
% %hop
% FCstats = computeFCstats(totalData(:,2+3*vmc));
% neighborTable(1,hop_idx) = {FCstats{1}{:} FCstats{2}{2:end}};

% convert to a matrix array
allData = cell2mat(totalData(:,1));
stopData = cell2mat(totalData(:,1+vmc));
crawlData = cell2mat(totalData(:,1+2*vmc));
hopData = cell2mat(totalData(:,1+3*vmc));

allAngle = cell2mat(totalData(:,2));
stopAngle = cell2mat(totalData(:,2+vmc));
crawlAngle = cell2mat(totalData(:,2+2*vmc));
hopAngle = cell2mat(totalData(:,2+3*vmc));


%%% trigTable %%%

n = '$N$ (angles)';
M1 = '$\vert \mathbf{M}_1 \vert$';
phi1 = '$\phi_1$';
Ms1 = '$-M^s_1$';
Mc1 = '$M^c_1$';
M2 = '$\vert \mathbf{M}_2 \vert$';
phi2 = '$\phi_2$';
Ms2 = '$-M^s_2$';
Mc2 = '$M^c_2$';
vars = {n, M1, phi1, Ms1, Mc1, M2, phi2, Ms2, Mc2};
types = cell(1,numel(vars)); types(:) = {'double'};
trigTable = table('Size',[4 numel(vars)],...
                  'VariableTypes',types,'VariableNames',vars,...
                  'RowNames',{'All','Stationary','Walking','Hopping'});
for i = 1:4
    trigTable(i,:) = neighborTable(1,idxes(i,:));
end

%%% Stat tests %%%
% p = circ_otest(alpha)

pall = circ_otest(allAngle);
pstop = circ_otest(stopAngle);
pcrawl = circ_otest(crawlAngle);
phop = circ_otest(hopAngle);

fprintf([
    'p-values for omnibus test for uniformity:\n '...
    'all states: %e | stop: %e | crawl: %e | hop: %e \n'...
        ],...
    pall, pstop, pcrawl, phop)

save(tableDataFile,'trigTable',...
                   '-append')
else
   load(tableDataFile, 'trigTable') 
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

%%% Computing discrete Fourier coefficients

function FCstats = computeFCstats(angles)
%%% input: 
%       angles = distribution of angles
%       mode = the Fourier mode you're interested in
%%% output:
    %   FCstats = 1 by 2 cell = { nAngles,  {(mode1) magM, phi, -Ms, Mc},
    %                                       {(mode2) magM, phi, -Ms, Mc}  }
    
    FCstats = cell(1,2);
    
    if iscell(angles)
        angles = cell2mat(angles);
    end
    angles = rmmissing(angles);
    nAngles = length(angles);
    N = 50; %How many random samples to draw for bootstrapping
    
    FCstats{1} = nAngles;

    % number of data points in each sample = number of data points in original data set
    % [Mcomplex magM_p phi_p] = circ_moment(alpha, [weights], p=mode, centerFlag=false, dim=longest)
%     bootMs = bootstrp(N, @(samples) {   circ_moment(samples,[],1),...
%                                         circ_moment(samples,[],2)   }, angles);
%                                     %FourierCoef(samples,'c',mode) FourierCoef(samples,'s',mode)], angles);
    for mode = 1:2
        [M, magM, phi] = circ_moment(angles,[],mode);
        
%         bootMagM = abs( [bootMs{:,mode}] );
%         bootPhi = angle( [bootMs{:,mode} ] );
% 
%         meanMagM = mean( bootMagM );
%         stdMagM = std( bootMagM );
% 
%         meanPhi1s = circ_mean( bootPhi );
%         stdPhi1 = circ_std( bootPhi );
        
%         if mode == 1
%             Msc = -imag(M);
%         elseif mode == 2
%             Msc = real(M);
%         end
        Ms = imag(M);
        Mc = real(M);

        FCstats{mode+1} = {magM, phi, -Ms, Mc};
    end
end

% function M = FourierCoef(angles, trigFlag, mode)
%     N = length(angles);
%     if trigFlag == "s"
%         M = (1/N)*sum(sin(mode*angles));
%     elseif trigFlag == "c"
%         M = (1/N)*sum(cos(mode*angles));
%     else
%         error("Enter c or s for cosine or sine function.")
%     end
% end
