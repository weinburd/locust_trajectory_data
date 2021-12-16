function [LSC, FP, FN] = linkAccuracy(auto, manual, costUnmatched)
%%% Inputs:
%%%     auto = matrix of data from automatic tracking
%%%     manual = matrix of data from manual tracking
%%%     **Must have same number of timesteps = size(input, 3)
%%%     **Timesteps must line up perfectly
%%%     costUnmatched = the cost of not matching a link
%%%           This is the cost associated with leaving a link in 'auto' NOT matched to
%%%           one in 'manual'. Since our costs are Euclidian distances in R^4, this
%%%           basically says the following:
%%%           For a link X in 'auto', if there is no link in 'manual' that is within a
%%%           distance of costUnassigned centimeters from X, then leave X unassigned.
%%% Outputs:
%%%     JSC = Jacob Similarity Coefficient =  TP/(TP + FP + FN);

% isolate the position data
auto = auto(:,1:2,:);
manual = manual(:,1:2,:);

% Now we assemble our link matrices
autoLinks = auto(:,1:2,1:end-1) ;
autoLinks = cat(2,autoLinks, auto(:,1:2,2:end));

manualLinks = manual(:,1:2,1:end-1) ;
manualLinks = cat(2,manualLinks, manual(:,1:2,2:end));

%%Here we assemble a Cost matrix of distances from the auto to the manual
%%points. For example, Cost(i,j) is the Euclidean distance from point i to point j

%%The Cost matrix has the number of rows equal to the number of tracks in
%%the automatic tracks and the number of columns equal to the number of
%%tracks in the manual track. Its depth is the number of frames.
% We assign costs by treating a link as a point in R^4
% and using the Euclidian distance between links as the cost

Cost = zeros(size(autoLinks,1),size(manualLinks,1),size(autoLinks, 3));
for j = 1:size(autoLinks, 3)
    for k = 1:length(autoLinks)
        Cost(k,:,j) = vecnorm( [autoLinks(k,1,j)-transpose(manualLinks(:,1,j));...
                                autoLinks(k,2,j)-transpose(manualLinks(:,2,j));
                                autoLinks(k,3,j)-transpose(manualLinks(:,3,j));...
                                autoLinks(k,4,j)-transpose(manualLinks(:,4,j))]...
                                ,2,1)';
    end
end

%%Now we set all NaN values to a super high value in the cost matrix, so
%%they will not be matched with anything. So long as they are higher than
%%the unassignment cost, I think this will be fine. This is just because
%%matchpairs can't handle NaN values. Match is a matrix where each row
%%shows the indices of paired points in automatic and manual track. So if a
%%row is [3,4], then the point in the 3rd row of autoLinks corresponds to
%%the point in the 4th row of testman.
Cost(isnan(Cost))=1000;

%unassignment = .5; % this is a distance (in cm)... beyond which points won't be matched?
unassignment = costUnmatched;

%for m = 0:.01:1.25
    %unassignment = m;
    Match = NaN(length(manualLinks),2,size(autoLinks,3));
    for j = 1:size(autoLinks, 3)
        singleMatch = matchpairs(Cost(:,:,j),unassignment);
        Match(1:size(singleMatch,1),:,j) = singleMatch;
    end

    %%Now I'll make a matrix to find the number of points present in both the
    %%manual and automatic tracks at each time step. Row 1 will be the number
    %%in the automatic, Row 2 will be the number in the manual, Row 3 will
    %%be the number of paired spots, Row 4 will be the number of false
    %%positives (Row 1 minus Row 3), andRow 5 will be the number of false
    %%negatives (Row 2 minus Row 3).
    PointsPresent = zeros(5,1,size(autoLinks,3));
    for j = 1:size(autoLinks, 3)
        PointsPresent(1,1,j) = sum(~isnan(autoLinks(:,1,j)));
        PointsPresent(2,1,j) = sum(~isnan(manualLinks(:,1,j)));
        PointsPresent(3,1,j) = sum(~isnan(Match(:,1,j)));
        PointsPresent(4,1,j) = PointsPresent(1,1,j)-PointsPresent(3,1,j);
        PointsPresent(5,1,j) = PointsPresent(2,1,j)-PointsPresent(3,1,j);
    end
    
    %%JS will be our version of the Jaccard Similarity Coefficient. It will be
    %%equal to the true positives, i.e., the total number of matched points
    %%within the unassignment distance (TP), divided by the sum of the true
    %%positives plus the number of false positives (FP) plus the number of
    %%false negatives (FN).
    TP = sum(PointsPresent(3,1,:),3);
    FP = sum(PointsPresent(4,1,:),3);
    FN = sum(PointsPresent(5,1,:),3);
    LSC = TP/(TP + FP + FN);
end