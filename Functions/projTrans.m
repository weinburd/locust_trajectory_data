function [trans, scale, fieldDims, newR_A] = projTrans(cornersPix, fieldDimsPix)
%%% Compute a projective transformation to apply to the data %%%
    % the width and length of the actual board in cm
    wid = 120.5985767; hght = 59.56373328; %see Video Details GSheet
    rect = [ wid 0; 0 0; 0 hght; wid hght]; %a rectangle
    
    % the four corners of the board in the video, depends on the video
    % NOTE: we can have negative numbers or numbers greater than the # pixels;
    % in some cases we had infer the location of a corner when the board was
    % only partially in the video field.
    %any quadrilateral [  x_tr y_tr; x_tl y_tl; x_bl y_bl; x_br y_br]
    %                  [ top right; top left; bottom left; bottom right]
    quad = cornersPix;
    
    tform = fitgeotrans(quad, rect, 'projective');
    trans = tform.T;
    
    %%% This part is in the function dataClean.m
    % Pquad = [quad [1;1;1;1]];
    % Prect = Pquad*trans;
    % newrect = Prect(:,1:2)./Prect(:,3);
    
    %%% Compute an average scaling of the area in the video %%%
    w = fieldDimsPix(2)-fieldDimsPix(1); 
    h = fieldDimsPix(4)-fieldDimsPix(3);
    oldR = [w 0; 0 0; 0 h; w h];
    oldR_A = quadArea(oldR);
    
    projR = [oldR ones(4,1)];
    newProjR = projR*trans;
    newR = newProjR(:,1:2)./newProjR(:,3);
    newR_A = quadArea(newR);
    
    scale = sqrt(oldR_A/newR_A); % pix/cm
    
    %now transform the video field dimensions
%     projDims = [fieldDimsPix 1];
%     projNewDims = projDims*trans;
%     fieldDims = [projNewDims(1) projNewDims(2)]/projNewDims(3);
    
    %fieldDims = [max(newR(:,1))-min(newR(:,1)) max(newR(:,2))-min(newR(:,2))];
    % These are the DIMENSIONS of the smallest rectangle containing the
    % transformed rectangle newR. But our Particle in Cell method for 
    % counting neighbors uses a box with these dimensions AND one corner at
    % the origin. Thus, we need to pass it a bigger box!
    fieldDims = [min(newR(:,1)) max(newR(:,1)) min(newR(:,2)) max(newR(:,2))];
    
    function quadA = quadArea(corners)
       
        tr = corners(1,:); tl = corners(2,:); bl = corners(3,:); br = corners(4,:);
        diagonal = tr - bl;
        left = tl-bl;
        theta1 = angle( (left(1)+1i*left(2))/(diagonal(1)+1i*diagonal(2)) );
        theta1 = - theta1; % due to the image's left-handed coord system
        right = br-tr;
        theta2 = angle( (right(1)+1i*right(2))/(-diagonal(1)-1i*diagonal(2)) );
        theta2 = -theta2; % due to the image's left-handed coord system
        h1 = norm(left)*sin(theta1);
        h2 = norm(right)*sin(theta2);
        
        quadA = 1/2*norm(diagonal)*(h1+h2); % area of the board in pix^2 
    end
end