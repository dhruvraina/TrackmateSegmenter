% TrackMate segmenter:
% Author: d.raina
% last edit: 11.April.17

%PseudoCode:

%Read in trackmate 'matlab' file. Use each point as a "seed" point to
%perform nuclear segmentation. Read out median intensities. Currently this
%works by using matlab's 'MultiThresh' applied locally. KMSeg is also an
%option to implement!

close all; clear all
%ToDo: File i/o ??  ::c_dir should be passed through the datloader-11/04/17

%% Inputs:
c_dir = '/Users/draina/Desktop/Draq5_Test/';
file = [c_dir 'InvertedTrack_matlab.xml'];

[trackList trackMeta] = importTrackMateTracks(file);

mainImg = tiffread2([c_dir '5mMDRAQ_Washout_200_m2 copy.tif']);
scalingFactor = 0.361;  %:: Scaling factor is derived from image Metadata. TBD: automatically extract this number.


%% For Inverted Tracks:
invertmode = 'y';
framestart = 1;
frameend = length(mainImg);
if invertmode =='y'
    
    %Make a backup of the original variable for debugging
    templist = trackList;
    
    %Convert the frame numbers into reverse order
    gg = cellfun(@(x) frameend-x(:,1), trackList, 'UniformOutput', 0);
    
    %Replace reversed frame numbers in original trackList
    for cc = 1:length(trackList)
        trackList{cc,1}(:,1)= gg{cc,1}(:,1);
    end
    
    %Re-sort accordind to new frame numbers:
    trackList = cellfun(@(x) sortrows(x,1), trackList, 'UniformOutput', 0);
    
end
keyboard


%% Framewise Segmentation and Datextraction:
for ctr1 = 1:length(mainImg)
    % ctr1 = 23;
    %Find all tracks which are active in this frame:
    c_img      = mainImg(ctr1).data;
    c_imgTrax1 = cellfun(@(x)find(x==ctr1),trackList(:,1), 'UniformOutput', false);
    c_imgTrax  = cellfun(@(x)~isempty(x), c_imgTrax1, 'UniformOutput', false);
    [c_imgTrax, q] = find(cell2mat(c_imgTrax));
    
    
    %%Debug:
    c_imgTrax2 = cellfun(@(x)find(x==ctr1),templist(:,1), 'UniformOutput', false);
    c_imgTrax3  = cellfun(@(x)~isempty(x), c_imgTrax2, 'UniformOutput', false);
    [c_imgTrax3, q] = find(cell2mat(c_imgTrax3));
    %%
    
    
    %% List XY Coords:
    for ctr2 = 1:length(c_imgTrax)
        tempIdx      = find(trackList{c_imgTrax(ctr2),1}(:,1)==ctr1);
        c_xy(ctr2,:) = [trackList{c_imgTrax(ctr2),1}(tempIdx,:)  c_imgTrax(ctr2)];
        c_xy(ctr2,2) = c_xy(ctr2,2)/scalingFactor;
        c_xy(ctr2,3) = c_xy(ctr2,3)/scalingFactor;
        
        %Trackmate sometimes sends along zeros as XY coords. Correcting for this:
        c_xy(c_xy(ctr2,2)==0,2)=1; 
        c_xy(c_xy(ctr2,3)==0,3)=1; 
    end
    
    %[xx, yy]  = meshgrid(1:512); %For circular mask
    % totalmask = zeros(size(c_img));
    
        debugoff = 0;
    switch debugoff
        
        case(0)
            %Debugging: Overwrite XY points raw image:
            temp_xy = [c_xy(:,2:3) c_xy(:,5)];  %Passing only XY and TrackID
            temp_xy = round(temp_xy+0.6);
            im_d1 = size(c_img,1);
            im_d2 = size(c_img,2);
            idxlist = sub2ind([im_d1 im_d2], temp_xy(:,2), temp_xy(:,1));
            lbl = c_img;
            lbl(idxlist) = 400;
            
        case(1)
            %% Watershed Thresholding:
            temp_xy = [c_xy(:,2:3) c_xy(:,5)];  %Passing only XY and TrackID
            try
                [lbl badtrax{ctr1}] = locwater2(c_img, temp_xy);
            catch
                keyboard
            end
            
            
            
            %Watershedding needs to be fixed
            %IMMEDIATE: xylocs don't track properly outside the image!! When the
            %spot moves outside the image, x,y is set to zero by TrackMate!!
            
            %1. Fix the fgm marker building to work more reliably (frame 220 for instance)
            %2. Force badpts to adopt a circular ROI, or try some sort of local
            %thresholding as implemented previously
            
            
            %% Data Extraction:
            try
                [tvec_med tvec_mod tvec_trax] = tmcalc(c_img, lbl);
            catch
                keyboard
            end
            
            %% Assigning into data structure:
            for ctr3 = 1:length(tvec_trax)
                resvec_seg(tvec_trax(ctr3)).median(ctr1) = tvec_med(ctr3);
                resvec_seg(tvec_trax(ctr3)).mode(ctr1) = tvec_mod(ctr3);
            end
    end
    
    clear tvec_med tvec_mod tvec_trax temp_xy c_xy c_imgTrax c_imgTrax1
    
    %Debug if loop
    
    %% Debugging:
    ff = figure('visible', 'off')
    hold on
    %subplot(1,2,1), imagesc(lbl), subplot(1,2,2), imagesc(c_img)
    imagesc(lbl)
    print(ff,[c_dir 'segimg' num2str(ctr1)], '-dpng','-r200')
    close gcf
    
    
    
end




