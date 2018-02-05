% TrackMate segmenter:
% Author: d.raina
% last edit: 11.April.17

%PseudoCode:

%Read in trackmate 'matlab' file. Use each point as a "seed" point to
%perform nuclear segmentation. Read out median intensities. Currently this
%works by using matlab's 'MultiThresh' applied locally. KMSeg is also an
%option to implement!

%ToDo: File i/o ??  ::c_dir should be passed through the datloader-11/04/17
c_dir = '/Users/draina/Desktop/Draq5_Test/';
file = [c_dir '5mM_invertedTime.xml'];

[trackList trackMeta] = importTrackMateTracks(file);

mainImg = tiffread2([c_dir '5mMDRAQ_Washout_200_m2.tif']);
scalingFactor = 0.361;  %:: Scaling factor is derived from image Metadata. TBD: automatically extract this number. 








%% Main:

for ctr1 = 1:length(mainImg)
    keyboard %Which ctr1 = FrameNo?
    c_img      = mainImg(ctr1).data;
    c_imgTrax1 = cellfun(@(x)find(x==ctr1),trackList(:,1), 'UniformOutput', false);
    c_imgTrax  = cellfun(@(x)~isempty(x), c_imgTrax1, 'UniformOutput', false);
    [c_imgTrax, q] = find(cell2mat(c_imgTrax));
    
    
    %List XY Coords:
    for ctr2 = 1:length(c_imgTrax)
        tempIdx      = find(trackList{c_imgTrax(ctr2),1}(:,1)==ctr1);
        c_xy(ctr2,:) = [trackList{c_imgTrax(ctr2),1}(tempIdx,:)  c_imgTrax(ctr2)];
        c_xy(ctr2,2) = c_xy(ctr2,2)/scalingFactor;
        c_xy(ctr2,3) = c_xy(ctr2,3)/scalingFactor;
    end
    
    
    [xx, yy]  = meshgrid(1:512); %For circular mask
    totalmask = zeros(size(c_img));
    
    figure, imagesc(c_img)
    hold on
    
    temp_xy = [c_xy(:,2:3) c_xy(:,5)];  %Passing only XY and TrackID
    [lbl badtrax] = locwater(c_img, temp_xy);
    
    %% local thresholding:
    %make a solid block of 1's 100 x 100 px around xy point
    pxwid = 100; %Must be an even number
    
    for ctr3 = 1:length(c_imgTrax)
        
        %points of the square:
        sx1 = round(c_xy(ctr3,2)- pxwid/2);
        sx2 = round(c_xy(ctr3,2)+ pxwid/2);
        
        sy1 = round(c_xy(ctr3,3)- pxwid/2);
        sy2 = round(c_xy(ctr3,3)+ pxwid/2);
        
        %Error handling for points with index>size(image)
        if sx1<1
            sx1 = 1;
        end
        if sx2>size(c_img,1)
            sx2 = size(c_img,1);
        end
        if sy1<1
            sy1 = 1;
        end
        if sy2>size(c_img, 2)
            sy2 = size(c_img,2);
        end
        
        %Square Mask:
        t_sqmask = zeros(size(c_img));
        t_sqmask(sy1:sy2, sx1:sx2) = 1;
        
        %Clean up image for multithresh:
        t_crop = t_sqmask.*double(c_img);
        %t_cropop = imopen(t_crop, strel('disk', 5));
        gfilt = fspecial('disk', 3);
        t_cropop = imfilter(t_crop, gfilt, 'replicate');
        mthresh = multithresh(t_cropop,5);
        mquant = imquantize(t_cropop, mthresh);
        %figure, imagesc(mquant);
        %Erode and reconstruct image
        
        x = imopen(t_crop, strel('disk', 10));
        
        figure, imagesc(x)
        
        
        
        
        %Circle around XY -> default segmenting when thresh is poor
        x = c_xy(ctr3, 2);
        y = c_xy(ctr3, 3);
        circ = sqrt((xx-x).^2 + (yy-y).^2)<=5; %Draw a circle with a radius of 5 centered around the current xy coords
        qfinder = nonzeros(mquant.*circ);
        modecheck = mode(qfinder);  %Find modal Multithresh values
        t_mask = zeros(size(c_img));
        t_maskfin = t_mask;
        
        %% Check Quality of Thresh
        %If at least 80pc. of the circ ROI is a single mode then thresh is
        %good, else its poor
        if length(find(qfinder==modecheck))>0.8*length(qfinder)
            t_mask(mquant==modecheck) = 1;
            
            %Clean Up:
            t_mask = imerode(t_mask, strel('disk', 3));
            t_mask = imfill(t_mask, 'holes');
            t_maskl = bwlabel(t_mask);
            modecheck2 = mode(nonzeros(circ.*t_maskl));
            t_maskfin(t_maskl==modecheck2) = 1;
            
        else
            t_maskfin = circ;
        end
        % figure, imagesc(t_mask)
        totalmask(t_maskfin==1) = ctr3;
        
        %% Plot Circle centered around current xy on any 'hold' image
        r = 4;
        thet = 0:pi/50:2*pi;
        xunit = r * cos(thet) + x;
        yunit = r * sin(thet) + y;
        plot(xunit, yunit, 'LineWidth',2)
        
    end
    
    figure, imagesc(totalmask)
end
