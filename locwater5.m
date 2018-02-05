%local watershed

function varargout = locwater(image, locxy)
generic = 1; %If generic on, it can do something else with outputs.
badtrax = 0; %only if generic==1

localthresh = 'thresh';

im = mat2gray(image);
locxy = round(locxy);
im = imtophat(im, strel('disk', 60));                                        %correct uneven background
immax = max(im(:));

im_d1 = size(im,1);
im_d2 = size(im,2);

%Get XY Points for tracks in image:
idxlist = sub2ind([im_d1 im_d2], locxy(:,2), locxy(:,1));

%XY points image for resegmenter (Don't dilate this!):
xyimg = zeros(im_d1, im_d2);
xyimg(idxlist) = locxy(:,3);
%figure, imagesc(xyimg)

%Initialize final image:
lblfinal = ones(im_d1, im_d2);
lbloffset = 1000;
lblfinal = lblfinal.*lbloffset; %This prevents problems with TRACK_ID==0 and 1

%% Local Thresholding + Watershed:

pxwid = 100;    %Crop Window for local thresholding
for ctr3 = 1:size(locxy, 1);
    
    sx1 = round(locxy(ctr3,1)- pxwid/2);
    sx2 = round(locxy(ctr3,1)+ pxwid/2);
    
    sy1 = round(locxy(ctr3,2)- pxwid/2);
    sy2 = round(locxy(ctr3,2)+ pxwid/2);
    
    %Error handling for points with index>size(image)
    if sx1<1
        sx1 = 1;
    end
    if sx2>size(im,1)
        sx2 = size(im,1);
    end
    if sy1<1
        sy1 = 1;
    end
    if sy2>size(im, 2)
        sy2 = size(im,2);
    end
    
    %Square Mask:
    t_sqmask = zeros(size(im));
    t_sqmask(sy1:sy2, sx1:sx2) = 1;
    
    %Croppped Image:
    im_crop = t_sqmask.*double(im);
    
    %Cropped XY-point mask:
    xy_crop= zeros(im_d1, im_d2);
    xy_crop(idxlist(ctr3))=1;
    xy_crop = imdilate(xy_crop, strel('disk', 3));
    
    
    
    %Local Thresholding:
    
    switch localthresh

        %Reconstruction
        case('reconstruction')
            im_iobr = imreconstruct(xy_crop, im_crop);
            im_thresh = im2bw(im_iobr, graythresh(nonzeros(im_iobr))*1.2);
            
        %Multithresh:
        case('quantization')
            threshl = multithresh(nonzeros(im_crop),5);
            im_thresh = imquantize(im_crop, threshl);
            im_thresh = im2bw(im_thresh, 2);
            
        %Simple Thresh
        case('thresh')
            bw = im2bw(im_crop, graythresh(nonzeros(im_crop)*0.9));
            im_thresh = imfill(bw, 'holes');
            
    end
    
    %Break apart thinly connected objects:
    im_mask_erode = imerode(im_thresh, strel('disk', 5));
    im_fil = imfill(im_mask_erode, 'holes');
    im_mask_erode = bwmorph(im_fil, 'thin', 2);
    im_mask_erode = bwmorph(im_mask_erode, 'open', 6);
    
    %Force the eroded mask to include the xy point(ensures its retained during watershedding)
    im_mask_erode(imdilate(xy_crop, strel('disk', 3))==1)=1;
    
            
    %Retain only the mask of the current cell (Intersection between XY Point and Mask):
    [maskerode_lbl n] = bwlabel(im_mask_erode);
    mask_ctrack = unique(maskerode_lbl(xy_crop>0));
    maskerode_lbl(maskerode_lbl~=mask_ctrack)=0;

    %Dilate to form the outer watershed ridgeline
    im_mask_dil = imdilate(maskerode_lbl, strel('disk', 5));
    
    %Canny Filter to enhance ridgeline:
    canfilt = edge(im_crop, 'canny', 0.5).*20;
    
    %Watershedding:
    im_maskedIM = im_mask_dil.*im;
    im_maskedIM = im_maskedIM + canfilt; %Enhances boundaries of nuclei 
    dist = bwdist(~im_maskedIM);                                                    %Distance transform for border lines
    impmin = imimposemin(dist, im_mask_erode | bwperim(im_mask_dil));                     %Impose min on borders and foreground
    ww = watershed(impmin);
    ww = mat2gray(ww);

    %Catch XY coords just outside the segments:
    ww_dil = imdilate(ww, strel('disk', 2));

    %only retain points of intersection
    intersct = ww_dil(idxlist(ctr3));
    ww_dil(ww_dil~=intersct)=lbloffset;
    ww_dil(ww_dil==intersct)=locxy(ctr3,3);      %label the ww_dil according to the track list
    
    lblfinal(ww_dil==locxy(ctr3,3)) = locxy(ctr3,3);
   
end  %LOOP END HERE
% 
keyboard
tt = lblfinal;
tt(lblfinal==1000) = 0.1;
figure, imagesc(tt)


%% simple
n = mat2gray(im);
can = edge(n,'canny', 0.5).*20;
sob = edge(n,'Sobel', 0.5, 'both');

mask = im2bw(n, graythresh(n));
mask = imfill(mask, 'holes');

eroded = imerode(mask, strel('disk', 2));

    erode2 = bwmorph(eroded, 'open', 6);
    %Force xyimg points
        erode2(imdilate(xyimg, strel('disk', 3))>0)=1;
    
    dilated = imdilate(erode2, strel('disk', 10));
    xy2 = xyimg;
    xy2(xyimg>0)=1;
    xy2 = imdilate(xy2, strel('disk',3)).*20;
    maskedIM = dilated.*im;
    maskedIM = maskedIM + can + xy2; %Enhances boundaries of nuclei 
    dist = bwdist(~maskedIM);                                                    %Distance transform for border lines
    impmin = imimposemin(dist, erode2| bwperim(im_mask_dil));                     %Impose min on borders and foreground
    ww = watershed(impmin);
    ww = mat2gray(ww);
figure, imagesc(ww)
%%



%% Eliminate overlapped ROIs
%Still need to think about how to implement this. Ideally you want to run
%it within the loop, check if there are any masks already made in the
%places you want to mask, and then accordingly make adjustments. But, idk.
%figure, imagesc(lblfinal)
% 
% %% Resegmenter if two xypts are in the same watershed basin:
% resegmenter = lblfinal(xyimg>0);
% resegmenter = double(resegmenter(diff(resegmenter)==0));
% 
% %If two spots in a mask:
% if ~isempty(resegmenter)
%     lblcopy = lblfinal;
%     lblfinal = resegtool(im_d1, im_d2, xyimg, resegmenter, im, lblfinal, idxlist);
% end
% 
% 
% 



%% Area Filter
ww_rp = regionprops(lblfinal, 'Area');
allcurrenttrack = unique(nonzeros(lblfinal));
allcurrenttrack = allcurrenttrack(allcurrenttrack<lbloffset);


area_filt = median([ww_rp(allcurrenttrack).Area]);                          %Area median of intersection between watershed label and xy points list
arfilt_H = area_filt*3;
arfilt_L = area_filt/5;



%% Force circular mask for poor segmentation, find bad and good points.
if generic ==0
    
    missing_idx = ~ismember(locxy(:,3), allcurrenttrack);
    missing_track = locxy(missing_idx, 3);
    
    
    filt_track = find(([ww_rp(:).Area]>arfilt_H) | ([ww_rp(:).Area]<arfilt_L));
    
    
    
    
    %Fix this last step, filtering and forcing circle mask The problem is
    %with the assumption of equal length for filtered vex.
    
       
    
    %Finding bad and good points:
    badpts_mask_all = find(([ww_rp(:).Area]>arfilt_H) | ([ww_rp(:).Area]<arfilt_L));
    
    badpts_locxy_idx = ismember(allcurrenttrack, badpts_mask_all);
    goodpts_locxy_all = allcurrenttrack(badpts_locxy_idx==0);
    
    %Clean up all MISSING and FILTERED points (badpts):
    lblfinal(~ismember(lblfinal, goodpts_locxy_all)) = 0;
    
    %Force circular mask:
    badpts = nonzeros(badpts_locxy_idx.*allcurrenttrack);
    if ~(isempty(badpts))
        circ_mask = zeros(im_d1, im_d2);
        circ_idx = nonzeros(idxlist.*badpts_locxy_idx);
        
        for ctr4 = 1:length(circ_idx)
            circ_mask(circ_idx(ctr4)) = badpts(ctr4);
        end
        
        circ_mask = imdilate(circ_mask, strel('disk', 5));
        lblfinal(circ_mask>0) = circ_mask(circ_mask>0);
    end
    
    
    %List out useful tracks:
    finlist = nonzeros(unique(lblfinal));
    
    %This is rewriting the trackIDs, not normally necessary, but leaving it
    %in here for edge cases.
    locxy(:,4) = finlist;
    
    %Outputs:
    goodtrax = locxy(locxy(:,3)==locxy(:,4),4);
    badtrax = locxy(locxy(:,3)~=locxy(:,4),4);
    
    %Changing lbls to reflect TrackIDs
    %cheating and using a for loop:
    lblfinalcpy = lblfinal;
    for aa = 1:size(locxy,1)
        lblfinal(lblfinalcpy==locxy(aa,4)) = locxy(aa,3);
    end
    
end


varargout{1} = lblfinal;
varargout{2} = badtrax;

end

