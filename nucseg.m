%local watershed

function varargout = nucseg(image, locxy)
generic = 1; %If generic on, it can do something else with outputs.
badtrax = 0; %only if generic==1

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


%% Simple Watershed:

%Foreground Marker:
fgm = imdilate(xyimg, strel('disk', 3));

%Image Mask Perim. for Background Marker:
mask = im2bw(im, graythresh(im));
mask = imfill(mask, 'holes');
eroded = imerode(mask, strel('disk', 3));        %remove touching objects
erode2 = bwmorph(eroded, 'open', 6);             %remove touching objects
erode2(fgm>0)=1;                                 %Force xyimg points to remain:
dilated = imdilate(erode2, strel('disk', 10));   %ridge lines should be 'outside' the object to segment
bgm = bwperim(dilated);
bgm = ~bgm;

%Sobel filter for GradMag
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(im, hy, 'replicate');
Ix = imfilter(im, hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);

%Watershed:
gradmag2 = imimposemin(gradmag, ~bgm | fgm);
ww = watershed(gradmag2);

mask_seglist = unique(ww(ww&fgm));
track_seglist = unique(fgm(ww&fgm));

for ctr4 = 1:length(track_seglist)
    ctrack = track_seglist(ctr4);
    t_img = zeros(im_d1, im_d2);
    t_img(xyimg==ctrack)=ctrack;
    t_intersect = mode(ww(ww&t_img));
    lblfinal(ww==t_intersect)=ctrack;
    clear t_img 
end


%% Force circular mask for poor segmentation, find bad and good points.
if generic ==0
 %% Area Filter
ww_rp = regionprops(lblfinal, 'Area');
allcurrenttrack = unique(nonzeros(lblfinal));
allcurrenttrack = allcurrenttrack(allcurrenttrack<lbloffset);


area_filt = median([ww_rp(allcurrenttrack).Area]);                          %Area median of intersection between watershed label and xy points list
arfilt_H = area_filt*3;
arfilt_L = area_filt/5;
   
    missing_idx = ~ismember(locxy(:,3), allcurrenttrack);
    missing_track = locxy(missing_idx, 3);
    
    
    filt_track = find(([ww_rp(:).Area]>arfilt_H) | ([ww_rp(:).Area]<arfilt_L));
    %%
    
    
    
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

