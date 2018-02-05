%local watershed

function varargout = locwater(image, locxy)
generic = 0; %If generic on, it can do something else with outputs.

im = mat2gray(image);
locxy = round(locxy);
im = imtophat(im, strel('disk', 60));                                        %correct uneven background

im_d1 = size(im,1);
im_d2 = size(im,2);

%Get XY Points for tracks in image:
idxlist = sub2ind([im_d1 im_d2], locxy(:,2), locxy(:,1));
xypts = zeros(im_d1,im_d2);
xypts(idxlist) = 1;
xypts = imdilate(xypts, strel('disk', 2));


lblfinal = zeros(im_d1, im_d2);

pxwid = 100;
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
        
        %Clean up image for multithresh:
        im_crop = t_sqmask.*double(im);


xy_crop= zeros(im_d1, im_d2);
xy_crop(idxlist(ctr3))=1;


%Local Thresholding:
im_iobr = imreconstruct(xy_crop, im_crop);
im_thresh = im2bw(im_iobr, graythresh(nonzeros(im_crop)));
im_fil = imfill(im_thresh, 'holes');
%figure, imagesc(im_fil)


%Running imerode in succession removes thinly connected objects (replace c.
%bwmorph)
im_mask = imerode(im_fil, strel('disk', 1));
im_mask = imerode(im_mask, strel('disk', 1));
im_mask = imerode(im_mask, strel('disk', 1));
im_mask2 = imdilate(im_mask, strel('disk', 6));


%Watershedding:
im_mask3 = im_mask2.*im;
dist = bwdist(~im_mask3);                                                    %Distance transform for border lines
impmin = imimposemin(dist, im_mask | bwperim(im_mask2));                     %Impose min on borders and foreground
ww = watershed(impmin);


%Catch XY coords just outside the segments:
ww_dil = imdilate(ww, strel('disk', 3));
figure, imagesc(ww_dil)

%only retain points of intersection
intersct = ww_dil(idxlist(ctr3));
ww_dil(ww_dil~=intersct)=0;
ww_dil(ww_dil==intersct)=locxy(ctr3,3);      %label the ww_dil according to the track list

lblfinal = lblfinal + double(ww_dil);

end  %LOOP END HERE



%% Resegmenter if two xypts are in the same watershed basin:
resegmenter = sort(intersct);
resegmenter = double(resegmenter(diff(resegmenter)==0));

%If more than three spots in a mask: Error out
if ~isempty(double(resegmenter(diff(resegmenter)==0)))  %If more than two spots accidentally fall into the same bin
    keyboard
end

%If two spots in a mask:
if ~isempty(resegmenter)
  reseg_temp = zeros(im_d1, im_d2);
  lblfinal = resegtool(reseg_temp, resegmenter, ww, lblfinal, idxlist);
end

keyboard
%% Area Filter
ww_rp = regionprops(lblfinal, 'Area');
area_filt = median([ww_rp(nonzeros(intersct)).Area]);                          %Area median of intersection between watershed label and xy points list
arfilt_H = area_filt*3;
arfilt_L = area_filt/4;


if generic ==0
    
   
    keyboard
    %This code fails: its not properly reading the trackIDs from the lbl
    %and therefore returning only ONE nuclues per mask. For example, run it
    %on ctr1 = 42. This will also give an example of an instance where you
    %lose one point and therefore have to force segmentation. 
    
%     So TBD: 
%     1. Fix below code to properly change this 
%     2. force the segmentation according to a circulat ROI
%     
    
    
    
    
    %Finding bad and good points:
    badpts_mask_all = find(([ww_rp(:).Area]>arfilt_H) | ([ww_rp(:).Area]<arfilt_L));
    
   
    
    
    badpts_locxy_idx = ismember(intersct, badpts_mask_all);
    goodpts_locxy_all = intersct(badpts_locxy_idx==0);
    goodpts_locxy_idx = double(intersct).*~badpts_locxy_idx;
    locxy(:,4) = goodpts_locxy_idx;
        
    
    %Outputs:
    lblfinal(~ismember(lblfinal, goodpts_locxy_all)) = 0;
    goodtrax = nonzeros(~badpts_locxy_idx.*locxy(:,3));
    badtrax = nonzeros(badpts_locxy_idx.*locxy(:,3));

    %Changing lbls to reflect TrackIDs
    %cheating and using a for loop:
    lblfinalcpy = lblfinal;
   % ww_dil = ww_dil_cpy;
    for aa = 1:size(locxy,1)
        lblfinal(lblfinalcpy==locxy(aa,4)) = locxy(aa,3);
    end
    
end


varargout{1} = lblfinal;
varargout{2} = badtrax;

end