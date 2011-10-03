function erode_mask_cereb_flr(in_file, erosionval)
erosionval = str2num(erosionval);
nii=load_untouch_nii(in_file);
% nii a structure, .img is the part of the structure where the data is
% stored
I=nii.img;

%squeeze empties the non used dimension of the array, we take out the third
%dimension
IM=squeeze(I(:,:,1));
%initalizing a matrix, copying to IMG1 for vertical scanning, IMG used for
%horizontal scanning. 
IMG=zeros(size(IM,1), size(IM,2)); IMG1 = IMG;
%for every slice
for z=1:42;
    %remove the z dimension
    IM=squeeze(I(:,:,z));
    for y=1:512
        %erosionval number of voxels to be removed around image
        %x for horizontal plane scanning
        x=squeeze(IM(y,:));
        %p for vertical plane scanning
        p=squeeze(IM(:,y));
        %k is all the voxels where x == 1
        k=find(x==1);
        k1=find(p==1);
        
        %if x is not 1; k is an array of all locations where x = 1, we are
        %scanning for that value. Only if k is not empty execute:
        if isempty(k) == 0
            startx=k(1); endx=k(end);
            %from the values startx to erosion threshold, set to 0, also on
            %the other side with endx
            x(1:startx+erosionval)=0; x(endx-erosionval:end)=0;
            %set x value in original IMG to changed vector x
            IMG(y,:) = x;
        else
            %if there are no ones, replace x with original value
            IMG(y,:) = x;
        end
        
        %the same operation in vertical plane, replace all of the rows with
        %the eroded threshold otherwise replace with original.  
        if isempty(k1) == 0
            startx=k1(1); endx=k1(end);
            p(1:startx+erosionval)=0; p(endx-erosionval:end)=0;
            IMG1(:,y) = p;
        else
            IMG1(:,y) = p;
        end        
    end
    %clever!
    IM = min(IMG,IMG1);
    I(:,:,z) = IM;
end
%setting .img part of the structure with I
nii.img = I;
%save within the same directory with filename eroded_image
save_untouch_nii(nii, 'out_mask.nii');



%figure(1); imagesc(I); colormap gray; axis image;
exit;
 
      
