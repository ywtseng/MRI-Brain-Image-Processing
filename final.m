%{
    When using, need:
    1. Job file (default file name: final_job.m)
    2. Folder contains DICOM image 
    3. Folder conatins functions (default folder name: functions)

    Result:
    1. seg_result folder:
        grey.stl
        white.stl
        other segmentation .nii files

    Need Image_Toolbox for: 
    1. Resize image (imreize)
    2. Find Largest Connected Component (getLargestCc)
    
    Need change manually:
    1. jobfile (variable): final_job.m  (file path)
    2. final_job.m (file): matlabbatch{4}.spm.spatial.preproc.opts.tpm (spm8 mask path)
    3. file number (variable): dir_num

    Customize:
    1. Cut: specify [cutx, cuty, cutz] in Part2: cut
    2. Resize: specify [res] in Part3: Resize

%}
close all
clear all
clc
a = license ('checkout','Image_Toolbox');

%-------------------------%
% Part 1: Segmentation %
%-------------------------%
% Convert DICOM into Nifti and segmentation
start = tic;
tic
spm('defaults', 'fmri');
spm_jobman('initcfg');
jobfile = {'C:\Users\LAB-690\Desktop\104project\project\final_job.m'};
spm_jobman('serial', jobfile, '');
fprintf('\nDone: Convert DICOM into Nifti and Segmentation\n')
toc


%----------------------%
% Part 2: Image Proc   %
%----------------------%
% Get output file
tic
tmp = dir;
list = {tmp.name};
dir_num = 3;
grey_matter_fname = list{dir_num};
white_matter_fname = list{dir_num+1};
csf_fname = list{dir_num+2};
clear tmp list
%----------------------%
%----------------------%
% Read nii 
gnii = load_nii(grey_matter_fname);
wnii = load_nii(white_matter_fname);
cnii = load_nii(csf_fname);
g_img = gnii.img;
w_img = wnii.img;
c_img = cnii.img;
img_num = gnii.hdr.dime.dim(2);

if (a)
    % resize
    pix_size =  gnii.hdr.dime.pixdim(2)/gnii.hdr.dime.pixdim(3);
    for i=1:img_num
        ptg(:,:,i) = imresize(squeeze(g_img(i,:,:)), pix_size);
        ptw(:,:,i) = imresize(squeeze(w_img(i,:,:)), pix_size);
        ptc(:,:,i) = imresize(squeeze(c_img(i,:,:)), pix_size);
    end
    g_img = permute(ptg, [3 1 2]);
    w_img = permute(ptw, [3 1 2]);
    c_img = permute(ptc, [3 1 2]);
end

% threshold
threshold = 0.5;

g_img(g_img<threshold) = 0;
w_img(w_img<threshold) = 0;
c_img(c_img<threshold) = 0;
% clear overlapping voxel
g_timg = g_img;
g_timg(g_img>w_img & g_img>c_img) = 1;
g_timg(g_img<=w_img & g_img<=c_img) = 0;
w_timg = w_img;
w_timg(w_img>=g_img & w_img>c_img & w_img ~= 0 ) = 1;
w_timg(w_img<g_img & w_img<=c_img) = 0;

% cut
dimx = size(g_timg, 3);
dimy = size(g_timg, 2);
dimz = size(g_timg, 1);
%{
       x dim: 1->max bottom to top
       y dim: 1->max back to front
       z dim: 1->196 right ear to left ear

       ori:     x=1, y=1, z=1
       1(cut a quarter):    x=dimx*2/3, y=dimy/2, z=dimz/2
       2(left brain):   x=dimx*2/3, y=1, z=dimz;
       3(upper brain):   x=dimx*11/18, y=1, z=dimz   
%}
cutx = round(1);
cuty = round(1);
cutz = round(1);

g_timg(1:cutz, cuty:dimy, cutx:dimx) = 0;
w_timg(1:cutz, cuty:dimy, cutx:dimx) = 0;
fprintf('\nDone: Threshold=0.5 and Cut\n')
toc
%----------------------%
%----------------------%
% Find Largest Connected component (Need Image_Toolbox)
if (a)
    tic
    g_timg = logical(g_timg);
    [g_timg rp] = getLargestCc(g_timg,[],1);
    w_timg = logical(w_timg);
    [w_timg rp] = getLargestCc(w_timg,[],1);
    fprintf('\nDone: Find Largest Connected component\n')
    toc
end
%----------------------%
%----------------------%
%{
% Show .nii
tg = gnii;
tw = wnii;

if(a)
    tg.hdr.dime.pixdim(3) = 1;
    tg.hdr.dime.pixdim(4) = 1;
    tw.hdr.dime.pixdim(3) = 1;
    tw.hdr.dime.pixdim(4) = 1;
end
tg.img = g_timg;
tw.img = w_timg;
view_nii(tg);
view_nii(tw);
view_nii(gnii);
view_nii(wnii);

% Save .nii
save_nii(tg, 'stl_grey_matter.nii');
save_nii(tw, 'stl_white_matter.nii');
%}


%----------------------%
% Part 3: making STL  %
%----------------------%
% Resize to fit machine
res = 0.7;
dimx = round(size(g_timg, 1)*res);
dimy = round(size(g_timg, 2)*res);
dimz = round(size(g_timg, 3)*res);

g_timg = resize(g_timg, [dimx dimy dimz]);
w_timg = resize(w_timg, [dimx dimy dimz]);

% Make .stl
tic
gfv = isosurface(~g_timg, 0.5);
wfv = isosurface(~w_timg, 0.5);
fprintf('\nDone: Find isosurface\n')
toc
%----------------------%
%----------------------%
tic
[sgfv.faces, sgfv.vertices] = N_SmoothMesh(gfv.faces, gfv.vertices);
[swfv.faces, swfv.vertices] = N_SmoothMesh(wfv.faces, wfv.vertices);
fprintf('\nDone: Smooth isosurface\n')
toc

tic
stlwrite('grey_ori_resize.stl', sgfv);
stlwrite('white_ori_resize.stl', swfv);
toc

elapsed = toc(start);
fprintf('\nTotal time: %f .\n', elapsed);

%-----------end----------%