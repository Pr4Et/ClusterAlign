%Calculate alignment error based on correlation between each tilt image and reproduced projection
%Load mrc file and reconstruct tomograms
%Written by Shahar Seifer, Weizmann Institute of Science, 2021
%Requires ASTRA Toolbox: https://www.astra-toolbox.com/downloads/index.html
%Requires @MRCImage library from MatTomo, PEET project: https://bio3d.colorado.edu/imod/matlab.html
clear;
gpu = gpuDevice();
reset(gpu);
disp('Error calucation based on back-projection reconstruction and then reprojection and comparison');
hoppe=input('Adaptive aspect ratio? (Hoppe=1, otherwise 0): ');
phi=input('Enter phi angle in degrees (0: y-axis, 90: x-axis rotation): ')*pi/180;
psi=input('Enter psi angle in degrees (0: default): ')*pi/180;
if isempty(psi)
psi=0;
end
if isempty(phi)
phi=90;
end
rotation_xaxis=(abs(cos(phi))<0.7);

[filename,path1] = uigetfile('D:\results\*.txt','Fetch fit error file');
Chosen_Filename_fiterr=[path1 filename];
Afiterr = readmatrix(Chosen_Filename_fiterr,'delimiter','\t');
Afiterr_last=Afiterr(length(Afiterr(:,1)),:);
fiterr=Afiterr_last;

bin=1;
nZ=400;
shift_limit=150;
do_filt=1;
save_rec=0;% use astra_reconst.m to SIRT reconstruction. input('Save reconstruction? =1');
max_fiterr=15;

[filename,path] = uigetfile('D:\results\*.mrc','Fetch ALIGNED MRC file');
Chosen_Filename=[path filename];
flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
showHeader=1; %  If 1 - Print out information as the header is loaded.
mRCImage=MRCImage;%Instentiate MRCImage in mRCImage
mRCImage = open(mRCImage, Chosen_Filename, flgLoadVolume, showHeader);
tilt = getVolume(mRCImage, [], [], []);
ntilts = getNZ(mRCImage);
nX = getNX(mRCImage);
nY = getNY(mRCImage);
sizeXangstrom=getCellX(mRCImage);
sizeYangstrom=sizeXangstrom;%getCellY(mRCImage) not used since it varies, but the actual is of the center angle 0
[filename_angles,path1] = uigetfile('D:\results\*.rawtlt','Fetch angles file');
Chosen_Filename_angles=[path1 filename_angles];
%angles=readmatrix(Chosen_Filename_angles,'FileType','text');
fileID = fopen(Chosen_Filename_angles,'r');
formatSpec = '%g';
angles = fscanf(fileID,formatSpec);
fclose(fileID);

tilt_rec=zeros(size(tilt));
%remove slices that fiterr was not between 0 and 5
count=0;
for idx=1:min([length(angles) ntilts])
    if fiterr(idx)>=0 && fiterr(idx)<=max_fiterr
        count=count+1;
        angles_rec(count)=angles(idx);
        tilt_rec(:,:,count)=tilt(:,:,idx);
    end
end



%Note it could crash if you lack a good GPU 
%det_spacing_x=1; %assuming rotation axis
%det_spacing_y=1;%although in raw data pixels scale as cos(angles);
det_row_count=round(nY/bin);
det_col_count=round(nX/bin);
% Create a 3D parallel beam geometry.  See the API for more information.
% det_spacing_x: distance between two horizontally adjacent detectors
% det_spacing_y: distance between two vertically adjacent detectors
% det_row_count: number of detector rows in a single projection
% det_col_count: number of detector columns in a single projection
% angles: projection angles in radians, should be between -pi/4 and 7pi/4
% proj_geom: MATLAB struct containing all information of the geometry
sizeZangstrom=nZ*sizeXangstrom/nX;
vol_geom = astra_create_vol_geom(nY/bin,nX/bin,nZ);
rec_id = astra_mex_data3d('create', '-vol', vol_geom, 0);
proj_vectors=zeros(length(angles_rec),12);
for idx=1:length(angles_rec)
    %assume the sample is static in axes and the entire setup rotates around it
    %In the case of rotation around x-axis (phi=90 degrees) and psi=0 the rays are described as:  
    %rayX=0; 
    %rayY=sin(angles(idx)*pi/180);
    %rayZ=cos(angles(idx)*pi/180);
    %In general the rays are A(theta)*(0,0,1),  where A is the rotation
    %transformation around axis (cos(psi)sin(phi),cos(psi)cos(phi),sin(psi))
    rayX=(1-cos(angles_rec(idx)*pi/180))*cos(psi)*sin(psi)*sin(-phi)+sin(angles_rec(idx)*pi/180)*cos(psi)*cos(-phi);
    rayY=(1-cos(angles_rec(idx)*pi/180))*cos(psi)*sin(psi)*cos(-phi)-sin(angles_rec(idx)*pi/180)*cos(psi)*sin(-phi);
    rayZ=cos(angles_rec(idx)*pi/180)+(1-cos(angles_rec(idx)*pi/180))*(sin(psi))^2;

    dX=0;
    dY=0;
    dZ=0;
    vX=0;    %u is for row shift of one pixel in the detector actual:(0,1)
    if hoppe
        vY=cos(angles_rec(idx)*pi/180)*cos(angles_rec(idx)*pi/180);
        vZ=-sin(angles_rec(idx)*pi/180)*cos(angles_rec(idx)*pi/180);
    else
        vY=cos(angles_rec(idx)*pi/180);
        vZ=-sin(angles_rec(idx)*pi/180);
    end
    uX=1;    %v is for column shift of one pixel in the detector actual (1,0)
    uY=0;
    uZ=0;
    proj_vectors(idx,:)=[rayX rayY rayZ dX dY dZ uX uY uZ vX vY vZ];
end
proj_geom = astra_create_proj_geom('parallel3d_vec',  det_row_count, det_col_count, proj_vectors);
%
% Create a 3D parallel beam geometry specified by 3D vectors.
%   See the API for more information.
% det_row_count: number of detector rows in a single projection
% det_col_count: number of detector columns in a single projection
% vectors: a matrix containing the actual geometry. Each row corresponds
%          to a single projection, and consists of:
%          ( rayX, rayY, rayZ, dX, dY, dZ, uX, uY, uZ, vX, vY, vZ )
%          ray: the ray direction
%          d  : the center of the detector 
%          u  : the vector from detector pixel (0,0) to (1,0)  -this is for x axis shift (correction in relation to manual, but consistent with order of u and v)
%          v  : the vector from detector pixel (0,0) to (0,1)  -The order in the imag/mrc file is (column,row), so this is for primary y' axis 
%          
% proj_geom: MATLAB struct containing all information of the geometry


proj_data_mat=zeros(nX/bin,length(angles_rec),nY/bin);
mask_data_mat=zeros(nX/bin,length(angles_rec),nY/bin);
%compensate for acquisition in adaptive aspect ratio
for idx=1:length(angles_rec)
    numrows=round(nY/bin);
    numcols=round(nX/bin);
    imag=double(imresize(tilt_rec(:,:,idx),[numcols numrows])); 
    simag1=mean(imag,1);%average of every line (parallel to x axis, the rotation), so it is easy to detect margins from alignments
    del_lines=simag1(2:numrows)-simag1(1:numrows-1);
    del_lines_usual=median(abs(del_lines));
    simag1c=mean(imag,2);%average of every col
    del_linesc=simag1c(2:numcols)-simag1(1:numcols-1);
    del_lines_usualc=median(abs(del_linesc));
    vect=3:round(numrows/4);
    if max(abs(del_lines(vect)))>5*del_lines_usual
        margin_row1=min(vect(abs(del_lines(vect))==max(abs(del_lines(vect)))))+1;
        temp1=double(imag(:,2:margin_row1-1));
        temp2=double(imag(:,margin_row1:numrows));
        if std(temp1(:))>0.25*std(temp2(:))
            margin_row1=1;
        end
    else
        margin_row1=1;
    end
    vect=round(3*numrows/4):numrows-2;
    if max(abs(del_lines(vect)))>5*del_lines_usual
        margin_row2=max(vect(abs(del_lines(vect))==max(abs(del_lines(vect)))));
        temp1=double(imag(:,margin_row2+2:numrows));
        temp2=double(imag(:,margin_row1:margin_row2));
        if std(temp1(:))>0.25*std(temp2(:))
            margin_row2=numrows;
        end
    else
        margin_row2=numrows;
    end

    vectc=3:round(numcols/4);
    if max(abs(del_linesc(vectc)))>5*del_lines_usualc
        margin_col1=min(vectc(abs(del_linesc(vectc))==max(abs(del_linesc(vectc)))))+1;
        temp1=double(imag(2:margin_col1-1,:));
        temp2=double(imag(margin_col1:numcols,:));
        if std(temp1(:))>0.25*std(temp2(:))
            margin_col1=1;
        end
    else
        margin_col1=1;
    end
    vectc=round(3*numcols/4):numcols-2;
    if max(abs(del_linesc(vectc)))>5*del_lines_usualc
        margin_col2=max(vectc(abs(del_linesc(vectc))==max(abs(del_linesc(vectc)))));
        temp1=double(imag(margin_col2+2:numcols,:));
        temp2=double(imag(margin_col1:margin_col2,:));
        if std(temp1(:))>0.25*std(temp2(:))
            margin_col2=numcols;
        end
    else
        margin_col2=numcols;
    end

    %normalize image intensity (this is bright field)
    croped_imag=imag(margin_col1:margin_col2,margin_row1:margin_row2);
    equal_imag=croped_imag;
    equal_imag=equal_imag-imgaussfilt(equal_imag,50);
    equal_imag=(10000.0*(equal_imag-mean(equal_imag(:)))/std(equal_imag(:)));
    imag=zeros(size(imag));
    imag(margin_col1:margin_col2,margin_row1:margin_row2)=equal_imag;
    imfull=imag;
    if ~hoppe && rotation_xaxis
    margin_up=max(margin_row1,1+round(round(numrows/2)*(1-cos(angles_rec(idx)*pi/180))));
    margin_down=min(margin_row2,round(round(numrows/2)*(1+cos(angles_rec(idx)*pi/180))));
    else
    margin_up=margin_row1;
    margin_down=margin_row2;
    end

    if ~hoppe && ~rotation_xaxis
    margin_left=max(margin_col1,1+round(round(numcols/2)*(1-cos(angles_rec(idx)*pi/180))));
    margin_right=min(margin_col2,round(round(numcols/2)*(1+cos(angles_rec(idx)*pi/180))));
    else
    margin_left=margin_col1;
    margin_right=margin_col2;
    end

    %Tapping unimportant regions
    if margin_up>1
        for n=1:margin_up
            imfull(:,n)=imag(:,n)*double(n/margin_up);
        end
    end
    if margin_down<numrows
        for n=margin_down:numrows
            imfull(:,n)=imag(:,n)*double((numrows-n)/(numrows-margin_down));
        end
    end
    if margin_left>1
        for n=1:margin_left
            imfull(n,:)=imag(n,:)*double(n/margin_left);
        end
    end
    if margin_right<numcols
        for n=margin_right:numcols
            imfull(n,:)=imag(n,:)*double((numcols-n)/(numcols-margin_right));
        end
    end
    figure(1)
    balanced_imshow(imfull);
    mask=zeros(size(imag));
    mask(:,margin_up:margin_down)=1;
    proj_data_mat(:,idx,:)=permute(imfull,[1 3 2]); % order should be: column(=x), angle, rows=y
    mask_data_mat(:,idx,:)=permute(mask,[1 3 2]); % order should be: column(=x), angle, rows=y
end
proj_id = astra_mex_data3d('create', '-proj3d', proj_geom, proj_data_mat);
%mask_id = astra_mex_data3d('create', '-proj3d', proj_geom, mask_data_mat);
%Then store back the manipolated according to the id of singoram

%cfg = astra_struct('SIRT3D_CUDA');
%cfg = astra_struct('CGLS3D_CUDA');
cfg = astra_struct('BP3D_CUDA');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = proj_id;
%cfg.option.SinogramMaskId=mask_id; %only with sirt3d and helps only at the margins of the volume slices
% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);
% Run 150 iterations of the algorithm
if length(cfg.type)==11
    if min(cfg.type=='SIRT3D_CUDA') || min(cfg.type=='CGLS3D_CUDA')
        astra_mex_algorithm('iterate', alg_id,150);
    end
end
if length(cfg.type)==9
    if min(cfg.type=='BP3D_CUDA')
        astra_mex_algorithm('run', alg_id);
    end
end
% Get the result
rec = astra_mex_data3d('get', rec_id);%maybe 'get_single'
if length(cfg.type)==11
    if min(cfg.type=='SIRT3D_CUDA')
        errorP=astra_mex_algorithm('get_res_norm', alg_id)/((nX/bin)*(nY/bin)*length(angles_rec));
        fprintf('Difference in ASTRA between source projections and reconstructions: %d \n',errorP);
    end
end


% Reconstruct projection data from reconstructed volume
if 1==1
for n=1:round(nZ/6)
    rec(:,:,n)=rec(:,:,n)*(n/(nZ/6));
    rec(:,:,nZ+1-n)=rec(:,:,nZ+1-n)*(n/(nZ/6));
end
end

proj_vectors2=zeros(length(angles),12);
proj_data_mat2=zeros(nX/bin,length(angles),nY/bin);
for idx=1:length(angles)
    rayX=(1-cos(angles(idx)*pi/180))*cos(psi)*sin(psi)*sin(-phi)+sin(angles(idx)*pi/180)*cos(psi)*cos(-phi);
    rayY=(1-cos(angles(idx)*pi/180))*cos(psi)*sin(psi)*cos(-phi)-sin(angles(idx)*pi/180)*cos(psi)*sin(-phi);
    rayZ=cos(angles(idx)*pi/180)+(1-cos(angles(idx)*pi/180))*(sin(psi))^2;

    dX=0;
    dY=0;
    dZ=0;
    vX=0;    %u is for row shift of one pixel in the detector actual:(0,1)
    if hoppe
        vY=cos(angles(idx)*pi/180)*cos(angles(idx)*pi/180);
        vZ=-sin(angles(idx)*pi/180)*cos(angles(idx)*pi/180);
    else
        vY=cos(angles(idx)*pi/180);
        vZ=-sin(angles(idx)*pi/180);
    end
    uX=1;    %v is for column shift of one pixel in the detector actual (1,0)
    uY=0;
    uZ=0;
    proj_vectors2(idx,:)=[rayX rayY rayZ dX dY dZ uX uY uZ vX vY vZ];

    imag=double(imresize(tilt(:,:,idx),[numcols numrows])); 
    proj_data_mat2(:,idx,:)=permute(imag,[1 3 2]); % order should be: column(=x), angle, rows=y

end
proj_geom2 = astra_create_proj_geom('parallel3d_vec',  det_row_count, det_col_count, proj_vectors2);


[proj_id, proj_data_fromrec] = astra_create_sino3d_cuda(rec, proj_geom2, vol_geom);

err_vector=zeros(1,length(angles));
for idx=1:length(angles)
    Imagem=permute(proj_data_fromrec(:,idx,:),[1 3 2]);
    Imagen=permute(proj_data_mat2(:,idx,:),[1 3 2]);
    shift_vect=r_mn(Imagem,Imagen,shift_limit,do_filt);
    err=sqrt(sum(shift_vect.^2));
    err_vector(idx)=err;
    fprintf('tilt slice=%d  Shift=%g [pix]\n',idx,err)
    pause(1);
end

if contains(Chosen_Filename,'jali.')
    if strcmp(strrep(Chosen_Filename,'jali.mrc',''),Chosen_Filename)
        Chosen_Filename_base=strrep(Chosen_Filename,'.mrc','.ex');
    else
        Chosen_Filename_base=strrep(Chosen_Filename,'jali.mrc','');
    end
else
    if strcmp(strrep(Chosen_Filename,'ali.mrc',''),Chosen_Filename)
        Chosen_Filename_base=strrep(Chosen_Filename,'.mrc','.ex');
    else
        Chosen_Filename_base=strrep(Chosen_Filename,'ali.mrc','');
    end
end

writematrix(err_vector,[Chosen_Filename_base '.clusteralign_err.csv']) 

fileID = fopen(Chosen_Filename_angles,'r');
formatSpec = '%g';
angles = fscanf(fileID,formatSpec);
fclose(fileID);

plotfile=strrep(Chosen_Filename_fiterr,'.txt','.tif');
Chosen_Filename_reproj=[Chosen_Filename_base '.clusteralign_err.csv'];
Areproj = err_vector;%readmatrix(Chosen_Filename_reproj);

%Afiterr_last(Afiterr_last<=0)=nan;
close all;
figure(1);
plot(angles,Afiterr_last,'r*',angles,Areproj','gh');
%ylim([0,4.5]);
xlabel('Angle [deg]');
ylabel('Error [pix]');
%legend('Fit error','Projection shift');
print(gcf,plotfile,'-dtiff');

disp(sprintf('RMS Fit error per slice (+-15deg): %g',sqrt(mean(Afiterr_last(abs(angles')<=15 & ~isnan(Afiterr_last))).^2)));
disp(sprintf('RMS alignement error (+-15deg): %g',sqrt(mean(Areproj(abs(angles)<=15).^2))));
disp(sprintf('Min. alignement error (+-1deg): %g',mean(Areproj(abs(angles)<=1))));

%Save to new MRC names rec_...
if save_rec
    newFilename=[path 'BP_rec_' filename];
    newmRCImage = MRCImage;%Instentiate MRCImage object
    newmRCImage.filename=newFilename;
    newmRCImage = setVolume(newmRCImage, rec); %enter to newmRCImage, do statistics, and fill many details to the header
    newmRCImage.header.cellDimensionX = sizeXangstrom;
    newmRCImage.header.cellDimensionY = sizeYangstrom;
    newmRCImage.header.cellDimensionZ = sizeZangstrom;
    %MODE=2:  32-bit signed real
    save(newmRCImage, newFilename);
    close(newmRCImage);
    disp(sprintf('Saved to file: %s ',newFilename));
end
    astra_mex_data3d('clear');%free GPU memory

%####################################################
function r_mn=r_mn(Imagem,Imagen,shift_limit,do_filt)
    if do_filt==1
        Imagem=imgaussfilt(Imagem-imgaussfilt(Imagem,100),3);
        Imagen=imgaussfilt(Imagen-imgaussfilt(Imagen,100),3);
    end
    figure(2);
    subplot(1,2,1);
    balanced_imshow(Imagem(round(0.15*size(Imagem,1)):floor(0.85*size(Imagem,1)),round(0.15*size(Imagem,2)):floor(0.85*size(Imagem,2))));
    subplot(1,2,2);
    balanced_imshow(Imagen(round(0.15*size(Imagen,1)):floor(0.85*size(Imagen,1)),round(0.15*size(Imagen,2)):floor(0.85*size(Imagen,2))));

    tempx=floor(0.3*size(Imagem,1));  % x are the row number, y is the col number (as observed with balanced_imshow). The rows progress along the first ordinate in Imagem/n.
    tempy=floor(0.3*size(Imagem,2));
    tempux=size(Imagem,1)-tempx;%floor(0.85*size(Imagem,1));
    tempuy=size(Imagem,2)-tempy;%floor(0.7*size(Imagem,2));
    view_in=Imagem(tempx:tempux,tempy:tempuy);
    correlationOutput = normxcorr2(view_in,Imagen);
    [maxCorrValue, maxIndex] = max(abs(correlationOutput(:)));
    [xpeak, ypeak] = ind2sub(size(correlationOutput),maxIndex(1));%find(correlationOutput==max(correlationOutput(:)));  xpeak is the row number
    yoffset = ypeak-tempuy;
    xoffset = xpeak-tempux;
    if abs(yoffset)>shift_limit || abs(xoffset)>shift_limit
        correlationOutput = normxcorr2(imgaussfilt(view_in,10),imgaussfilt(Imagen,10));
        [maxCorrValue, maxIndex] = max(abs(correlationOutput(:)));
        [xpeak, ypeak] = ind2sub(size(correlationOutput),maxIndex(1));%find(correlationOutput==max(correlationOutput(:)));
        yoffset = ypeak-tempuy;
        xoffset = xpeak-tempux;
        if abs(yoffset)>shift_limit || abs(xoffset)>shift_limit
            r_mn=[NaN NaN];
        else
            r_mn=[xoffset yoffset];
        end
        disp('Only rough shift estimate')
        return;
    end
    %refine to subpixel
    sample16=correlationOutput(xpeak-7:xpeak+8,ypeak-7:ypeak+8);
    Intsample16=fftInterpolate(sample16,[512 512]);
    [maxCorrValue2, maxIndex2] = max(abs(Intsample16(:)));
    [xpeak2, ypeak2] = ind2sub(size(Intsample16),maxIndex2(1));%find(Intsample16==max(Intsample16(:)));
    yoffset2=yoffset+(ypeak2-256+30)/32;
    xoffset2=xoffset+(xpeak2-256+31)/32;
    r_mn=[xoffset2 yoffset2];
end


