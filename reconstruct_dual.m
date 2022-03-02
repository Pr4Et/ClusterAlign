function reconstruct_dual(cosine_sample,phi_deg1,psi_deg1,Chosen_Filename_ali1,Chosen_Filename_angles1,phi_deg2,psi_deg2,Chosen_Filename_ali2,Chosen_Filename_angles2,varargin)
if false
    cosine_sample=0;
    phi_deg1=90;
    psi_deg1=0;
    phi_deg2=90;
    psi_deg2=0;
    Chosen_Filename_ali1='D:\results\AD68_1_pt12a.mrc';
    Chosen_Filename_ali2='D:\results\AD68_1_pt12b.mrc';
    Chosen_Filename_angles1='D:\results\AD68_1_pt12.rawtlt';
    Chosen_Filename_angles2='D:\results\AD68_1_pt12.rawtlt';
    varargin='';
end
%(C) Written by Shahar Seifer, Weizmann Institute of Science, 2022
%Load aligned mrc tomograms and reconstruct 3D by SIRT3D in Astra toolbox -
%1 is for a and 2 is for b tilt series in dual axis tomogram
%Requires ASTRA Toolbox, https://www.astra-toolbox.com/downloads/index.html#
%Requires @MRCImage library from MatTomo, PEET project: https://bio3d.colorado.edu/imod/matlab.html
%varagin includes: Chosen_Filename_fiterr,bin,nZ
%##### check parameters ####
if ~isempty(varargin)
    Chosen_Filename_fiterr=varargin{1};
end
if length(varargin)>1
    bin=varargin{2};
    if bin>8
        disp(sprintf('Note that you have chosen bin=%g',bin));
    end
else
    bin=4; %binning (bin=2 is suitable for 8Gb GPU for output volume of 1024X1024X400 and input of 2048X2048 images)
end
if length(varargin)==3
    nZ=varargin{3};
    if nZ<10
        disp(sprintf('Note that you have chosen thickness of %g',nZ));
    end
else
    nZ=400; %thickness (in pixels after binning)  
end

max_fiterr=15;
%###########################
gpu = gpuDevice();
reset(gpu);
%cosine_sample=input('Adaptive aspect ratio? (cosine_sample=1, otherwise 0): ');
%phi=input('Enter phi angle in degrees (0: y-axis, 90: x-axis rotation): ')*pi/180;
%psi=input('Enter psi angle in degrees (0: default): ')*pi/180;
if isempty(psi_deg1)
psi_deg1=0;
end
if isempty(phi_deg1)
phi_deg1=90;
end
if isempty(psi_deg2)
psi_deg2=0;
end
if isempty(phi_deg2)
phi_deg2=0;
end

rotateb=0;
if abs(phi_deg2-phi_deg1)<30
    if phi_deg2>45
        phi_deg2=phi_deg2-90;
        rotateb=-90; % + is clockwise, so we rotate b image counterclockwise
    else
        phi_deg2=phi_deg2+90;
        rotateb=+90; % + is clockwise
    end
end

phi1=phi_deg1*pi/180;
psi1=psi_deg1*pi/180;
phi2=phi_deg2*pi/180;
psi2=psi_deg2*pi/180;

rotation_xaxis1=(abs(cos(phi1))<0.7); %rough estimate for tapering margins only
rotation_xaxis2=(abs(cos(phi2))<0.7); %rough estimate for tapering margins only


%[filename,path] = uigetfile('D:\results\*.mrc','Fetch MRC file');
%Chosen_Filename_ali=[path filename];
flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
showHeader=1; %  If 1 - Print out information as the header is loaded.
mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
mRCImage = open(mRCImage, Chosen_Filename_ali1, flgLoadVolume, showHeader);
tilt1 = double(getVolume(mRCImage, [], [], []));
ntilts1 = getNZ(mRCImage);
nX = getNX(mRCImage);
nY = getNY(mRCImage);
sizeXangstrom=getCellX(mRCImage);
sizeYangstrom=sizeXangstrom;%getCellY(mRCImage) not used since it varies, but the actual is of the center angle 0
%[filename_angles,path1] = uigetfile('D:\results\*.rawtlt','Fetch angles file');
%Chosen_Filename_angles=[path1 filename_angles];
fileID = fopen(Chosen_Filename_angles1,'r');
formatSpec = '%g';
angles1 = fscanf(fileID,formatSpec);
fclose(fileID);

mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
mRCImage = open(mRCImage, Chosen_Filename_ali2, flgLoadVolume, showHeader);
tilt2 = double(getVolume(mRCImage, [], [], []));
ntilts2 = getNZ(mRCImage);
nX = getNX(mRCImage);
nY = getNY(mRCImage);
sizeXangstrom=getCellX(mRCImage);
sizeYangstrom=sizeXangstrom;%getCellY(mRCImage) not used since it varies, but the actual is of the center angle 0
%[filename_angles,path1] = uigetfile('D:\results\*.rawtlt','Fetch angles file');
%Chosen_Filename_angles=[path1 filename_angles];
fileID = fopen(Chosen_Filename_angles2,'r');
formatSpec = '%g';
angles2 = fscanf(fileID,formatSpec);
fclose(fileID);

%Correct rotation
for n=1:ntilts2
    imag=tilt2(:,:,n);
    imag=imrotate(imag,rotateb,'bilinear','crop'); %rotate + or -90 deg  
    tilt2(:,:,n)=imag;
end
%align b with respect to a, based on center slices
vect=1:length(angles1);
ind1=min(vect(angles1==min(abs(angles1))));
vect=1:length(angles2);
ind2=min(vect(angles2==min(abs(angles2))));
imag1=tilt1(:,:,ind1);
imag2=tilt2(:,:,ind2);
shift_limit=500;
do_filt=1;
rshift=r_mn_rect(imag2,imag1,shift_limit,do_filt);
%figure(3)
%balanced_imshow(imag1);
%figure(4)
%balanced_imshow(imtranslate(imag2,rshift));
for n=1:ntilts2
    imag=tilt2(:,:,n);
    imag=imtranslate(imag,rshift); %align translational shift in b image  
    tilt2(:,:,n)=imag;
end



%try
%    Afiterr_all = readmatrix(Chosen_Filename_fiterr,'delimiter','\t');
%    fiterr=Afiterr_all(length(Afiterr_all(:,1)),1:min([length(angles) ntilts]));
%catch
fiterr1=zeros(1,length(angles1));
fiterr2=zeros(1,length(angles2));
%end


%remove slices that fiterr was not between 0 and 5
count=0;
for idx=1:min([length(angles1) ntilts1])
    if fiterr1(idx)>=0 && fiterr1(idx)<=max_fiterr
        count=count+1;
        angles1(count)=angles1(idx);
        tilt1(:,:,count)=tilt1(:,:,idx);
    end
end
angles1=angles1(1:count);
count1=count;
count=0;
for idx=1:min([length(angles2) ntilts2])
    if fiterr2(idx)>=0 && fiterr2(idx)<=max_fiterr
        count=count+1;
        angles2(count)=angles2(idx);
        tilt2(:,:,count)=tilt2(:,:,idx);
    end
end
angles2=angles2(1:count);
angles=zeros(1,length(angles1)+length(angles2));
angles(1:length(angles1))=angles1;
angles(1+length(angles1):length(angles1)+length(angles2))=angles2;

%Note it could crash if you lack a good GPU 
det_row_count=round(nY/bin);
det_col_count=round(nX/bin);
sizeZangstrom=nZ*sizeXangstrom/nX;
vol_geom = astra_create_vol_geom(nY/bin,nX/bin,nZ);
rec_id = astra_mex_data3d('create', '-vol', vol_geom, 0);
proj_vectors=zeros(length(angles1)+length(angles2),12);
for idx=1:length(angles)
    if idx<=count1
        rotation_xaxis=rotation_xaxis1;
        phi=phi1;
        psi=psi1;
    else
        rotation_xaxis=rotation_xaxis2;
        phi=phi2;
        psi=psi2;
    end
    %transformation around axis (cos(psi)sin(phi),cos(psi)cos(phi),sin(psi))
    rayX=(1-cos(angles(idx)*pi/180))*cos(psi)*sin(psi)*sin(-phi)+sin(angles(idx)*pi/180)*cos(psi)*cos(-phi);
    rayY=(1-cos(angles(idx)*pi/180))*cos(psi)*sin(psi)*cos(-phi)-sin(angles(idx)*pi/180)*cos(psi)*sin(-phi);
    rayZ=cos(angles(idx)*pi/180)+(1-cos(angles(idx)*pi/180))*(sin(psi))^2;
    dX=0;
    dY=0;
    dZ=0;
    if cosine_sample
        if rotation_xaxis
            vX=0;    %for row shift of one pixel in the detector actual:(0,1)
            vY=cos(angles(idx)*pi/180)*cos(angles(idx)*pi/180);
            vZ=-sin(angles(idx)*pi/180)*cos(angles(idx)*pi/180);
        else
            uY=0;    %for row shift of one pixel in the detector actual:(1,0)
            uX=cos(angles(idx)*pi/180)*cos(angles(idx)*pi/180);
            uZ=-sin(angles(idx)*pi/180)*cos(angles(idx)*pi/180);
        end
    else
        if rotation_xaxis
            vX=0;    %for row shift of one pixel in the detector actual:(0,1)
            vY=cos(angles(idx)*pi/180);
            vZ=-sin(angles(idx)*pi/180);
        else
            uX=cos(angles(idx)*pi/180);    
            uY=0;
            uZ=-sin(angles(idx)*pi/180);
        end
    end
    if rotation_xaxis
        uX=1;    %for column shift of one pixel in the detector actual (1,0)
        uY=0;
        uZ=0;
    else
        vY=1;
        vX=0;
        vZ=0;
    end
    proj_vectors(idx,:)=[rayX rayY rayZ dX dY dZ uX uY uZ vX vY vZ];
end
proj_geom = astra_create_proj_geom('parallel3d_vec',  det_row_count, det_col_count, proj_vectors);
%
% Create a 3D parallel beam geometry specified by 3D vectors.
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

proj_data_mat=zeros(nX/bin,length(angles),nY/bin);
%compensate for acquisition in adaptive aspect ratio
for idx=1:length(angles)
    %if cosine_sample
    %numrows=round((nY/bin)*cos(angles(idx)*pi/180));
    %else
    %numrows=round(nY/bin);
    %end
    if idx<=count1
        rotation_xaxis=rotation_xaxis1;
    else
        rotation_xaxis=rotation_xaxis2;
    end

    numrows=round(nY/bin);
    numcols=round(nX/bin);
    if idx<=count1
        imag=imresize(tilt1(:,:,idx),[numcols numrows]);
    else
        imag=imresize(tilt2(:,:,idx-count1),[numcols numrows]);
    end
    simag1=mean(imag(round(0.35*numcols):round(0.65*numcols),:),1);%average of every line (parallel to x axis, the rotation), so it is easy to detect margins from alignments
    del_lines=simag1(2:numrows)-simag1(1:numrows-1);
    del_lines_usual=median(abs(del_lines));
    simag1c=mean(imag(:,round(0.35*numrows):round(0.65*numrows)),2);%average of every col
    del_linesc=simag1c(2:numcols)-simag1c(1:numcols-1);
    del_lines_usualc=median(abs(del_linesc));
    vect=3:round(0.4*numrows);
    if max(abs(del_lines(vect)))>10*del_lines_usual
        margin_row1=min(vect(abs(del_lines(vect))==max(abs(del_lines(vect)))))+1;
        temp1=double(imag(round(0.35*numcols):round(0.65*numcols),2:margin_row1-2));
        temp2=double(imag(round(0.35*numcols):round(0.65*numcols),margin_row1:numrows));
        if std(temp1(:))>0.25*std(temp2(:))
            margin_row1=1;
        else 
            margin_row1=margin_row1+1;
        end
    else
        margin_row1=1;
    end
    vect=round(0.6*numrows):numrows-2;
    if max(abs(del_lines(vect)))>10*del_lines_usual
        margin_row2=max(vect(abs(del_lines(vect))==max(abs(del_lines(vect)))));
        temp1=double(imag(round(0.35*numcols):round(0.65*numcols),margin_row2+2:numrows));
        temp2=double(imag(round(0.35*numcols):round(0.65*numcols),margin_row1:margin_row2));
        if std(temp1(:))>0.25*std(temp2(:))
            margin_row2=numrows;
        else
            margin_row2=margin_row2-1;
        end
    else
        margin_row2=numrows;
    end
    vectc=3:round(0.4*numcols);
    if max(abs(del_linesc(vectc)))>10*del_lines_usualc
        margin_col1=min(vectc(abs(del_linesc(vectc))==max(abs(del_linesc(vectc)))))+1;
        temp1=double(imag(2:margin_col1-2,round(0.35*numrows):round(0.65*numrows)));
        temp2=double(imag(margin_col1:numcols,round(0.35*numrows):round(0.65*numrows)));
        if std(temp1(:))>0.25*std(temp2(:))
            margin_col1=1;
        else
            margin_col1=margin_col1+1;
        end
    else
        margin_col1=1;
    end
    vectc=round(0.6*numcols):numcols-2;
    if max(abs(del_linesc(vectc)))>10*del_lines_usualc
        margin_col2=max(vectc(abs(del_linesc(vectc))==max(abs(del_linesc(vectc)))));
        temp1=double(imag(margin_col2+2:numcols,round(0.35*numrows):round(0.65*numrows)));
        temp2=double(imag(margin_col1:margin_col2,round(0.35*numrows):round(0.65*numrows)));
        if std(temp1(:))>0.25*std(temp2(:))
            margin_col2=numcols;
        else
            margin_col2=margin_col2-1;
        end
    else
        margin_col2=numcols;
    end
    
    %normalize image intensity 
    croped_imag=imag(margin_col1:margin_col2,margin_row1:margin_row2);
    equal_imag=croped_imag;
    %equal_imag=double(adapthisteq(equal_imag,'clipLimit',0.05,'Distribution','rayleigh'));%Contrast-limited adaptive histogram equalization (CLAHE)
    equal_imag=equal_imag-imgaussfilt(equal_imag,round(numcols/10));
    equal_imag=(10000.0*(equal_imag-mean(equal_imag(:)))/std(equal_imag(:)));
    imag=zeros(size(imag));
    imag(margin_col1:margin_col2,margin_row1:margin_row2)=equal_imag; 

    imfull=imag;
    if ~cosine_sample && rotation_xaxis
    margin_up=max(margin_row1,1+round(round(numrows/2)*(1-cos(angles(idx)*pi/180))));
    margin_down=min(margin_row2,round(round(numrows/2)*(1+cos(angles(idx)*pi/180))));
    else
    margin_up=margin_row1;
    margin_down=margin_row2;
    end

    if ~cosine_sample && ~rotation_xaxis
    margin_left=max(margin_col1,1+round(round(numcols/2)*(1-cos(angles(idx)*pi/180))));
    margin_right=min(margin_col2,round(round(numcols/2)*(1+cos(angles(idx)*pi/180))));
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


    %figure(1)
    %balanced_imshow(imfull);
    %mask=zeros(size(imag));
    %mask(:,margin_up:margin_down)=1;
    proj_data_mat(:,idx,:)=permute(imfull,[1 3 2]); % order should be: column(=x), angle, rows=y
    %mask_data_mat(:,idx,:)=permute(mask,[1 3 2]); % order should be: column(=x), angle, rows=y
end
proj_id = astra_mex_data3d('create', '-proj3d', proj_geom, proj_data_mat);
% Set up the parameters for a reconstruction algorithm using the GPU

cfg = astra_struct('SIRT3D_CUDA');
%cfg = astra_struct('CGLS3D_CUDA');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = proj_id;
%cfg.option.SinogramMaskId=mask_id; %only with sirt3d and helps only at the margins of the volume slices
% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);
% Run 50 iterations of the algorithm
astra_mex_algorithm('iterate', alg_id,150);
% Get the result
rec = astra_mex_data3d('get', rec_id);%maybe 'get_single'
%errorP=astra_mex_algorithm('get_res_norm', alg_id)/((nX/bin)*(nY/bin)*length(angles));
%Save to new MRC names rec_...
newFilename=strrep(Chosen_Filename_ali1,'.ali.','_dual.rec_SIRT.');
newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=newFilename;
newmRCImage = setVolume(newmRCImage, rec); %enter to newmRCImage, do statistics, and fill many details to the header
newmRCImage.header.cellDimensionX = sizeXangstrom;
newmRCImage.header.cellDimensionY = sizeYangstrom;
newmRCImage.header.cellDimensionZ = sizeZangstrom;
%MODE=2:  32-bit signed real
save(newmRCImage, newFilename);
close(newmRCImage);
%disp(sprintf('Saved to file: %s ',newFilename));
%fprintf('Difference between source projections and reconstructions: %d \n',errorP);

astra_mex_data3d('clear');%free GPU memory

end %of main function


%####################################################
function r_mn=r_mn_rect(Imagem,Imagen,shift_limit,do_filt)
    if do_filt==1
        Imagem=imgaussfilt(Imagem-imgaussfilt(Imagem,30),3);
        Imagen=imgaussfilt(Imagen-imgaussfilt(Imagen,30),3);
    end

    %figure(2);
    %subplot(1,2,1);
    %balanced_imshow(Imagem);
    %subplot(1,2,2);
    %balanced_imshow(Imagen);
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
            r_mn=[0 0];
        else
            r_mn=[yoffset xoffset];
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
    r_mn=[yoffset2 xoffset2];
end


