% Comments added by Ashwin Renganathan 4/1/2015
clear all
close all
clc

g_name  = 'lincoln.bmp';        % Given grayscale image name
c_name  = 'lincoln_marked.bmp'; % color scribbled grayscale name
out_name= 'lincol_res.bmp';    % Final colorized name

%set solver=1 to use a multi-grid solver 
%and solver=2 to use an exact matlab "\" solver
solver=2; 

gI=double(imread(g_name))/255;    % normalized input image, varies from 0:1
cI=double(imread(c_name))/255;    % normalized scribbled, varies from 0:1

%ashwin's mod
% gI(:,:,1) = gI;
% gI(:,:,2) = gI(:,:,1);
% gI(:,:,3) = gI(:,:,1);
% cI(:,:,1) = cI;
% cI(:,:,2) = cI(:,:,1);
% cI(:,:,3) = cI(:,:,1);

% Question - Why do we need this?
abs(gI - cI);
colorIm=(sum(abs(gI-cI),3)>0.01); % take abs difference of above images, add them up in the channel direction and keep only > 0.01 values
                                  % Note: the result is a matrix with '1'
                                  % everywhere the value is non zero
                                  % why add them? -  may be because it
                                  % lowers the dimensionality?
colorIm=double(colorIm);          % convert above matrix to double precision

sgI=rgb2ntsc(gI);                 % convert an rgb format image to ntsc which has YIQ where Y is luminance, I & Q are Chorminance. 
                                  % Note: the channels 2 and 3 are all '0's
scI=rgb2ntsc(cI);                 % Note: the channels 2 and 3 are not '0's 
                                  % and that is why we are interested in
                                  % them
   
ntscIm(:,:,1)=sgI(:,:,1); % Copy all the 1st channel from input image (Luminance)
ntscIm(:,:,2)=scI(:,:,2); % Copy all the 2nd channel from scribbled (Chrominance)
ntscIm(:,:,3)=scI(:,:,3); % Copy all the 3rd channel from scribbled (Chrominance)

% This part copies a select portion of each image and for that purpose 
% identifies the indices that form the bound within which it has to be
% copied

% Question -  What is the logic behind this?_ *
max_d = floor(log(min(size(ntscIm,1),size(ntscIm,2)))/log(2)-2);
iu    = floor(size(ntscIm,1)/(2^(max_d-1)))*(2^(max_d-1));
ju    = floor(size(ntscIm,2)/(2^(max_d-1)))*(2^(max_d-1));

id=1; jd=1;
colorIm = colorIm(id:iu,jd:ju,:);
ntscIm  = ntscIm(id:iu,jd:ju,:);

if (solver == 1) % Multigrid SOlver
  nI=getVolColor(colorIm,ntscIm,[],[],[],[],5,1);
  nI=ntsc2rgb(nI);
else             % LU Factorization
  nI=getColorExact(colorIm,ntscIm);
  % getColorExact(colorIm, ntscIm)
  % Arguments - colorIm=(sum(abs(gI-cI),3)>0.01);
  % Arguments - ntscIm which is the output of rgb2ntsc,
  %             but contains only outputs in the id:iu range    
end

figure, imshow(nI)

imwrite(nI,out_name)
   
  

%Reminder: mex cmd
%mex -O getVolColor.cpp fmg.cpp mg.cpp  tensor2d.cpp  tensor3d.cpp
