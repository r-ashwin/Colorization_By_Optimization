function [nI,snI]=getColorExact(colorIm,ntscIm)
% Comments added by Ashwin Renganathan 4/1/2015

% colorIm - contains absolute diff between input and scribbled image
%           and summing up all 3 channels.
% ntscIm  - has YIQ form. Has intensity from input-image and chrominance
%           from the scribble image
% nI      - the final colorized image
% snI     - 


n=size(ntscIm,1); % Rows in the stripped YIQ file
m=size(ntscIm,2); % Cols in the stripped YIQ file
imgSize=n*m;      % image size

nI(:,:,1) = ntscIm(:,:,1);          % just copying "intensity channel" of ntscIm onto the output nI
indsM = reshape([1:imgSize],n,m);   % reshape a vector of length imgSize to a matrix if size nxm
%this is how reshape works
% reshape([1:6],2,3)
%
% ans =
%
%     1     3     5
%     2     4     6

lblInds = find(colorIm);              % find indices of non-zero elements

wd=1; 

len=0;
consts_len=0;
col_inds=zeros(imgSize*(2*wd+1)^2,1); % initialize a matrix of zeros. Size (n*m*9,1) 
row_inds=zeros(imgSize*(2*wd+1)^2,1); % initialize a matrix of zeros. Size (n*m*9,1) 
vals=zeros(imgSize*(2*wd+1)^2,1);     % initialize a matrix of zeros. Size (n*m*9,1) 
gvals=zeros(1,(2*wd+1)^2);            % initialize a matrix of zeros. Size (1,9) 


for j=1:m 
   for i=1:n
      consts_len=consts_len+1;
      
      if (~colorIm(i,j))                        % iterating through each element of colorIm, returns 1 if == 0, 0 otherwise
                                                % we are basically
                                                % attempting to access the
                                                % parts of the image that
                                                % has no scribblings
        tlen=0;
        for ii=max(1,i-wd):min(i+wd,n)          % iterate from i-1 to i+1. Just 3 elements. but at boundaries 2 elements  
           for jj=max(1,j-wd):min(j+wd,m)       % iterate from j-1 to j+1. Just 3 elements
            
              if (ii~=i)||(jj~=j)               % while ii neq 1 or n, jj new 1 or n
                 len=len+1; tlen=tlen+1;        % increment len and tlen
                 row_inds(len)= consts_len;     % 
                 col_inds(len)=indsM(ii,jj);    % indsM just contains integers from 1:imgSize 
                 gvals(tlen)=ntscIm(ii,jj,1);   % Copy ii,jj element of 1st channel of ntscIm into gvals
              end
           end
        end
        t_val=ntscIm(i,j,1); % Copy i,j element of 1st channel of ntscIm into t_val
        gvals(tlen+1)=t_val; % copy t_val to gvals(tlen+1)
        c_var=mean((gvals(1:tlen+1)-mean(gvals(1:tlen+1))).^2); % something like mean([y - <y>]^2)
        csig=c_var*0.6; 
        mgv=min((gvals(1:tlen)-t_val).^2);
        
        if (csig<(-mgv/log(0.01)))
	   csig=-mgv/log(0.01);
        end
        if (csig<0.000002)
	   csig=0.000002;
        end

        gvals(1:tlen)=exp(-(gvals(1:tlen)-t_val).^2/csig);
        gvals(1:tlen)=gvals(1:tlen)/sum(gvals(1:tlen));
        vals(len-tlen+1:len)=-gvals(1:tlen);
      end

        
      len=len+1;
      row_inds(len)= consts_len;
      col_inds(len)=indsM(i,j);
      vals(len)=1; 

   end
end

       
vals=vals(1:len);
col_inds=col_inds(1:len);
row_inds=row_inds(1:len);


A=sparse(row_inds,col_inds,vals,consts_len,imgSize);
b=zeros(size(A,1),1);


for t=2:3
    curIm=ntscIm(:,:,t);
    b(lblInds)=curIm(lblInds);
    new_vals=A\b;   
    nI(:,:,t)=reshape(new_vals,n,m,1);    
end



snI=nI;
nI=ntsc2rgb(nI);

