function [overlap, labelling] = segment_gc(unary, current_labels, pairc, MDL,im) % ~,
%Graphcut based segmentation that minimizes the energy
%\sum_{p\in P} { Unary(p,label(p)) + \Sum_{q\in N_p} pairc Delta(label(p)\neq label(q)) } +\MDL*length(unique(label))
% where p is a pixel of the image
%
% (INPUTS):
%
% unary: Nrows x Ncols x Nactlabls0 array where unary(i,j,i_lab) is the
% cost of assigning the label with index i_lab to the pixel (i,j). Here,
% this cost comes from average (over frames) photometric error evaluated
% using the motion parameters that correspond to the label.
%
% lablfunc0: intialization of the labeling function as a Nrows x Ncols
% array: lablfunc0(i,j) is the label of pixel (i,j) at initialization.
%
%
% pair,MDL: weighting constants of the energy to be minimized.
%

unary = double(unary);
pairc = double(pairc);
current_labels = double(current_labels);
MDL = double(MDL);
% 


if nargin<5
    im=[];
end

% % % % % test:
% pairc, MDL,im
% pause;
% % % % % end of test

width=size(im,1);
height=size(im,2);
index=1:(width*height);
index=reshape(index,width,height);
index=padarray(index,[1 1]);
pair=zeros(8,width*height);
pair(1,:)=reshape(index(1:width,2:(height+1)),1,width*height);
pair(2,:)=reshape(index(2:(width+1),1:(height+0)),1,width*height);
pair(3,:)=reshape(index(3:(width+2),2:(height+1)),1,width*height);
pair(4,:)=reshape(index(2:(width+1),3:(height+2)),1,width*height);% Principle axis
pair(5,:)=reshape(index(1:(width),1:(height)),1,width*height);
pair(6,:)=reshape(index(3:(width+2),1:(height)),1,width*height);
pair(7,:)=reshape(index(1:(width),3:(height+2)),1,width*height);
pair(8,:)=reshape(index(3:(width+2),3:(height+2)),1,width*height); %Diagonals
pair = sort(pair);
pair=pair(end:-1:1,:);
if (isempty(im))
    pcost=pairc*ones(width*height,4);
else
    im = double(im);
    im(:,:,1) = medfilt2(im(:,:,1),[3 3]);
    im(:,:,2) = medfilt2(im(:,:,2),[3 3]);
    im(:,:,3) = medfilt2(im(:,:,3),[3 3]);
    figure();
    t=double(reshape(im,size(im,1)*size(im,2),3));
    t=t(1:size(pair,2),:);
    t2=repmat(t, [1 1 8]);
    t3=zeros(size(t2));
    for i=1:8
        t3(:,:,i)=t(max(1,pair(i,:)),:);
    end
    pcost=-sum(abs(t3-t2),2);
    pcost=reshape(pcost,size(pcost,1),8);
    
    pcost = exp(pcost/20);    
    
end
%unary=reshape(unary,width*height,size(unary,1));

% % % % test
% unary;any(isnan(unary(:))),pause;
% pair;any(isnan(pair(:))),pause;
% pcost';any(isnan(pcost(:))),pause;
% % current_labels(:)'-1;pause;
% % 0.01*MDL*ones(size(unary,2)),pause;
% % % % 
imshow(reshape(mean(pcost')/2,width,height)*5)
% [~,labelling]=expand(unary,pair,pcost'*1000*pairc,current_labels(:)'-1,0*ones(size(unary,2),1),Inf);%1000
% figure();
% t=randperm(max(labelling)+1);imagesc(reshape(t(labelling+1),width,height)');
% figure();
% out=6.2037e+03;max(unary(sub2ind(1:length(labelling),labelling+1)));
out=300;
%[overlap,labelling]=allgc(unary,pair,pair,pcost'*100,current_labels(:)-1,1,1000*ones(size(unary,2),1),out);
[overlap,labelling,out]=allgc(unary,pair,pair,(1+pcost)'*10,current_labels(:)-1,1,150*ones(size(unary,2),1),220);
 plotoverlap(overlap,out,im);
 figure();
 t=randperm(max(labelling)+1);imagesc(reshape(t(labelling+1),width,height)');
labelling=labelling'+1;

labelling = reshape(labelling,width,height);
%
% (OUTPUTS):
%
% lablfunc, actlabls: result of the segmentation, represented by the
% labeling function and active-labels vector correspondingly. The format is
% the same as in the input variables lablfunc0, actlabls0. Note that the
% solution might have less active labels than the initialization, which
% means that the length of actlabls might be smaller than that of
% actlabls0.



