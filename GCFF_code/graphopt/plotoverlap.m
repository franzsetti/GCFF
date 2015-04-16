function a= plotoverlap(ov,out,im)
ov=reshape(sum(ov,2),size(im,1),size(im,2));
out=reshape(sum(out,2)~=0,size(im,1),size(im,2));
edge=2*(ov>1)+(ov==1);
edge=repmat(edge/2, [1 1 3]);
out=repmat(out, [1 1 3]);

 im(:,:,1)=im(:,:,1).*(ov>=1);
 im(:,:,2)=im(:,:,2).*(ov>=1);
 im(:,:,3)=im(:,:,3).*(ov>=1);
 edge2=im/255;
 im(:,:,1)=im(:,:,1).*(ov>1);
 im(:,:,2)=im(:,:,2).*(ov>1);
 im(:,:,3)=im(:,:,3).*(ov>1);
 imshow([edge,edge2,out,im/255]);
a=true;