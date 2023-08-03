function R=gray2rgb(img1,img2)
tic
 clc;
 warning off;
 imt=imread(img1);
 ims=imread(img2);
 [sx sy sz]=size(imt);
 [tx ty tz]=size(ims);
 if sz~=1
     imt=rgb2gray(imt);
 end
 if tz~=3
     disp ('img2 must be a color image (not indexed)');
 else
     imt(:,:,2)=imt(:,:,1);
     imt(:,:,3)=imt(:,:,1);

     nspace1=rgb2ycbcr(ims);
     nspace2= rgb2ycbcr(imt);

    ms=double(nspace1(:,:,1));
     mt=double(nspace2(:,:,1));
     m1=max(max(ms));
     m2=min(min(ms));
     m3=max(max(mt));
     m4=min(min(mt));
     d1=m1-m2;
     d2=m3-m4;
     dx1=ms;
     dx2=mt;
     dx1=(dx1*255)/(255-d1);
     dx2=(dx2*255)/(255-d2);
     [mx,my,mz]=size(dx2);
     for i=1:mx
         for j=1:my
              iy=dx2(i,j);
              tmp=abs(dx1-iy);
              ck=min(min(tmp));
              [r,c] = find(tmp==ck);
              ck=isempty(r);
              if (ck~=1)            
                  nimage(i,j,2)=nspace1(r(1),c(1),2);
                  nimage(i,j,3)=nspace1(r(1),c(1),3);
                  nimage(i,j,1)=nspace2(i,j,1);           
             end
          end
      end
     rslt=ycbcr2rgb(nimage)
     R=uint8(rslt);
     toc
 end