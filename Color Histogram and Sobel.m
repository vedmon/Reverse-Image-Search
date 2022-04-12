clc
clear all;
count = 0;
delete('D:\Photos\ABC\*.jpg');
srcFiles = dir('D:\Photos\ImageDatabase\*.jpg');  % the folder in which ur images exists
   
    
    %%%%%COLOR DETECTION BY HISTOGRAMS%%%%%

for i = 1 : length(srcFiles)
    filename = strcat('D:\Photos\ImageDatabase\',srcFiles(i).name);
    %I = imread(filename);
    %figure, imshow(I);
    
    Im1 = imread('D:\Photos\IMG_4134.jpg');
    Im2 = imread(filename);
    %Im2 = imread('F:\Photos\SPBM2016\IMG_4134.jpg');
    %imshow(Im1);
    %imshow(Im2);

    Red1 = Im1(:, :, 1);    Red2 = Im2(:, :, 1);
    Green1 = Im1(:, :, 2);  Green2 = Im2(:, :, 2);
    Blue1 = Im1(:, :, 3);   Blue2 = Im2(:, :, 3);

    HnRed1 = imhist(Red1)./numel(Red1);         HnRed2 = imhist(Red2)./numel(Red2); 
    HnGreen1 = imhist(Green1)./numel(Green1);   HnGreen2 = imhist(Green2)./numel(Green2); 
    HnBlue1 = imhist(Blue1)./numel(Blue1);      HnBlue2 = imhist(Blue2)./numel(Blue2);

    FRed = sum((HnRed1 - HnRed2).^2);
    FGreen = sum((HnGreen1 - HnGreen2).^2);
    FBlue = sum((HnBlue1 - HnBlue2).^2);

    %F = Alpha*FBlue + Beta*FRed + Gamma*FGreen
    F = 0.2989*FRed + 0.5870*FGreen + 0.1140*FBlue;
    
    if  F <= 0.003
        count = count+1;
        disp(count);
        disp(F);
        str = sprintf('D:\\Photos\\ABC\\%d.jpg',count); % e.g. "1.png"
        disp(str);
        %str = 'F:\Photos\ABC\dog.jpg'
        imwrite(Im2,str,'jpg');
        %imwrite(Im1,'F:\Photos\ABC\'baseFileName.jpg','jpg');
    end
    
    %if count <= 1
    %    disp(F);
    %end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    
    %%%%%EDGE DETECTION BY SOBEL METHOD%%%%%

%close all;
%clc;
img=imread('D:\Photos\SPBM2016\IMG_4134.jpg');
%img=Im1;
B=rgb2gray(img);
%subplot(2,2,1)
%imshow(B)
pause(2)
I=double(B);

for i=1:size(I,1)-2
for j=1:size(I,2)-2
%Sobel mask for x-direction:
mx=((2*I(i+2,j+1)+I(i+2,j)+I(i+2,j+2))-(2*I(i,j+1)+I(i,j)+I(i,j+2)));
%Sobel mask for y-direction:
my=((2*I(i+1,j+2)+I(i,j+2)+I(i+2,j+2))-(2*I(i+1,j)+I(i,j)+I(i+2,j)));

B(i,j)=sqrt(mx.^2+my.^2);
end
end
%subplot(2,2,2)
%imshow(B); title('Sobel gradient');
pause(2)
%Define a threshold value


Thresh=100;
B=max(B,Thresh);
B(B==round(Thresh))=0;
B=uint8(B);

%subplot(2,2,1)
%imshow(~B);title('Edge detected Image')
%subplot(2,2,2)
%imhist(~B);
%c=imhist(~B);
i1=~B;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
match=1;
srcFiles2 = dir('D:\Photos\ABC\*.jpg');
arr = 1 : length(srcFiles2);
count1=1;
for i = 1 : length(srcFiles2)


%close all;
%clc;
%img=imread('F:\Photos\SPBM2016\IMG_4136.jpg');
filename = strcat('D:\Photos\ABC\',srcFiles2(i).name);
img = imread(filename);
B=rgb2gray(img);
%subplot(2,2,3)
%imshow(B)
%pause(2)
I=double(B);

for i=1:size(I,1)-2
for j=1:size(I,2)-2
%Sobel mask for x-direction:
mx=((2*I(i+2,j+1)+I(i+2,j)+I(i+2,j+2))-(2*I(i,j+1)+I(i,j)+I(i,j+2)));
%Sobel mask for y-direction:
my=((2*I(i+1,j+2)+I(i,j+2)+I(i+2,j+2))-(2*I(i+1,j)+I(i,j)+I(i+2,j)));

B(i,j)=sqrt(mx.^2+my.^2);
end
end
%subplot(2,2,2)
%imshow(B); title('Sobel gradient');
%pause(2)
%Define a threshold value


Thresh=100;
B=max(B,Thresh);
B(B==round(Thresh))=0;
B=uint8(B);

%subplot(2,2,3)
%imshow(~B);title('Edge detected Image')
%subplot(2,2,4)
%imhist(~B);
%d=imhist(~B);
%difference=pdist('c','d')
%i1=imread();
i1=i1(:,:,1);
[c1,n]=imhist(i1);
c1=c1/size(i1,1)/size(i1,2);
i2=~B;
i2=i2(:,:,1);
[c2,n2]=imhist(i2);
c2=c2/size(i2,1)/size(i2,2);
d=pdist2(c1',c2')

if d < match
    match = d;
end
arr(count1) = d
count1 = count1 + 1;
end

result = match

for i = 1 : length(srcFiles2)
    filename = strcat('D:\Photos\ABC\',srcFiles2(i).name);
    if arr(i) == result
        img = imread(filename);
        imshow(img); title('Result')
    end
end

