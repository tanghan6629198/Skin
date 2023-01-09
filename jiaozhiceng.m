

close all;clear all;clc;
warning off;
%% 读OCT数据
% path='C:\Users\汤瀚\Desktop\';
% imgfile=strcat(path,'goodskin.png'); %skin_pro  ,goodskin  ,sk22_raw
% 
% img=imread(imgfile);
% img = imfilter(img,fspecial('gaussian',[3 3],3));
% figure();imagesc(img);colormap('gray');title('强度数据');  
% imagesc(img); axis image; colormap('gray'); hold on;
% axis on
 handle = OCTFileOpen('F:\UV-VC-BC\th20220428\3-1\th20220428_0037_Mode2D.oct');
img = OCTFileGetIntensity(handle);
img3 = img;
img = imfilter(img,fspecial('gaussian',[3 3],3));
figure();imagesc(img);colormap('gray');title('强度数据');  
imagesc(img); axis image; colormap('gray'); hold on;
axis on;
%% 裁剪图片
img2 = img;
r = 285;
centerY = 500;
% up = inputdlg('Top of image');
% low = inputdlg('Bottom of image');
% up = transpose(str2num(cell2mat(up)));
% low = transpose(str2num(cell2mat(low)));
up =450;
low=800;

img = img(up:low,centerY-r:centerY+r);
img2 = img2(up:low,centerY-r:centerY+r);
img3 = img3(up:low,centerY-r:centerY+r);
figure();imagesc(img);colormap('gray');title('强度数据');  
imagesc(img); axis image; colormap('gray'); hold on;
axis on;
% %%%%%%%%%%%%%%%%%z方向像素大小%%%%%%%%%%%%%%%%
% RefractiveIndex = OCTFileGetProperty(handle, 'RefractiveIndex');
deta = 0.003493;
% Redeta = 0.003493/str2num(RefractiveIndex);
Redeta = 0.003493/1.38;
%%%%%%%%%%%%%%%%%左右加两列强度0%%%%%%%%%%%%%%  
szImg = size(img);
imgNew = zeros([szImg(1) szImg(2)+2]);
imgNew(:,2:1+szImg(2)) = img;
szImgNew = size(imgNew);
%%%%%%%%%%%%%%%%强度归一化图像%%%%%%%%%%%%%%%%%
IntensityImg = nan(szImgNew);
IntensityImg =(imgNew-min(imgNew(:)))/(max(imgNew(:))-min(imgNew(:)));%%强度归一化
% figure();imagesc(IntensityImg);colormap('gray');title('强度');
neighborIter = [1 1 1 0 0 -1 -1 -1;...%%方向
                1 0 -1 1 -1 1 0 -1];
minWeight = 1E-5;
%% 上下表面
% [pathX,pathY] = OCTGetUpskin( imgNew,szImgNew,IntensityImg );
IntensityImg2 = IntensityImg;
% [pathX2,pathY2] = OCTGetLowskin( pathX,imgNew,szImgNew,IntensityImg2 );
%% 角质层
% get  vertical gradient image 获得垂直梯度图像
gradImg = nan(szImgNew);
for i = 1:size(imgNew,2)
    gradImg(:,i) = -1*gradient(imgNew(:,i),2); %%求每个像素点垂直梯度，2是均匀间距
end
gradImg = (gradImg-min(gradImg(:)))/(max(gradImg(:))-min(gradImg(:)));%%亮→暗 
figure();imagesc(gradImg);colormap('gray');title('亮→暗');
gradImg2 = gradImg*-1+1;%暗→亮
[pathX,pathY] = OCTGetUpskin( imgNew,szImgNew,gradImg2 );
[pathX2,pathY2] = OCTGetLowskin( pathX,imgNew,szImgNew,IntensityImg2 );
[pathX3,pathY3] = OCTGetCuticle(pathX,pathX2,imgNew,szImgNew,gradImg);

pathX = pathX(1,2:2*r+2);
pathY = 1:2*r+1;
pathX2 = pathX2(1,2:2*r+2);
pathY2 = 1:2*r+1;
pathX3 = pathX3(1,2:2*r+2);
pathY3 = 1:2*r+1;
%% 分层
imagesc(img2); axis image; colormap('gray'); hold on;
plot(pathY,pathX,'r-','linewidth',1); hold on;
plot(pathY2,pathX2,'g-','linewidth',1); hold on;
plot(pathY3,pathX3,'b-','linewidth',1); hold on;

%% 求皮肤厚度和角质层厚度
if pathX(1,1)<pathX2(1,1)
    Upskin = mean(pathX,2)
    Lowskin = mean(pathX2,2)
else
    Upskin = mean(pathX2,2)
    Lowskin = mean(pathX,2)
end

Cuticle = mean(pathX3,2)

if pathX(1,1)<pathX2(1,1)
    CuticleThickness = abs((pathX-pathX3)*Redeta*1000);
else
    CuticleThickness = abs((pathX2-pathX3)*Redeta*1000);
end


thickness = roundn(abs((Lowskin-Upskin)*Redeta*1000),-4)  
Cuticlethickness = roundn(abs((Cuticle-min(Upskin,Lowskin))*Redeta*1000),-4)      %μm
q = 0;
for i =1:571
    q = q + (CuticleThickness(1,i)- Cuticlethickness)^2;
end
Rq = sqrt(q/571)  % 均方根高度


% Cuticlethickness = roundn(abs((Cuticle-min(Upskin,Lowskin))*Redeta*1000) + 3.493/1.38 ,-4) 
% roundn(abs((mean(pathX3(1:50),2)-mean(pathX2(1:50),2))*Redeta*1000),-4)

%% ROI区域纹理特征分析
% imgI = 0;
% %拉平
%     for i = 1:100
% %         depth =  max(pathX2(1,begin+i),pathX(1,begin+i)) - pathX3(1,235+i) ;
% %         imgI(50-depth:50,i) = img3(pathX3(1,235+i): max(pathX2(1,begin+i),pathX(1,begin+i)),235+i);
%         begin = 35;
%         depth = max(pathX2(1,begin+i),pathX(1,begin+i)) - pathX3(1,begin+i) ;
%         imgI(50-depth:50,i) = img3(pathX3(1,begin+i):max(pathX2(1,begin+i),pathX(1,begin+i)),begin+i);
%         imgI = imgI(1:48,:);
% 
% 
%     end
%     figure();imagesc(imgI); axis image; colormap('gray'); hold on;
%     
%     I = imgI(29:48,:);
%     Offset = [0,3];  %0° d=3
%     N = 8; %灰度等级
% [ASM,ENT,COR] = graycomatrix(Offset,I,N);



% Upskin = mean(pathX,2)
% Lowskin = mean(pathX2,2)
% Cuticle = mean(pathX3,2);
% thickness = abs((Lowskin-Upskin)*Redeta)
% Upskin = mean(pathX(1,1:200));
% Cuticle = mean(pathX3(1,1:200))
% Cuticlethickness = abs((Cuticle-Upskin)*Redeta)


% A150 = img2(:,81);
% A150 = A150';
% A300 = img2(:,231);
% A300 = A300';
% A450 = img2(:,381);
% A450 = A450';
% A600 = img2(:,531);
% A600 = A600';
% plot(A150);
% plot(A300);
% plot(A450);
% plot(A600);