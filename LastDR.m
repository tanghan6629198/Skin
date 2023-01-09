clc
clear 
warning off
%*******************************************************

%功能：批量计算3D.data数据散射系数，先分角质层再统计3D散射系数
%完成度：完成
%码农：Tanghan
%时间：2022.5.9
%Matlab版本：2020a

%*******************************************************

 %% 进行3D数据处理

dataFilePath =  uigetdir('Select unziped 2D data path');     %选择数据文件夹
tic
if isequal(dataFilePath,0)    % 判断是否正常选择文件夹
   disp('User selected Cancel') ;
else
    disp(['User selected ', dataFilePath]);
end

resultDir = strcat(dataFilePath,'\','result_SC');
if exist(resultDir,'dir') == 0
    mkdir(resultDir);
end

dataSubPath = dir(dataFilePath);
wb1 = waitbar(0,'处理进度...');
toc

for iwb = 1:length(dataSubPath)
    
    
    if(isequal(dataSubPath(iwb).name,'.')||... % 去除系统自带的两个隐文件夹
        isequal(dataSubPath(iwb).name,'..')||...
        ~dataSubPath(iwb).isdir) % 去除遍历中不是文件夹的
        continue;
    end
    dataDir = dir([dataFilePath '\' dataSubPath(iwb).name '\*.data']);
    for j =1:length(dataDir) % 遍历所有文件
        fileName = [dataFilePath '\' dataSubPath(iwb).name '\' dataDir(j).name]
        disp(fileName);

%% 改进后的DR算法，先拟合然后DR
%   miu_Origin 三维分割结果
%286-286 最中心Aline
%取中间401*401区域计算角质层和散射系数，最终统计半径200区域
%% 基本参数设置
A_scan=1024;  % 采样点数
B_scan=1000;  %A-Scan的数量
line=1000;
B_Length = 10; %mm
Redeta = 0.003493/1.38;
%% 读取数据
% dataFilePath =  uigetdir('Select unziped 3D data path');     %选择要处理的文件
% dataDir = dir([dataFilePath '\*.data']);
%  for j =1:length(dataDir) % 遍历所有文件
% 	fileName = [dataFilePath '\' dataDir(j).name]
%     disp(fileName);  
%  end

shft = B_scan*A_scan;  %一个bscan总的采样点数
fid=fopen(fileName,'r','n');
wb = waitbar(0,'读取3D数据中...');
disp('读取3D数据中...');
data3D = zeros(A_scan,B_scan,line,'single');
data3D_1 = zeros(A_scan,B_scan,line,'single');
A = zeros(A_scan,B_scan,'single');
for i=1:line   
    A = fread(fid,[A_scan,B_scan],'float32','n'); 
    if isempty(A)==1    %表示逻辑运算符“非”，也就是取反
        break;  %当A是空矩阵时，跳出
    end
    fseek(fid, 4*i*shft, 'bof');  %指针移入第j*shft个元素
    data3D(:,:,i)=A;
       for dataj = 1:A_scan
          data3D_1(dataj,:,i) = smooth(data3D(dataj,:,i),11);   %滑动平均
       end
    waitbar(i/line,wb,['读取3D数据中...' num2str(100*i/line) '%']);
end
for h=1:B_scan
    for k=1:A_scan
        data3D_1(k,h,:) = smooth(data3D_1(k,h,:),11);
    end
end
clear A;
close(wb);
disp('读取3D数据完成');
%% 裁剪有效区域
% 确定 AScan 的有效范围
section2D = zeros(A_scan,B_scan,'single')*NaN;
section2D(:) = data3D_1(:,:,line/2);
figure(1);
clf;
imagesc(section2D);

rangeMin = 400;
rangeMax = 700;
range = rangeMin:rangeMax;  %皮肤信号范围

close(figure(1));
centerX = B_scan/2;
centerY = line/2;
r=200;
r2 = 200*200;
data3D_1 = data3D_1(rangeMin:rangeMax,centerX-r:centerX+r,centerY-r:centerY+r);
p1 = [centerX-r,centerX+r];   %x轴的坐标
p2 = [rangeMax,rangeMin];   %y轴的坐标

IndexK = 0;
%% 分割角质层 上下边界
for i = 1:2*r+1
    img = data3D_1(:,:,i);
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
[pathX1,pathY] = OCTGetUpskin( imgNew,szImgNew,IntensityImg );

IntensityImg2 = IntensityImg;
[pathX2,pathY2] = OCTGetLowskin( pathX1,imgNew,szImgNew,IntensityImg2 );

%% 角质层
% get  vertical gradient image 获得垂直梯度图像
    gradImg = nan(szImgNew);
    for L = 1:size(imgNew,2)
        gradImg(:,L) = -1*gradient(imgNew(:,L),2); %%求每个像素点垂直梯度，2是均匀间距
    end
    gradImg = (gradImg-min(gradImg(:)))/(max(gradImg(:))-min(gradImg(:)));%%亮→暗 
    figure();imagesc(gradImg);colormap('gray');title('亮→暗');
    close
    gradImg2 = gradImg*-1+1;%暗→亮
%     [pathX,pathY] = OCTGetUpskin( imgNew,szImgNew,gradImg2 );
%     [pathX2,pathY2] = OCTGetLowskin( pathX,imgNew,szImgNew,IntensityImg2 );
    [pathX3,pathY3] = OCTGetCuticle(pathX1,pathX2,imgNew,szImgNew,gradImg);
    pathX1 = pathX1(1,2:2*r+2);
    pathX2 = pathX2(1,2:2*r+2);
    pathX3 = pathX3(1,2:2*r+2);
    pathX(:,:,i) =[pathX1;pathX3;pathX2];
end
%% 计算散射值
wbSkin3dProcess = waitbar(100,'计算皮肤初始散射值...');
disp('计算皮肤初始散射值...');
detaz = 0.003493/1.38;
FOV = 2.59; %MM
for k =1:2*r+1  % 每幅B-Scan
    [n,m] = size(data3D_1(:,:,k));  %获取B-Scan的 大小   n是AScan的范围rangeFix-rangeMin+1；m=571
    for indexi = 1:n                     %强度信号转化成功率信号
        for indexj = 1:m
            data3D_1(indexi,indexj,k) = 10^(data3D_1(indexi,indexj,k)/20);   %dB转化成 强度
        end
    end    
end
%% 强度归一化

for k =1:2*r+1   % 每幅B-Scan
    fft_intensity = data3D_1(:,:,k);
      % % %     加入灵敏因子计算
    z = detaz*(1:1024);
    s_z = exp(-z.^2/25);  %S(z)--灵敏度因子---考虑了灵敏度衰减所引起的影响，该方法还采用了成像深度Z和宽度dert的高斯模型，dert=5；
    s_z_fix = s_z(1,rangeMin:rangeMax)';  %竖着一列 rangeMaxFix-rangeMin+1 个像素点
    miuOrigin = fft_intensity./s_z_fix;
    for choseAScan = 1:2*r+1 %每条A-scan
            sur = pathX(1,choseAScan,k);
            cul = pathX(2,choseAScan,k);
            bot = pathX(3,choseAScan,k);
           
            while miuOrigin(bot,choseAScan) >= miuOrigin(cul,choseAScan)
                bot = bot-1;
            end
            miuBegin = cul;
            miuEnd = bot;
            if miuEnd > miuBegin
            Intensity2 = miuOrigin(miuBegin:miuEnd,choseAScan);  
            I1 = Intensity2;
            Intensity2 = log(Intensity2);
            [fitX,~] = size(Intensity2);
            fitX = ((1:fitX)*detaz)'; % 纵向的像素大小
            fixA = polyfit(fitX,Intensity2,1); %获得拟合的两个数值
            a = fixA(1)/-2 ; % 对斜率/2 获得拟合散射值
            CF(choseAScan,k) = a;
            szI1 = size(I1,1);
            miu1 = zeros(1,szI1);


            miu1(szI1) = a;
            for k1 = szI1:-1:2
                sumI = sum(I1(k1:szI1));
                summiu1 = sum(miu1(k1:szI1));
                miu1(k1-1) = I1(k1)/(sumI*2*detaz)*(1-exp(-2*detaz*summiu1));
            end
            
            
            miu3D(1:szI1,choseAScan,k) = miu1;
            miu3D_mean(choseAScan,k) = mean(miu1);
            else
                CF(choseAScan,k) = 0;
                miu3D(1,choseAScan,k) = 0;
                miu3D_mean(choseAScan,k) = 0;
            end

    end
    waitbar(k/line,wbSkin3dProcess,['计算皮肤初始散射值...' num2str(100*k/line) '%']);
end
 close(wbSkin3dProcess);
 disp('计算皮肤散射值完成');   
 
%% 取半径200px的统计区域
miu_CF=0;miu_DR=0;count = 0;sumALL = 0;
for indexB = 1:2*r+1
	for indexA = 1:2*r+1
        if((indexA-201)^2+(indexB-201)^2)>r2
          	CF(indexA,indexB)=nan;
            miu3D_mean(indexA,indexB)=nan;
        else
                if CF(indexA,indexB) < 0
                    CF(indexA,indexB) = 0;
                end
                    if miu3D_mean(indexA,indexB) < 0
                        miu3D_mean(indexA,indexB) = 0;
                    end     
        end
	end
end
%% 计算CF和DR的平均散射系数

for indexB = 1:2*r+1
	for indexA = 1:2*r+1
        if ~isnan(CF(indexA,indexB))
            count = count+1;
            miu_CF = miu_CF + CF(indexA,indexB);
            miu_DR = miu_DR + miu3D_mean(indexA,indexB);
            sumALL = sumALL + sum(miu3D(:,choseAScan,k));
        end
	end
end
miu_CF = miu_CF/count
miu_DR = miu_DR/count  
sumALL
numBeginBigger = 0; st = 0;
for indexB = 1:2*r+1
	for indexA = 1:2*r+1
        if ~isnan(miu3D_mean(indexA,indexB))
            st = st+(miu3D_mean(indexA,indexB)-miu_DR).^2;
            if miu3D_mean(indexA,indexB)>=miu_DR
                numBeginBigger = numBeginBigger + 1;
            end
        end
	end
end
stDR = sqrt(1/(count)*st)
numBeginBigger = numBeginBigger/count

%% DR作图
miu3D_mean2 = miu3D_mean;
fig15 = figure(15);
x=(-r:r)*0.01; %x轴的像素大小以及分布
y=(-r:r)*0.01; %y轴的像素大小以及分布
h=imagesc([x x],[y y],miu3D_mean);
set(h,'alphadata',~isnan(miu3D_mean));  % NaN 的颜色显示白色而不是值最小的颜色
set(gca,'FontSize',32);%,'FontSize',16
set(gca,'LineWidth',1.5);
colorbar;
colormap(parula);
caxis([0,5]);
set(gca,'xtick',-2:1:2);   
set(gca,'ytick',-2:1:2);  
axis([-2,2,-2,2]);
% xlabel('mm^-^1','FontName','Times New Roman','FontSize',32,'color','k');%x轴坐标
xlabel('mm','FontName','Times New Roman','FontSize',32,'color','k');%x轴坐标
ylabel('mm','FontName','Times New Roman','FontSize',32,'color','k');%y轴坐标
set(gca,'FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',20);
title('Skin scattering distribution','FontName','Times New Roman','FontSize',28,'color','k'); %标题
set(fig15,'position',[250 250 1000 800]);
% print(fig15,'-dbitmap',strcat(dataFilePath,'\','Begin_End_Imagesc_DR.bmp'));
print(fig15,'-dbitmap',strcat(resultDir,'\',dataSubPath(iwb).name,'Begin_End_Imagesc_DR.bmp'));
savefig(fig15,strcat(resultDir,'\',dataSubPath(iwb).name,'Begin_End_Imagesc_DR.fig'),'compact');
% resultDir,dataSubPath(iwb).name
% %%衰减系数频数统计图
% %排除0.1噪声

for lasti = 1:2*r+1
    for lastj =1:2*r+1
        if  miu3D_mean2(lasti,lastj) == 0
            miu3D_mean2(lasti,lastj) = nan;
        end
    end
end
 fig6 = figure(6);clf;
h = histogram(miu3D_mean2,'Normalization','probability','binwidth',0.1);
set(gca,'FontSize',16);%,'FontSize',16
% set(gca,'LineWidth',1.5);
xlabel('um^{-1}','FontName','Times New Roman','FontSize',24,'color','k');%x轴坐标
ylabel('Frequency','FontName','Times New Roman','FontSize',24,'color','k');%y轴坐标
title('Depth frequency distribution','FontName','Times New Roman','FontSize',24,'color','k'); %标题
xlim([0,10]);
ylim([0,0.1]);
set(fig6,'position',[100 100 800 600]);
% print(fig6,'-dbitmap',strcat(dataFilePath,'\','Frequency_DR.bmp'));
print(fig6,'-dbitmap',strcat(resultDir,'\',dataSubPath(iwb).name,'Frequency_DR.bmp'));
% savefig(fig6,strcat(savePath,'\',name,'Frequency_DR.fig'),'compact');
savefig(fig6,strcat(resultDir,'\',dataSubPath(iwb).name,'Frequency_DR.fig'),'compact');
%% 单列A-scan
 fig7 = figure(7);clf;
    th = data3D_1(:,:,201);
      % % %     加入灵敏因子计算
    z = detaz*(1:1024);
    s_z = exp(-z.^2/25);  %S(z)--灵敏度因子---考虑了灵敏度衰减所引起的影响，该方法还采用了成像深度Z和宽度dert的高斯模型，dert=5；
    s_z_fix = s_z(1,rangeMin:rangeMax)';  %竖着一列 rangeMaxFix-rangeMin+1 个像素点
    miuOrigin = th./s_z_fix;
    IntensityOrigin = miuOrigin(:,201);
plot(IntensityOrigin);
print(fig7,'-dbitmap',strcat(resultDir,'\',dataSubPath(iwb).name,'Acan.bmp'));
miu3D286 = miu3D(:,:,286);

% % % 作图%
fig12 = figure(12);clf;
IntensityOrigin_depth = 400:1:700;
IntensityOrigin_depth = IntensityOrigin_depth*detaz;
plot(IntensityOrigin_depth(1,100:301),IntensityOrigin(100:301,1),'LineWidth',1.2,'color','k')
hold on;
plot(IntensityOrigin_depth(1,191:222),IntensityOrigin(191:222,1),'LineWidth',1.2,'color','r')
hold on;
xlabel('Depth/mm','FontName','Times New Roman','FontSize',14,'color','k');%x轴坐标
ylabel('OCT Signal (A.U.)','FontName','Times New Roman','FontSize',14,'color','k');%y轴坐标
axis([1.4 1.7 0 6000]) 
% title('功率对数信号','FontName','Times New Roman','FontSize',14,'color','k');
%% 存图
miu_Origin1 = pathX(:,:,201);
x1 = miu_Origin1(1,:);
y = 1:2*r+1;
x2 = miu_Origin1(2,:);
x3 = miu_Origin1(3,:);

miu_Origin2 = pathX(:,:,101);
y1 = miu_Origin2(1,:);
y2 = miu_Origin2(2,:);
y3 = miu_Origin2(3,:);

miu_Origin3 = pathX(:,:,301);
z1 = miu_Origin3(1,:);
z2 = miu_Origin3(2,:);
z3 = miu_Origin3(3,:);
%
 fig8 = figure(8);clf;
data3D1 = data3D(rangeMin:rangeMax,300:700,500);
data3D1 = imfilter(data3D1,fspecial('gaussian',[3 3],3));
imagesc(data3D1); axis image; colormap('gray'); hold on;
% print(fig8,'-dbitmap',strcat(dataFilePath,'\','fenge1.bmp'));
print(fig8,'-dbitmap',strcat(resultDir,'\',dataSubPath(iwb).name,'fenge1.bmp'));
 fig9 = figure(9);clf;
data3D1 = data3D(rangeMin:rangeMax,300:700,500);
data3D1 = imfilter(data3D1,fspecial('gaussian',[3 3],3));
imagesc(data3D1); axis image; colormap('gray'); hold on;
plot(y,x1,'r-','linewidth',1); hold on;
plot(y,x2,'g-','linewidth',1); hold on;
plot(y,x3,'b-','linewidth',1); hold on;
% print(fig9,'-dbitmap',strcat(dataFilePath,'\','fenge2.bmp'));
print(fig9,'-dbitmap',strcat(resultDir,'\',dataSubPath(iwb).name,'fenge2.bmp'));
 fig10 = figure(10);clf;
data3D_2 = data3D(rangeMin:rangeMax,300:700,400);
data3D_2 = imfilter(data3D_2,fspecial('gaussian',[3 3],3));
imagesc(data3D_2); axis image; colormap('gray'); hold on;
plot(y,y1,'r-','linewidth',1); hold on;
plot(y,y2,'g-','linewidth',1); hold on;
plot(y,y3,'b-','linewidth',1); hold on;
% print(fig10,'-dbitmap',strcat(dataFilePath,'\','fenge3.bmp'));
print(fig10,'-dbitmap',strcat(resultDir,'\',dataSubPath(iwb).name,'fenge3.bmp'));
 fig11 = figure(11);clf;
data3D_3 = data3D(rangeMin:rangeMax,300:700,600);
data3D_3 = imfilter(data3D_3,fspecial('gaussian',[3 3],3));
imagesc(data3D_3); axis image; colormap('gray'); hold on;
plot(y,z1,'r-','linewidth',1); hold on;
plot(y,z2,'g-','linewidth',1); hold on;
plot(y,z3,'b-','linewidth',1); hold on;
% print(fig11,'-dbitmap',strcat(dataFilePath,'\','fenge4.bmp'));
print(fig11,'-dbitmap',strcat(resultDir,'\',dataSubPath(iwb).name,'fenge4.bmp'));
close all

    
    progress = (iwb-1)/(length(dataSubPath)-1);
    waitbar(progress,wb1,['处理进度...' num2str(100*progress) '%']);
    disp(['处理进度...' num2str(100*progress) '%']);
    end
end
toc
close(wb1);
% close all

