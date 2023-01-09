 


clc
clear 

%*******************************************************

%功能：进行程序的优化处理 完成度：待完成 码农：Zhao Ruihang 时间：2020.04.04 Matlab版本：2019a

%*******************************************************

 %% 进行3D数据处理

dataFilePath =  uigetdir('Select unziped 2D data path');     %选择数据文件夹
tic
if isequal(dataFilePath,0)    % 判断是否正常选择文件夹
   disp('User selected Cancel') ;
else
    disp(['User selected ', dataFilePath]);
end

resultDir = strcat(dataFilePath,'\','result');
if exist(resultDir,'dir') == 0
    mkdir(resultDir);
end

dataSubPath = dir(dataFilePath);
wb = waitbar(0,'处理进度...');
toc

for i = 1:length(dataSubPath)
    
    
    if(isequal(dataSubPath(i).name,'.')||... % 去除系统自带的两个隐文件夹
        isequal(dataSubPath(i).name,'..')||...
        ~dataSubPath(i).isdir) % 去除遍历中不是文件夹的
        continue;
    end
    dataDir = dir([dataFilePath '\' dataSubPath(i).name '\*.data']);
    for j =1:length(dataDir) % 遍历所有文件
        fileName = [dataFilePath '\' dataSubPath(i).name '\' dataDir(j).name]
        disp(fileName);

%进行计算皮肤样本厚度、粗糙度的程序
          [surfaceLocation,bottomLocation,depth,r] = skinDataProcessV4(fileName); % 只进行下表面相邻峰值取舍的优化算法 
%             [surfaceLocation,depth,r] = skinDataProcessV2(fileName);
%           [surfaceLocation,bottomLocation,surfaceIntensityA,bottomIntensityA,depth,r] = skinDataProcessV5(fileName); % 只进行边界识别的优化算法 
%           [surfaceLocation,bottomLocation,surfaceIntensityA,bottomIntensityA,depth,r] = skinDataProcessV2(fileName); % 边界优化+上下峰值计算
          save(strcat(dataFilePath,'\',dataSubPath(i).name,'\','surfaceLocation','.mat'),'surfaceLocation');
          save(strcat(dataFilePath,'\',dataSubPath(i).name,'\','bottomLocation','.mat'),'bottomLocation');
          save(strcat(dataFilePath,'\',dataSubPath(i).name,'\','depth','.mat'),'depth');
          save(strcat(dataFilePath,'\',dataSubPath(i).name,'\','r','.mat'),'r');
          surfaceLocation=load([dataFilePath '\' dataSubPath(i).name '\' 'surfaceLocation.mat']);
          bottomLocation=load([dataFilePath '\' dataSubPath(i).name '\' 'bottomLocation.mat']);
          depth=load([dataFilePath '\' dataSubPath(i).name '\' 'depth.mat']);
          r=load([dataFilePath '\' dataSubPath(i).name '\' 'r.mat']);

%           [filtSurface,SurfaceTrue] = saveSkinDataImageV3(surfaceLocation.surfaceLocation,bottomLocation.bottomLocation,depth.depth,r.r,resultDir,dataSubPath(i).name);
          [filtSurface,SurfaceTrue] = saveSkinDataImageV2_5(surfaceLocation.surfaceLocation,bottomLocation.bottomLocation,depth.depth,r.r,resultDir,dataSubPath(i).name);
%             [filtSurface,SurfaceTrue] = saveSkinDataImageV2_2(surfaceLocation.surfaceLocation,bottomLocation.bottomLocation,depth.depth,r.r,resultDir,dataSubPath(i).name);
           save(strcat(dataFilePath,'\',dataSubPath(i).name,'\','SurfaceTrue','.mat'),'SurfaceTrue');
%             saveSkinDataImageV3(surfaceLocation.surfaceLocation,depth.depth,r.r,resultDir,dataSubPath(i).name);
%                saveSkinDataImageV8(resultDir,dataSubPath(i).name,dataFilePath);  % 论文粗糙度的计算
    
    progress = (i-1)/(length(dataSubPath)-1);
    waitbar(progress,wb,['处理进度...' num2str(100*progress) '%']);
    disp(['处理进度...' num2str(100*progress) '%']);
    end
end
toc
close(wb);
% close all