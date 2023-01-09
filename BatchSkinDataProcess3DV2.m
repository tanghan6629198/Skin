 


clc
clear 

%*******************************************************

%���ܣ����г�����Ż����� ��ɶȣ������ ��ũ��Zhao Ruihang ʱ�䣺2020.04.04 Matlab�汾��2019a

%*******************************************************

 %% ����3D���ݴ���

dataFilePath =  uigetdir('Select unziped 2D data path');     %ѡ�������ļ���
tic
if isequal(dataFilePath,0)    % �ж��Ƿ�����ѡ���ļ���
   disp('User selected Cancel') ;
else
    disp(['User selected ', dataFilePath]);
end

resultDir = strcat(dataFilePath,'\','result');
if exist(resultDir,'dir') == 0
    mkdir(resultDir);
end

dataSubPath = dir(dataFilePath);
wb = waitbar(0,'�������...');
toc

for i = 1:length(dataSubPath)
    
    
    if(isequal(dataSubPath(i).name,'.')||... % ȥ��ϵͳ�Դ����������ļ���
        isequal(dataSubPath(i).name,'..')||...
        ~dataSubPath(i).isdir) % ȥ�������в����ļ��е�
        continue;
    end
    dataDir = dir([dataFilePath '\' dataSubPath(i).name '\*.data']);
    for j =1:length(dataDir) % ���������ļ�
        fileName = [dataFilePath '\' dataSubPath(i).name '\' dataDir(j).name]
        disp(fileName);

%���м���Ƥ��������ȡ��ֲڶȵĳ���
          [surfaceLocation,bottomLocation,depth,r] = skinDataProcessV4(fileName); % ֻ�����±������ڷ�ֵȡ����Ż��㷨 
%             [surfaceLocation,depth,r] = skinDataProcessV2(fileName);
%           [surfaceLocation,bottomLocation,surfaceIntensityA,bottomIntensityA,depth,r] = skinDataProcessV5(fileName); % ֻ���б߽�ʶ����Ż��㷨 
%           [surfaceLocation,bottomLocation,surfaceIntensityA,bottomIntensityA,depth,r] = skinDataProcessV2(fileName); % �߽��Ż�+���·�ֵ����
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
%                saveSkinDataImageV8(resultDir,dataSubPath(i).name,dataFilePath);  % ���Ĵֲڶȵļ���
    
    progress = (i-1)/(length(dataSubPath)-1);
    waitbar(progress,wb,['�������...' num2str(100*progress) '%']);
    disp(['�������...' num2str(100*progress) '%']);
    end
end
toc
close(wb);
% close all