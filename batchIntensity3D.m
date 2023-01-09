clear all;

%****************‘***************************************

%功能：解压缩3D数据
%完成度：已完
%码农：Hu Jie/Yan Li
%时间：------
%Matlab版本：2019a

%*******************************************************
fastMode = 1;

if fastMode
   tempPath = 'D:\'; 
end

octFilePath = uigetdir('Select 3D data path');      %选择数据文件夹
currentPWD = pwd;
tic
intensityPath = [octFilePath '\Intensity3D\'];
if isequal(octFilePath,0)    % 判断是否正常选择文件夹
   disp('User selected Cancel') 
else
    disp(['User selected ', octFilePath])
end

intensityDataNames = string(missing);

mkdir(intensityPath)
octSubPath = dir(octFilePath);
cd('C:\WinRAR')
for i = 1:length(octSubPath)
    if(isequal(octSubPath(i).name,'.')||... % 去除系统自带的两个隐文件夹
        isequal(octSubPath(i).name,'..')||...
        ~octSubPath(i).isdir) % 去除遍历中不是文件夹的
        continue;
    end
    
    octDir = dir([octFilePath '\' octSubPath(i).name '\*3D.oct']);
    for j =1:length(octDir) % 遍历所有文件
        fileName = [octFilePath '\' octSubPath(i).name '\' octDir(j).name];
        zipFileName = strrep(fileName,'.oct','.rar');
        unzipFilePath = [intensityPath octSubPath(i).name '\'];
        unzipFileName = [unzipFilePath 'Intensity.data'];
        unzipFileNameRename = strrep(unzipFileName,'Intensity.data',strrep(octDir(j).name,'.oct','.data'));
        
        disp(fileName);
        
        movefile(fileName,zipFileName);
        if fastMode
            mkdir(unzipFilePath);
            [status, results] = dos(['winrar e -ibck' ' ' zipFileName  ' ' 'data\Intensity.data' ' '  tempPath]);  %解压文件  
            movefile([tempPath 'Intensity.data'],unzipFileName);
        else
            [status, results] = dos(['winrar e -ibck' ' ' zipFileName  ' ' 'data\Intensity.data' ' '  unzipFilePath]);  %解压文件  
        end
        disp(['unzip status = ' num2str(status)]);
        movefile(zipFileName,fileName);
        movefile(unzipFileName,unzipFileNameRename);
        if ismissing(intensityDataNames)    %记录解压出的数据文件路径
            intensityDataNames(1) = unzipFileNameRename;
        else
            intensityDataNames = [intensityDataNames unzipFileNameRename];
        end
    end
end

cd(currentPWD);
toc