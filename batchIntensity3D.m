clear all;

%****************��***************************************

%���ܣ���ѹ��3D����
%��ɶȣ�����
%��ũ��Hu Jie/Yan Li
%ʱ�䣺------
%Matlab�汾��2019a

%*******************************************************
fastMode = 1;

if fastMode
   tempPath = 'D:\'; 
end

octFilePath = uigetdir('Select 3D data path');      %ѡ�������ļ���
currentPWD = pwd;
tic
intensityPath = [octFilePath '\Intensity3D\'];
if isequal(octFilePath,0)    % �ж��Ƿ�����ѡ���ļ���
   disp('User selected Cancel') 
else
    disp(['User selected ', octFilePath])
end

intensityDataNames = string(missing);

mkdir(intensityPath)
octSubPath = dir(octFilePath);
cd('C:\WinRAR')
for i = 1:length(octSubPath)
    if(isequal(octSubPath(i).name,'.')||... % ȥ��ϵͳ�Դ����������ļ���
        isequal(octSubPath(i).name,'..')||...
        ~octSubPath(i).isdir) % ȥ�������в����ļ��е�
        continue;
    end
    
    octDir = dir([octFilePath '\' octSubPath(i).name '\*3D.oct']);
    for j =1:length(octDir) % ���������ļ�
        fileName = [octFilePath '\' octSubPath(i).name '\' octDir(j).name];
        zipFileName = strrep(fileName,'.oct','.rar');
        unzipFilePath = [intensityPath octSubPath(i).name '\'];
        unzipFileName = [unzipFilePath 'Intensity.data'];
        unzipFileNameRename = strrep(unzipFileName,'Intensity.data',strrep(octDir(j).name,'.oct','.data'));
        
        disp(fileName);
        
        movefile(fileName,zipFileName);
        if fastMode
            mkdir(unzipFilePath);
            [status, results] = dos(['winrar e -ibck' ' ' zipFileName  ' ' 'data\Intensity.data' ' '  tempPath]);  %��ѹ�ļ�  
            movefile([tempPath 'Intensity.data'],unzipFileName);
        else
            [status, results] = dos(['winrar e -ibck' ' ' zipFileName  ' ' 'data\Intensity.data' ' '  unzipFilePath]);  %��ѹ�ļ�  
        end
        disp(['unzip status = ' num2str(status)]);
        movefile(zipFileName,fileName);
        movefile(unzipFileName,unzipFileNameRename);
        if ismissing(intensityDataNames)    %��¼��ѹ���������ļ�·��
            intensityDataNames(1) = unzipFileNameRename;
        else
            intensityDataNames = [intensityDataNames unzipFileNameRename];
        end
    end
end

cd(currentPWD);
toc