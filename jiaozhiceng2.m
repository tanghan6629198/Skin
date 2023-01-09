%% 1
imgI = 0;
imgIbegin = 25;
%拉平
    for i = 1:100
        begin = 35;
        depth = max(pathX2(1,begin+i),pathX(1,begin+i)) - pathX3(1,begin+i) ;
        imgI(50-depth:50,i) = img3(pathX3(1,begin+i):max(pathX2(1,begin+i),pathX(1,begin+i)),begin+i);
        imgI = imgI(1:48,:);
    end
    figure();imagesc(imgI); axis image; colormap('gray'); hold on;
    
    I = imgI(imgIbegin:48,:);
    Offset = [0,3];  %0° d=3
    N = 8; %灰度等级
[ASM1,ENT1,COR1] = graycomatrix(Offset,I,N);
%% 2
imgI = 0;
%拉平
    for i = 1:100
        begin = 135;
        depth = max(pathX2(1,begin+i),pathX(1,begin+i)) - pathX3(1,begin+i) ;
        imgI(50-depth:50,i) = img3(pathX3(1,begin+i):max(pathX2(1,begin+i),pathX(1,begin+i)),begin+i);
        imgI = imgI(1:48,:);
    end
    figure();imagesc(imgI); axis image; colormap('gray'); hold on;
    
    I = imgI(imgIbegin:48,:);
    Offset = [0,3];  %0° d=3
    N = 8; %灰度等级
[ASM2,ENT2,COR2] = graycomatrix(Offset,I,N);
%% 3
imgI = 0;
%拉平
    for i = 1:100
        begin = 235;
        depth = max(pathX2(1,begin+i),pathX(1,begin+i)) - pathX3(1,begin+i) ;
        imgI(50-depth:50,i) = img3(pathX3(1,begin+i):max(pathX2(1,begin+i),pathX(1,begin+i)),begin+i);
        imgI = imgI(1:48,:);
    end
    figure();imagesc(imgI); axis image; colormap('gray'); hold on;
    
    I = imgI(imgIbegin:48,:);
    Offset = [0,3];  %0° d=3
    N = 8; %灰度等级
[ASM3,ENT3,COR3] = graycomatrix(Offset,I,N);
%% 4
imgI = 0;
%拉平
    for i = 1:100
        begin = 335;
        depth = max(pathX2(1,begin+i),pathX(1,begin+i)) - pathX3(1,begin+i) ;
        imgI(50-depth:50,i) = img3(pathX3(1,begin+i):max(pathX2(1,begin+i),pathX(1,begin+i)),begin+i);
        imgI = imgI(1:48,:);
    end
    figure();imagesc(imgI); axis image; colormap('gray'); hold on;
    
    I = imgI(imgIbegin:48,:);
    Offset = [0,3];  %0° d=3
    N = 8; %灰度等级
[ASM4,ENT4,COR4] = graycomatrix(Offset,I,N);
%% 5
imgI = 0;
%拉平
    for i = 1:100
        begin = 435;
        depth = max(pathX2(1,begin+i),pathX(1,begin+i)) - pathX3(1,begin+i) ;
        imgI(50-depth:50,i) = img3(pathX3(1,begin+i):max(pathX2(1,begin+i),pathX(1,begin+i)),begin+i);
        imgI = imgI(1:48,:);
    end
    figure();imagesc(imgI); axis image; colormap('gray'); hold on;
    
    I = imgI(imgIbegin:48,:);
    Offset = [0,3];  %0° d=3
    N = 8; %灰度等级
[ASM5,ENT5,COR5] = graycomatrix(Offset,I,N);

ASM = (ASM5+ASM4+ASM3+ASM2+ASM1)/5
ENT = (ENT5+ENT4+ENT3+ENT2+ENT1)/5;
COR = (COR5+COR4+COR3+COR2+COR1)/5