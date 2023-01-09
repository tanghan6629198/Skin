function [pathX3,pathY3] = OCTGetCuticle(pathX,pathX2,imgNew,szImgNew,gradImg)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
% T = 8;%DAY13    角质层阈值
% T = 3; %DAY5
T =8; %最后1天
wv = 0.2;
for i = 1:size(imgNew,2)
    if pathX(1,1)<pathX2(1,1)
        gradImg(1:pathX(1,i),i) = (1E-5);
        gradImg(pathX2(1,i):szImgNew(1),i) = (1E-5);
        gradImg(pathX(1,i):pathX(1,i)+T,i) = (1E-5);
%         gradImg(pathX(1,i)+4,i) = (1E-5);
%         gradImg(pathX(1,i)+5,i) = (1E-5);
%         
        gradImg(pathX2(1,i)-1,i) = (1E-5);
        gradImg(pathX2(1,i)-2,i) = (1E-5);
        gradImg(pathX2(1,i)-3,i) = (1E-5);
        gradImg(pathX2(1,i)-4,i) = (1E-5);
%         gradImg(pathX2(1,i)-5,i) = (1E-5);
%         gradImg(pathX2(1,i)-6,i) = (1E-5);%%皮肤较厚时
%         gradImg(pathX2(1,i)-7,i) = (1E-5);
%         gradImg(pathX2(1,i)-8,i) = (1E-5);
%         gradImg(pathX2(1,i)-9,i) = (1E-5);
%         gradImg(pathX2(1,i)-10,i) = (1E-5);%%皮肤较厚时
    else
        gradImg(1:pathX2(1,i),i) = (1E-5);
        gradImg(pathX(1,i):szImgNew(1),i) = (1E-5);
        gradImg(pathX2(1,i):pathX2(1,i)+T,i) = (1E-5);
%         gradImg(pathX2(1,i)+4,i) = (1E-5);
%         gradImg(pathX2(1,i)+5,i) = (1E-5);
        
        gradImg(pathX(1,i)-1,i) = (1E-5);
        gradImg(pathX(1,i)-2,i) = (1E-5);
        gradImg(pathX(1,i)-3,i) = (1E-5);
        gradImg(pathX(1,i)-4,i) = (1E-5);
        gradImg(pathX(1,i)-5,i) = (1E-5);
%         gradImg(pathX(1,i)-6,i) = (1E-5);%%皮肤较厚时
%         gradImg(pathX(1,i)-7,i) = (1E-5);
%         gradImg(pathX(1,i)-8,i) = (1E-5);
%         gradImg(pathX(1,i)-9,i) = (1E-5);
%         gradImg(pathX(1,i)-10,i) = (1E-5);%%皮肤较厚时
    end

end

%minimum weight 图中最小权重，稳定系统
minWeight = 1E-5;

%arry to store weights
adjMW3 = nan([numel(imgNew(:)),7]);

%arry to store point A locations
adjMX3 = nan([numel(imgNew(:)),7]);

%arry to store point B locations
adjMY3 = nan([numel(imgNew(:)),7]);

neighborIter = [1 0 -1 2 -2 1 -1;...
                1 1 1 1 1 0 0];
                   
szadjMW3 = size(adjMW3);
ind3 = 1; indR3 = 0;
while ind3 ~= szadjMW3(1)*szadjMW3(2) %this step can be made more efficient to increase speed.
    [i3, j3] = ind2sub(szadjMW3,ind3);    
    [iX3,iY3] = ind2sub(szImgNew,i3);    
    jX3 = iX3 + neighborIter(1,j3);
    jY3  = iY3 + neighborIter(2,j3);
    if j3 == 2
        if jX3 >=1 && jX3 <= szImgNew(1) && jY3 >=1 && jY3 <= szImgNew(2)
         %save weight
         % set to minimum if on the sides
         if jY3 == 1 || jY3 == szImgNew(2)
            adjMW3(i3,j3) = minWeight;
            
         % else, calculate the actual weight based on equation 1.
         else
            adjMW3(i3,j3) = 2 - gradImg(iX3,iY3) - gradImg(jX3,jY3) + minWeight;     
         end
        %save the subscript of the corresponding nodes
        adjMX3(i3,j3) = sub2ind(szImgNew,iX3,iY3);%%a的位置
        adjMY3(i3,j3) = sub2ind(szImgNew,jX3,jY3);%%b的位置
       end        
   else 
     if jX3 >=1 && jX3 <= szImgNew(1) && jY3 >=1 && jY3 <= szImgNew(2)
         %save weight
         % set to minimum if on the sides
         if jY3 == 1 || jY3 == szImgNew(2)
            adjMW3(i3,j3) = minWeight;
            
         % else, calculate the actual weight based on equation 1.
         else
            adjMW3(i3,j3) = 2 - gradImg(iX3,iY3) - gradImg(jX3,jY3) +wv + minWeight;     
         end
        %save the subscript of the corresponding nodes
        adjMX3(i3,j3) = sub2ind(szImgNew,iX3,iY3);%%a的位置
        adjMY3(i3,j3) = sub2ind(szImgNew,jX3,jY3);%%b的位置
     end
    end
    ind3 = ind3+1;
    
    %display progress
    if indR3 < round(10*ind3/szadjMW3(1)/szadjMW3(2))
        indR3 = round(10*ind3/szadjMW3(1)/szadjMW3(2));
%         display(sprintf('progress: %1.0f%% done, this may take a while...\n',100*indR3/10));
    end
    
end

%assemble the adjacency matrix 

keepInd3 = ~isnan(adjMW3(:)) & ~isnan(adjMX3(:)) & ~isnan(adjMY3(:));
adjMW3 = adjMW3(keepInd3);
adjMX3 = adjMX3(keepInd3);
adjMY3 = adjMY3(keepInd3);

%sparse matrices, based on eq 1 with the gradient,
adjMatrixW3 = sparse(adjMX3(:),adjMY3(:),adjMW3(:),numel(imgNew(:)),numel(imgNew(:)));

%% get shortest path 

% get layer 
[ dist3,path3{1} ] = graphshortestpath( adjMatrixW3, 1, numel(imgNew(:)) );
[pathX3,pathY3] = ind2sub(szImgNew,path3{1});

% get rid of first and last few points that is by the image borders
pathX3 =pathX3(gradient(pathY3)~=0);
pathY3 =pathY3(gradient(pathY3)~=0);
end

