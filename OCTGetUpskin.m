function [pathX,pathY] = OCTGetUpskin( imgNew,szImgNew,IntensityImg )
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
wv = 0.1;
% wv = 2;
%minimum weight 图中最小权重，稳定系统
minWeight = 1E-5;

%arry to store weights
adjMW = nan([numel(imgNew(:)),7]);

%arry to store point A locations
adjMX = nan([numel(imgNew(:)),7]);

%arry to store point B locations
adjMY = nan([numel(imgNew(:)),7]);

neighborIter = [1 0 -1 2 -2 1 -1;...
                1 1 1 1 1 0 0];

%fill in the above arrays according to Section 3.2
szadjMW = size(adjMW);
ind = 1; indR = 0;
while ind ~= szadjMW(1)*szadjMW(2) %this step can be made more efficient to increase speed.
    [i, j] = ind2sub(szadjMW,ind);    
    [iX,iY] = ind2sub(szImgNew,i);    
    jX = iX + neighborIter(1,j);
    jY = iY + neighborIter(2,j);
    if j == 2
        if jX >=1 && jX <= szImgNew(1) && jY >=1 && jY <= szImgNew(2)
         %save weight
         % set to minimum if on the sides
         if jY == 1 || jY == szImgNew(2)
            adjMW(i,j) = minWeight;
            
         % else, calculate the actual weight based on equation 1.
         else
            adjMW(i,j) = 2 - IntensityImg(iX,iY) - IntensityImg(jX,jY) - abs(IntensityImg(iX,iY)-IntensityImg(jX,jY)) + minWeight;     
         end
        %save the subscript of the corresponding nodes
        adjMX(i,j) = sub2ind(szImgNew,iX,iY);%%a的位置
        adjMY(i,j) = sub2ind(szImgNew,jX,jY);%%b的位置
        end
    else
     if jX >=1 && jX <= szImgNew(1) && jY >=1 && jY <= szImgNew(2)
         %save weight
         % set to minimum if on the sides
         if jY == 1 || jY == szImgNew(2)
            adjMW(i,j) = minWeight;
            
         % else, calculate the actual weight based on equation 1.
         else
            adjMW(i,j) = 2 - IntensityImg(iX,iY) - IntensityImg(jX,jY) - abs(IntensityImg(iX,iY)-IntensityImg(jX,jY)) +wv + minWeight;     
         end
        %save the subscript of the corresponding nodes
        adjMX(i,j) = sub2ind(szImgNew,iX,iY);%%a的位置
        adjMY(i,j) = sub2ind(szImgNew,jX,jY);%%b的位置
     end
    end
    ind = ind+1;
    
    %display progress
    if indR < round(10*ind/szadjMW(1)/szadjMW(2))
        indR = round(10*ind/szadjMW(1)/szadjMW(2));
%         display(sprintf('progress: %1.0f%% done, this may take a while...\n',100*indR/10));
    end
    
end

%assemble the adjacency matrix 

keepInd = ~isnan(adjMW(:)) & ~isnan(adjMX(:)) & ~isnan(adjMY(:));
adjMW = adjMW(keepInd);
adjMX = adjMX(keepInd);
adjMY = adjMY(keepInd);

%sparse matrices, based on eq 1 with the gradient,
adjMatrixW = sparse(adjMX(:),adjMY(:),adjMW(:),numel(imgNew(:)),numel(imgNew(:)));

%% get shortest path 上表面
%获得上表面
[ dist,path{1} ] = graphshortestpath( adjMatrixW, 1, numel(imgNew(:)) );
[pathX,pathY] = ind2sub(szImgNew,path{1});
%去掉图像边界的前几点和后几点
pathX =pathX(gradient(pathY)~=0);
pathY =pathY(gradient(pathY)~=0);
end

