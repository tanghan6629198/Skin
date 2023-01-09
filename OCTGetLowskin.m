function [pathX2,pathY2] = OCTGetLowskin( pathX1,imgNew,szImgNew,IntensityImg2 )
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明
wv = 0.6;
for i = 1:size(imgNew,2)
%     pathX(1,i) = a;
    IntensityImg2(pathX1(1,i),i) = (1E-5);
    IntensityImg2(pathX1(1,i)-1,i) = (1E-5);
    IntensityImg2(pathX1(1,i)-2,i) = (1E-5);
    IntensityImg2(pathX1(1,i)-3,i) = (1E-5);
    IntensityImg2(pathX1(1,i)-4,i) = (1E-5);
    IntensityImg2(pathX1(1,i)-5,i) = (1E-5);  
    IntensityImg2(pathX1(1,i)-6,i) = (1E-5);%%皮肤较厚时
    IntensityImg2(pathX1(1,i)-7,i) = (1E-5);     %   
    IntensityImg2(pathX1(1,i)-8,i) = (1E-5);     %   
    IntensityImg2(pathX1(1,i)-9,i) = (1E-5);     %   
    IntensityImg2(pathX1(1,i)-10,i) = (1E-5);    %   
    IntensityImg2(pathX1(1,i)-11,i) = (1E-5);    %   
    IntensityImg2(pathX1(1,i)-12,i) = (1E-5);    %   
    IntensityImg2(pathX1(1,i)-13,i) = (1E-5);    %   
    IntensityImg2(pathX1(1,i)-14,i) = (1E-5);    %   
    IntensityImg2(pathX1(1,i)-15,i) = (1E-5);%%皮肤较厚时
     
    IntensityImg2(pathX1(1,i)+1,i) = (1E-5);
    IntensityImg2(pathX1(1,i)+2,i) = (1E-5);
    IntensityImg2(pathX1(1,i)+3,i) = (1E-5);
    IntensityImg2(pathX1(1,i)+4,i) = (1E-5);
    IntensityImg2(pathX1(1,i)+5,i) = (1E-5);  
    IntensityImg2(pathX1(1,i)+6,i) = (1E-5);%%皮肤较厚时
    IntensityImg2(pathX1(1,i)+7,i) = (1E-5);     % 
    IntensityImg2(pathX1(1,i)+8,i) = (1E-5);     %
    IntensityImg2(pathX1(1,i)+9,i) = (1E-5);     %
    IntensityImg2(pathX1(1,i)+10,i) = (1E-5);    %
    IntensityImg2(pathX1(1,i)+11,i) = (1E-5);    %
    IntensityImg2(pathX1(1,i)+12,i) = (1E-5);    %
    IntensityImg2(pathX1(1,i)+13,i) = (1E-5);    %
    IntensityImg2(pathX1(1,i)+14,i) = (1E-5);    %
    IntensityImg2(pathX1(1,i)+15,i) = (1E-5);%%皮肤较厚时
end

%minimum weight 图中最小权重，稳定系统
minWeight = 1E-5;

%arry to store weights
adjMW2 = nan([numel(imgNew(:)),7]);

%arry to store point A locations
adjMX2 = nan([numel(imgNew(:)),7]);

%arry to store point B locations
adjMY2 = nan([numel(imgNew(:)),7]);

neighborIter = [1 0 -1 2 -2 1 -1;...
                1 1 1 1 1 0 0];
                     
szadjMW2 = size(adjMW2);
ind2 = 1; indR2 = 0;
while ind2 ~= szadjMW2(1)*szadjMW2(2) %this step can be made more efficient to increase speed.
    [i2, j2] = ind2sub(szadjMW2,ind2);    
    [iX2,iY2] = ind2sub(szImgNew,i2);    
    jX2 = iX2 + neighborIter(1,j2);
    jY2  = iY2 + neighborIter(2,j2);
    if j2 == 2
       if jX2 >=1 && jX2 <= szImgNew(1) && jY2 >=1 && jY2 <= szImgNew(2)
         %save weight
         % set to minimum if on the sides
         if jY2 == 1 || jY2 == szImgNew(2)
            adjMW2(i2,j2) = minWeight;
            
         % else, calculate the actual weight based on equation 1.
         else
            adjMW2(i2,j2) = 2 - IntensityImg2(iX2,iY2) - IntensityImg2(jX2,jY2) - abs(IntensityImg2(iX2,iY2) - IntensityImg2(jX2,jY2)) + minWeight;     
         end
        %save the subscript of the corresponding nodes
        adjMX2(i2,j2) = sub2ind(szImgNew,iX2,iY2);%%a的位置
        adjMY2(i2,j2) = sub2ind(szImgNew,jX2,jY2);%%b的位置
       end
    else
     if jX2 >=1 && jX2 <= szImgNew(1) && jY2 >=1 && jY2 <= szImgNew(2)
         %save weight
         % set to minimum if on the sides
         if jY2 == 1 || jY2 == szImgNew(2)
            adjMW2(i2,j2) = minWeight;
            
         % else, calculate the actual weight based on equation 1.
         else
            adjMW2(i2,j2) = 2 - IntensityImg2(iX2,iY2) - IntensityImg2(jX2,jY2) - abs(IntensityImg2(iX2,iY2) - IntensityImg2(jX2,jY2)) +wv + minWeight;     
         end
        %save the subscript of the corresponding nodes
        adjMX2(i2,j2) = sub2ind(szImgNew,iX2,iY2);%%a的位置
        adjMY2(i2,j2) = sub2ind(szImgNew,jX2,jY2);%%b的位置
     end
    end
    ind2 = ind2+1;
    
    %display progress
    if indR2 < round(10*ind2/szadjMW2(1)/szadjMW2(2))
        indR2 = round(10*ind2/szadjMW2(1)/szadjMW2(2));
%         display(sprintf('progress: %1.0f%% done, this may take a while...\n',100*indR2/10));
    end
    
end

%assemble the adjacency matrix 

keepInd2 = ~isnan(adjMW2(:)) & ~isnan(adjMX2(:)) & ~isnan(adjMY2(:));
adjMW2 = adjMW2(keepInd2);
adjMX2 = adjMX2(keepInd2);
adjMY2 = adjMY2(keepInd2);

%sparse matrices, based on eq 1 with the gradient,
adjMatrixW2 = sparse(adjMX2(:),adjMY2(:),adjMW2(:),numel(imgNew(:)),numel(imgNew(:)));

%% get shortest path 

% get layer 
[ dist2,path2{1} ] = graphshortestpath( adjMatrixW2, 1, numel(imgNew(:)) );
[pathX2,pathY2] = ind2sub(szImgNew,path2{1});

% get rid of first and last few points that is by the image borders
pathX2 =pathX2(gradient(pathY2)~=0);
pathY2 =pathY2(gradient(pathY2)~=0);

end

