    function [surfaceLocation,bottomLocation,depth,r] = skinDataProcessV4( filePath )
% function [depth_Line,surfaceLine,bottomLine] = skinDataProcessV2( filePath )
%*******************************************************
%不同皮肤样本参数修改： A_scan  B_scan  line   Depth
%功能：进行算法优化后的比较
%完成度：表面检测算法
%码农：Zhao Ruihang
%时间：2020.05.07

%*******************************************************

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% 基本参数设置
% 最新的参数
% A_scan=1024;  % 采样点数
% B_scan=1024;   % the times of A-scan in per B-scan   %A-Scan的数量
% line=1024;     % the times of B-scan    %B-Scan的数量S


A_scan=1024;  % 采样点数
B_scan=1000;   % the times of A-scan in per B-scan   %A-Scan的数量
line=1000;     % the times of B-scan    %B-Scan的数量
r_range = 285;
Depth = 18;      %%皮肤 厚度阈值 Day不同  要修改
%% 读取数据

shft = B_scan*A_scan;  %一个bscan总的采样点数
fid=fopen(filePath,'r','n');
wb = waitbar(0,'读取3D数据中...');
disp('读取3D数据中...');
data3D = zeros(A_scan,B_scan,line,'single');
A = zeros(A_scan,B_scan,'single')*NaN;
for i=1:line   
    A = fread(fid,[A_scan,B_scan],'float32','n'); 
    if isempty(A)==1    %表示逻辑运算符“非”，也就是取反
        break;  %当A是空矩阵时，跳出
    end
    fseek(fid, 4*i*shft, 'bof');  %指针移入第j*shft个元素
    data3D(:,:,i)=A;
    waitbar(i/line,wb,['读取3D数据中...' num2str(100*i/line) '%']);
end
clear A;
close(wb);
% for i = 1:line/2
%     for j = 1:B_scan/2
%         data3Dnew(:,j,i) = data3D(:,2*j-1,2*i-1);
%     end
% end
%% 裁剪有效区域
wb = waitbar(0,'裁剪有效区域...');
disp('裁剪有效区域...');

% 确定范围
section2D = zeros(A_scan,B_scan,'single')*NaN;
section2D(:) = data3D(:,:,line/2);
figure(1);
clf;
imagesc(section2D);

% answer = inputdlg({'RangeMin','RangeMax'},'Range');
% rangeMin = str2double(answer(1));
% rangeMax = str2double(answer(2));
rangeMin = 450;
rangeMax = 750;
close(figure(1));

range = rangeMin:rangeMax;  %皮肤信号范围

% Thrmo的边缘检测
% section90 = zeros(A_scan,B_scan,'single')*NaN;
% section180 = zeros(A_scan,B_scan,'single')*NaN;
section90 = data3D(:,:,line/2);
section180 = data3D(:,B_scan/2,:);
%% 边缘识别


centerX = B_scan/2;
centerY = line/2;

% r = 36; %半径0.25；直径 0.5
%  r = 71;  % 半径：0.5mm；直径为 1
% r = 398; %半径 2.5；直径 5
% r = 278;
r=r_range;     %如何跟据像素来进行范围的确定  %2.5mm
% r = 200;

r2 = r^2;
waitbar(1,wb,['裁剪有效区域...' num2str(20)  '%']);
nanLine = zeros(A_scan,1)*NaN('single');

for j = 1:B_scan
    for k = 1:line
        if((j-centerX)^2+(k-centerY)^2)>r2
          	data3D(:,j,k)=nanLine;
        end
    end
    waitbar(j/B_scan,wb,['裁剪有效区域...' num2str(20 + j/B_scan*80) '%']);
end

close(wb);

p1 = [centerX-r,centerX+r];   %x轴的坐标
p2 = [rangeMax,rangeMin];   %y轴的坐标

% %% 原始成图
% figure(1),clf;
% imagesc(section2D);
% hold on
% p1 = [range_begin90+20,range_end90-20];   %x轴的坐标
% p2 = [rangeMax,rangeMin];   %y轴的坐标
% plot([p1(1),p1(1)],[p2(1),p2(2)],'color','r','linewidth',2);   %起始端线条
% plot([p1(2),p1(2)],[p2(1),p2(2)],'color','r','linewidth',2);   %终止端线条
% hold off

% xlabel('mm','FontName','Times New Roman','FontSize',40,'color','k');%x轴坐标
% ylabel('mm','FontName','Times New Roman','FontSize',40,'color','k');%y轴坐标
% zlabel('Thickness/\mum','FontName','Times New Roman','FontSize',40,'color','k');%z轴坐标

%% 分割皮肤表面

% surfaceRelativeLocation = zeros(1,p2(2)-p2(1),'single')*NaN;
% bottomRelativeLocation = zeros(1,p2(2)-p2(1),'single')*NaN;
% surfaceIntensity = zeros(1,p2(2)-p2(1),'single')*NaN;
% bottomIntensity = zeros(1,p2(2)-p2(1),'single')*NaN;

% JudgeBig = zeros(1,p2(2)-p2(1),'single')*NaN;
% IntensityBig = zeros(1,p2(2)-p2(1),'single')*NaN;
% JudgeSmall = zeros(1,p2(2)-p2(1),'single')*NaN;
% IntensitySmall = zeros(1,p2(2)-p2(1),'single')*NaN;


surfaceLocation = zeros(2*r,p2(2)-p2(1),'single')*NaN;
bottomLocation = zeros(2*r,p2(2)-p2(1),'single')*NaN;
surfaceIntensityA = zeros(2*r,p2(2)-p2(1),'single')*NaN;
bottomIntensityA = zeros(2*r,p2(2)-p2(1),'single')*NaN;
depth = zeros(2*r,p2(2)-p2(1),'single')*NaN;

% surfacePercentageA = zeros(2*r,p2(2)-p2(1),'single')*NaN;
% bottomPercentageA = zeros(2*r,p2(2)-p2(1),'single')*NaN;

% indexAdd = range(1)-1;
wbSkin3dProcess = waitbar(100,'计算皮肤厚度...');
zeroline = zeros(rangeMax +1 - rangeMin,1);
IndexK = 0;
for k = centerY - r  : centerY  + r     %每个B-Scan依次进行处理     %选择B-scan，选择是哪条B-Scan
%     k;
    IndexK = IndexK + 1;
    Intensity = data3D(range,:,k); %提取的强度信号
    
    %先找到每条B-Scan的位置和数值大小
    %判断是不是极值点：二阶导数：>0那么就是极小值，反之为极大值
    %判断左右领域的导数值……左+右- 为极大值点。左-右+ 为极小值点。
    i_big = 1;
    i_small = 1;
    Index = 0;
    for j = centerX - r : centerX  + r     %每个A-Scan依次处理   %选择B-Scan的哪条A-scan
%         j;
        Index = Index + 1;
        Data = Intensity(:,j);
        DataDiff = diff(Intensity(:,j));
        
%         %A-Scan成像
%         figure(k),clf;
%         plot(Data);
        
        if (~isnan(Intensity(1,j))&&all(Intensity(:,j)~=zeroline))
%             problem = 1;  %初期程序进行排查的标志位
%每条A-scan重新置0，防止干扰            
            JudgeBig = zeros(1,p2(2)-p2(1),'single') ;
            JudgeSmall  = zeros(1,p2(2)-p2(1),'single');
            IntensityBig = zeros(1,p2(2)-p2(1),'single');
            IntensitySmall = zeros(1,p2(2)-p2(1),'single');
            IndexJudge = 0;
            m = 0;n = 0;
%             SeekBigIndex = 0;
            SeekBigIntensity = 0;
            maxPeak = 0;maxIndex = 0;
            max2Peak = 0;max2Index = 0;
            surfaceRelativeLocation = 0; 
            bottomRelativeLocation = 0;   
            surfaceIntensity = 0; 
            bottomIntensity = 0 ; 
            depthNum = 0;
%判断极大值和极小值的位置和强度值
            for rangeIndex = 1:rangeMax-rangeMin-1
                if DataDiff(rangeIndex)>0 && DataDiff(rangeIndex+1) <0   %极大值点
%                     IndexJudge = 1; %前期程序进行错误排查标志位
                    JudgeBig(i_big) = rangeIndex + 1;          %极大值的位置
                    IntensityBig(i_big) = Data(rangeIndex + 1);   %极大值的强度值
                    i_big = i_big + 1;
                else
                    if DataDiff(rangeIndex) <0 && DataDiff(rangeIndex+1) >0 %极小值点
%                         IndexJudge = 2; %前期程序进行错误排查的标志位
                        JudgeSmall(i_small) = rangeIndex + 1;  %极小值的位置
                        IntensitySmall(i_small) = Data(rangeIndex + 1);  %极小值的强度值
                        i_small = i_small +1;
                    end
                end
            end
%对A-Scan进行表面定位
%            [m,n] = size(JudgeBig);
%            SeekBigIndex = JudgeBig;         % 重新赋值极大值的位置，防止原始数据丢失--用处不大
           SeekBigIntensity = IntensityBig; % 重新赋值极大值的轻强度值，防止原始数据丢失
           [~,maxNumIndex] = max(IntensityBig); %获得最大值的强度值 以及  最大值是第几个峰值
           maxIndex = JudgeBig(maxNumIndex);          %获得最大值的相对所在位置  ：实际所在位置应该加上rangeMin
           SeekBigIntensity(maxNumIndex) = 0; %最大值的强度值赋为0，以查找下一个峰值

           
           flag = 1 ;  %来进入循环，以判断峰值的准确性
           while(flag)                
               %查找下一处峰值
               [~,max2NumIndex] = max(SeekBigIntensity); %获得次大值的强度值  以计  次大值是第几个峰值
               max2Index = JudgeBig(max2NumIndex);          %获得次大值的相对所在位置
               SeekBigIntensity(max2NumIndex) = 0;              %次大值的强度值赋为0
               
               
               Day2 = Depth;
               if abs(maxIndex - max2Index) <= Depth
                   Day2 = Depth - 5;   %皮肤厚度不均匀时使用
               end
               if abs(maxIndex - max2Index) > Day2  %10pixel = 25um.两个峰值不能小于25um--Day01的厚度--9Day为15--13Day为20
                   
                  %获得上表面和下表面的相对位置和对应的强度值 
                  surfaceRelativePeak = min(maxNumIndex,max2NumIndex);  %获得上表面的是第几个峰值
                  bottomRelativePeak = max(maxNumIndex,max2NumIndex); %获得下表面的是第几个峰值
                  surfaceRelativeLocation = JudgeBig(surfaceRelativePeak); %获取上表面峰值的相对位置
                  bottomRelativeLocation = JudgeBig(bottomRelativePeak);   %获取下表面峰值的相对位置
                  surfaceIntensity = IntensityBig(surfaceRelativePeak); %获取上表面的强度值
                  bottomIntensity = IntensityBig(bottomRelativePeak);   %获取下表面的强度值
                  
                  %对峰值前后进行自适应计算            
                  InspectRange =round( abs(surfaceRelativeLocation - bottomRelativeLocation) / 2); %检查距离 = 两个峰值的一半距离
%                     InspectRange = 5;
%上表面检查                 
                  %判断在检查距离内有哪几个峰值？
                  surfaceRangeMin = round(surfaceRelativeLocation - InspectRange);
                  surfaceRangeMax = round(surfaceRelativeLocation + InspectRange);
                  [~,n] = size(JudgeBig);
                  flagChange = 1;
                  for surfaceRangeIndex = 1:n
                      if flagChange == 1
                          if JudgeBig(surfaceRangeIndex) >= surfaceRangeMin   %如果峰值在检查距离的范围
                            surfaceInspectNumMin = surfaceRangeIndex  ;  %确定检查范围内的起始峰值是第几个峰值
                            flagChange = 2;
                          end
                      else
                          if flagChange == 2
                              if JudgeBig(surfaceRangeIndex) > surfaceRangeMax   %如果峰值在检查距离的范围外，则检查结束
                                surfaceInspectNumMax = surfaceRangeIndex - 1; %确定检查范围内的终止峰值是第几个峰值
                                break;
                              end
                          end
                      end                    
                  end
                  %在检查距离内的较大值峰值的强度值
                  [surfaceInspectIntensity,surfaceInspectRelativePeak] = max(SeekBigIntensity(surfaceInspectNumMin:surfaceInspectNumMax));
%                   surfaceInspectIntensity   %较大峰值的强度值
                  surfaceInspectRelativePeak = surfaceInspectRelativePeak + surfaceInspectNumMin - 1; %较大峰值是第几个峰值
                  surfaceInspectRelativeLocation = JudgeBig(surfaceInspectRelativePeak);      %较大峰值的相对位置  
                  
                  %与上表面的强度差值为
                  surfaceDifference = surfaceIntensity - surfaceInspectIntensity;
                  %较大峰值与上表面峰值的差值  所占的上表面峰值的百分比
                  surfacePercentage =surfaceDifference / surfaceIntensity * 100;
%                   surfacePercentage =vpa(surfaceDifference / surfaceIntensity*100);
%                   disp(['surfacePercentage=',char(surfacePercentage),'%'])   %输出百分比

                  %判断是否需要替换位置 
                  surfaceInspectBiggerNum = 1;
                  if (surfacePercentage <= 20)  && (surfaceRelativeLocation - surfaceInspectRelativeLocation >= 0)   %如果小于20% 并且在峰值之前，则进行替换
                      surfaceRelativeLocation = surfaceInspectRelativeLocation ;
                      surfaceIntensity = surfaceInspectIntensity;
                      surfaceInspectBiggerNum = surfaceInspectBiggerNum + 1;  %有多少需要替换的A-Scan
                  end
% % % 下表面检查
% %                   判断在检查距离内有哪几个峰值？
%                   bottomRangeMin = round(bottomRelativeLocation - InspectRange);
%                   bottomRangeMax = round(bottomRelativeLocation + InspectRange);
% %                   [~,n] = size(SeekBigIndex); %重复了
%                   flagChange = 1;
%                   for bottomRangeIndex = 1:n
%                       if flagChange == 1
%                           if JudgeBig(bottomRangeIndex) >= bottomRangeMin   %如果峰值在检查距离的范围
%                             bottomInspectNumMin = bottomRangeIndex;    %确定检查范围内的起始峰值是第几个峰值
%                             flagChange = 2;
%                           end
%                       else
%                           if flagChange == 2
%                               if JudgeBig(bottomRangeIndex) >= bottomRangeMax   %如果峰值在检查距离的范围外，则检查结束
%                                 bottomInspectNumMax = bottomRangeIndex - 1; %确定检查范围内的终止峰值是第几个峰值
%                                 break;
%                               end
%                           end
%                       end                    
%                   end
%                   %在检查距离内的较大值峰值的强度值
%                   [bottomInspectIntensity,bottomInspectNumIndex] = max(SeekBigIntensity(bottomInspectNumMin:bottomInspectNumMax));
% %                   bottomInspectIntensity   %较大峰值的强度值
%                   bottomInspectRelativePeak = bottomInspectNumIndex + bottomInspectNumMin - 1; %较大峰值是第几个峰值
%                   bottomInspectRelativeLocation = JudgeBig(bottomInspectRelativePeak);      %较大峰值的相对位置  
%                   
%                   %与上表面的强度差值为
%                   bottomDifference = bottomIntensity - bottomInspectIntensity;
%                   %较大峰值与下表面峰值的差值  所占的下表面峰值的百分比
%                   bottomPercentage =bottomDifference / bottomIntensity*100;
% %                   bottomPercentage =vpa(bottomDifference / bottomIntensity*100);
% %                   disp(['bottomPercentage=',char(bottomPercentage),'%'])              % 输出百分比
% 
%                   %判断是否需要替换位置 
%                   bottomInspectBiggerNum = 1;
%                   if (bottomPercentage <= 20)  & (bottomInspectRelativeLocation >= bottomRelativeLocation)   %如果小于20% 并且在峰值之后，则进行替换 %下峰值优化算法
% %                   if (bottomPercentage <= 20)  & (bottomInspectRelativeLocation <= bottomRelativeLocation)   %如果小于20% 并且在峰值之前，则进行替换  %上峰值优化算法
%                       bottomRelativeLocation = bottomInspectRelativeLocation ;
%                       bottomIntensity = bottomInspectIntensity;
%                       bottomInspectBiggerNum = bottomInspectBiggerNum + 1;   %有多少需要替换的A-scan
%                   end

                  flag = 0;  %结束死循环
               else
                   SeekBigIntensity(max2NumIndex) = 0;  %如果在10pixel内有次大值出现，则为误差，赋0
               end             
           end
           
           surfaceRelativeLocation = surfaceRelativeLocation+rangeMin-1;     %上表面的峰值位置
           bottomRelativeLocation = bottomRelativeLocation+rangeMin-1;      %下表面的峰值位置
           depthNum = (abs(surfaceRelativeLocation-bottomRelativeLocation))*0.003493/1.38*10^3;   %厚度值的确定
            %这是A-Scan的循环结束，可以在这总结   
            %出现的问题是j--导致前面的为0
%             surfacePercentage(Index) = surfacePercentage;
%             bottomPercentage(Index) = bottomPercentage;

        else   %如果数据位置为0或者数据所在行围为0
%             problem = 0; %前期程序进行排查
            surfaceRelativeLocation = nan; 
            bottomRelativeLocation = nan;   
            surfaceIntensity = nan; 
            bottomIntensity = nan ; 
            depthNum = nan;
%             surfacePercentage = nan;
%             bottomPercentage = nan;
        end
            %对数据进行处理
            surfaceLocation(Index,IndexK) = surfaceRelativeLocation; %获取上表面峰值的位置
            bottomLocation(Index,IndexK) = bottomRelativeLocation;  %获取下表面峰值的位置
            surfaceIntensityA(Index,IndexK) = surfaceIntensity; %获取上表面的强度值
            bottomIntensityA(Index,IndexK) = bottomIntensity;   %获取下表面的强度值
            depth(Index,IndexK) = depthNum;  %获取皮肤厚度值
%             surfacePercentageA(Index,IndexK) = surfacePercentage; %为了获取上表面的阈值范围
%             bottomPercentageA(Index,IndexK) =  bottomPercentage; %为了获取上表面的阈值范围
    end
        %一个B-Scan内所有A-scan循环结束
        % 对A-Scan图像进行极大值极小值进行标注
%         hold on
%         scatter(JudgeBig,IntensityBig,'*');  %描绘极大值点
%         scatter(JudgeSmall,IntensitySmall,'^'); %描绘极小值点
%         plot([surfaceRelativeLocation,surfaceRelativeLocation],[maxPeak,max2Peak],'color','r','linewidth',2);   %起始端线条
%         plot([bottomRelativeLocation,bottomRelativeLocation],[maxPeak,max2Peak],'color','r','linewidth',2);   %终止端线条     
%         hold off   
        
        %上下表面的相对位置和强度
%         surfaceIntensity
%         bottomIntensity
%         surfaceRelativeLocation
%         bottomRelativeLocation
        
%           %进行阈值的频数统计
%           tabulate(surfacePercentageA(:))
%           fig2 = figure(2),clf;
%           histfit(surfacePercentageA);
%           i = 1;
%           print(fig2,'-dbitmap',strcat('D:\系统自带\桌面\ProgramChange\Percentage','i','SurfacePercentage.bmp'));
%           i = i +1;
   
    waitbar(k/line,wbSkin3dProcess,['计算皮肤厚度...' num2str(100*k/line) '%']);
end
close(wbSkin3dProcess);

%             SuanZi=fspecial('average',[3 3]); %% 定义一个滤波器 
%             averagesurfaceLocation=imfilter(surfaceLocation,SuanZi,'replicate');
%             averagebottomLocation=imfilter(bottomLocation,SuanZi,'replicate');
% 

% % 进行特殊点的消除优化
% [x,y]=size(depth)
% for i = 15:x-15
%     for j =15:y-15
%         if depth(i,j)-depth(i-1,j-1)>15
%             depth(i,j) = round((depth(i,j-10)+depth(i,j+10))/2);
%         end
%     end
% end

            %查看B-Scan的上下表面标识结果
%             [m,n] = size(surfaceLocation); 
%             fig1 = figure(1);hold on
%             for juli = 1:m
%                 surfaceEnd(juli+centerX-r-1) = averagesurfaceLocation(juli,line/2-centerY+r+1 );
%                 bottomEnd(juli+centerX-r-1) = averagebottomLocation(juli,line/2-centerY+r+1 );   
%             end
%             plot(surfaceEnd,'r');
%             plot(bottomEnd,'r');
%             hold off


%验证数据是否传输成功
% A = surfaceLocation(400,400)
% B = surfaceIntensityA(400,400)
% C = bottomLocation(400,400)
% D = bottomIntensityA(400,400)
end
%输出百分比的模板
%a=vpa(0.45);
%disp(['a=',char(a),'%'])