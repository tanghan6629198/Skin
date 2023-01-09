 function [filtSurface,SurfaceTrue] = saveSkinDataImageV2_5(surfaceLocation,~,depth,r,savePath,name)


%*******************************************************

%功能：用depth求Ra
%码农：汤瀚
%时间：2022.1.6
%Matlab版本：2020a

%*******************************************************

%   UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%%  验证数据传输



%% 进行数据处理
depth2 = depth;
%% 平均面
mean_depth = 0;num = 0;
centerX = r+1;
centerY = r+1;
for i = 1:2*r+1
    for j = 1:2*r+1
        if~isnan(depth(i,j))
          	mean_depth = mean_depth + depth(i,j);
            num = num + 1;
        end
    end
end
mean_depth = mean_depth/(num);
%整体平均面
for i = 1:2*r+1
    for j = 1:2*r+1
        if~isnan(depth(i,j))
          	depth2(i,j) = mean_depth;
        end
    end
end
filtSurface = depth - depth2;
filtSurface = medfilt2(filtSurface,[3,3]); % 中值滤波去椒盐噪声

x=(-r:r)*0.01; %x轴的像素大小以及分布
y=(-r:r)*0.01; %y轴的像素大小以及分布
%% 拟合面
% [X,Y] = meshgrid(x,y);
% [fitresult, ~] = createFit(X, Y, depth);  %二元三次拟合以减去曲率
% fitSurface = fitresult(X,Y);
% % flatData=(s-fitSurface);
% % filtSurface = (s-fitSurface);  %未拟合的上表面-拟合后上表面   0 有正有负 去曲率
% filtSurface = (depth-fitSurface);  %拟合后上表面-未拟合的上表面   0 有正有负 去曲率
% filtSurface = medfilt2(filtSurface,[3,3]); % 中值滤波去椒盐噪声


%进行平滑处理--用以成图美观
w2=fspecial('average',[9 9]); 
averageDepth=imfilter(depth,w2,'replicate');  %平滑过后的厚度值

%%对样本厚度进行计算
Ra = 0 ; sum = 0;count = 0;thickness=0;st=0;
for i=1:2*r+1
    for j=1:2*r+1
        if ~isnan(averageDepth(i,j))
            sum =sum + filtSurface(i,j);  
            thickness = thickness+averageDepth(i,j);
            count = count+1;
        else
            filtSurface(i,j)=nan;  %绘图需要
        end
    end   
end
thickness = thickness/(count)  %平均厚度
% averageDepth存储的厚度信息
% %%
% meanDepth = mean(averageDepth);


%%
% 分布范围间距（RD）、大于均值占比(MP)以及频数峰值处厚度（PT）

num = 0;
for i=1:2*r+1
    for j=1:2*r+1
        if ~isnan(averageDepth(i,j))%(~isnan(filtSurface(i,j))&(filtSurface(i,j)~=0))%
            st = st+(averageDepth(i,j)-thickness).^2;
           if averageDepth(i,j)>=thickness
               num = num +1;
           end
        end
    end   
end
rateBig = num/count   
st = sqrt(1/(count)*st)  %厚度标准差
 
DepthMin = min(averageDepth);
DepthMin = min(DepthMin);
DepthMax = max(averageDepth);
DepthMax = max(DepthMax);
Depthm = DepthMax - DepthMin

%%粗糙度

filtSurface2 = filtSurface(91:481,91:481);

centerX = 196;
centerY = 196;
rr= 195;
rr2 = 195*195;
for j = 1:391
    for k = 1:391
        if((j-centerX)^2+(k-centerY)^2)>rr2
          	filtSurface2(k,j)=nan;
        end
    end
end
%%取中间200*200像素
for i = 1:401
        s401(:,i) = filtSurface(86:486,85+i);
end
sRa = zeros(390,390);
for i = 1:391
    for j = 1:391
        si = s401(i:i+10,j:j+10);
        si = abs(si);
        Rai = mean(mean(si));
        sRa(i,j) = Rai;
    end
end
for j = 1:391
    for k = 1:391
        if((j-centerX)^2+(k-centerY)^2)>rr2
          	sRa(k,j)=nan;
        end
    end
end
ssRa = 0;count=0;
for i=1:2*rr+1
    for j=1:2*rr+1
        if ~isnan(sRa(i,j))           
             ssRa = ssRa + sRa(i,j);
             count =  count+1;
        end
    end   
end
ssRa = ssRa/count
for i=1:2*rr+1
    for j=1:2*rr+1
        if ~isnan(sRa(i,j))
             SurfaceTrue(i,j) = abs(sRa(i,j)-ssRa);
        end
    end   
end
[a,b]=size(SurfaceTrue);
num = 0;Rsk=0;Rku=0;q=0;sk=0;ku=0;count = 0;
for i=1:a
    for j=1:b
        if SurfaceTrue(i,j) ~= 0
              q = q + SurfaceTrue(i,j).^2;  %来计算Rq，距离的平方/A
              sk = sk + SurfaceTrue(i,j).^3; %来计算Rsk，距离的三次方/A
              ku = ku + SurfaceTrue(i,j).^4; %来计算Rku，距离的四次方/A
              num=num+1;
        end
    end   
end
Rq = sqrt(1/(num)*q)  % 均方根高度
Rsk = (sk/num)/(Rq^3) % 偏斜度
Rku = (ku/num)/(Rq^4) % 尖锐度

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 统计区域等分

%%%%%%%%%%%%%%% 1.垂直
centerX = r+1;
centerY = r+1;
Bscan_1 = zeros(1,500,'single')*NaN;
Bscan_1(1,250) = filtSurface(centerX,centerY);
for i = 1:249
    if ~isnan(filtSurface(centerX-1,centerY))
    Bscan_1(1,250 - i) = filtSurface(centerX-i,centerY);
    else
        break
    end
end
for i = 1:250
    if ~isnan(filtSurface(centerX+1,centerY))
    Bscan_1(1,250 + i) = filtSurface(centerX+i,centerY);
    else
        break
    end
end
%%%%%%%%%%%%%%% 2.
centerX = r+1;
centerY = r+1;
Bscan_2 = zeros(1,500,'single')*NaN;
Bscan_2(1,250) = filtSurface(centerX,centerY);
for i = 1:249
    if ~isnan(filtSurface(centerX-2*i,centerY+i))
    Bscan_2(1,250 - i) = filtSurface(centerX-2*i,centerY+i);
    else
        break
    end
end
for i = 1:250
    if ~isnan(filtSurface(centerX+2*i,centerY-i))
    Bscan_2(1,250 + i) = filtSurface(centerX+2*i,centerY-i);
    else
        break
    end
end

%%%%%%%%%%%%%%% 3.
centerX = r+1;
centerY = r+1;
Bscan_3 = zeros(1,500,'single')*NaN;
Bscan_3(1,250) = filtSurface(centerX,centerY);
for i = 1:249
    if ~isnan(filtSurface(centerX- i,centerY+ i))
    Bscan_3(1,250 - i) = filtSurface(centerX - i,centerY + i);
    else
        break
    end
end
for i = 1:250
    if ~isnan(filtSurface(centerX+ i,centerY- i))
    Bscan_3(1,250 + i) = filtSurface(centerX + i,centerY - i);
    else
        break
    end
end

%%%%%%%%%%%%%%% 4.
centerX = r+1;
centerY = r+1;
Bscan_4 = zeros(1,500,'single')*NaN;
Bscan_4(1,250) = filtSurface(centerX,centerY);
for i = 1:249
    if ~isnan(filtSurface(centerX- i,centerY+ 2*i))
    Bscan_4(1,250 - i) = filtSurface(centerX - i,centerY + 2*i);
    else
        break
    end
end
for i = 1:250
    if ~isnan(filtSurface(centerX+ i,centerY- 2*i))
    Bscan_4(1,250 + i) = filtSurface(centerX + i,centerY - 2*i);
    else
        break
    end
end

%%%%%%%%%%%%%%% 5.水平
centerX = r+1;
centerY = r+1;
Bscan_5 = zeros(1,500,'single')*NaN;
Bscan_5(1,250) = filtSurface(centerX,centerY);
for i = 1:249
    if ~isnan(filtSurface(centerX,centerY+ i))
    Bscan_5(1,250 - i) = filtSurface(centerX ,centerY + i);
    else
        break
    end
end
for i = 1:250
    if ~isnan(filtSurface(centerX,centerY - i))
    Bscan_5(1,250 + i) = filtSurface(centerX ,centerY - i);
    else
        break
    end
end

%%%%%%%%%%%%%%% 6.
centerX = r+1;
centerY = r+1;
Bscan_6 = zeros(1,500,'single')*NaN;
Bscan_6(1,250) = filtSurface(centerX,centerY);
for i = 1:249
    if ~isnan(filtSurface(centerX + i,centerY+ 2*i))
    Bscan_6(1,250 - i) = filtSurface(centerX + i,centerY + 2*i);
    else
        break
    end
end
for i = 1:250
    if ~isnan(filtSurface(centerX - i,centerY - 2*i))
    Bscan_6(1,250 + i) = filtSurface(centerX - i,centerY - 2*i);
    else
        break
    end
end

%%%%%%%%%%%%%%% 7.
centerX = r+1;
centerY = r+1;
Bscan_7 = zeros(1,500,'single')*NaN;
Bscan_7(1,250) = filtSurface(centerX,centerY);
for i = 1:249
    if ~isnan(filtSurface(centerX + i,centerY+ i))
    Bscan_7(1,250 - i) = filtSurface(centerX + i,centerY + i);
    else
        break
    end
end
for i = 1:250
    if ~isnan(filtSurface(centerX - i,centerY - i))
    Bscan_7(1,250 + i) = filtSurface(centerX - i,centerY - i);
    else
        break
    end
end

%%%%%%%%%%%%%%% 8.
centerX = r+1;
centerY = r+1;
Bscan_8 = zeros(1,500,'single')*NaN;
Bscan_8(1,250) = filtSurface(centerX,centerY);
for i = 1:249
    if ~isnan(filtSurface(centerX + 2*i,centerY+ i))
    Bscan_8(1,250 - i) = filtSurface(centerX + 2*i,centerY + i);
    else
        break
    end
end
for i = 1:250
    if ~isnan(filtSurface(centerX - 2*i,centerY - i))
    Bscan_8(1,250 + i) = filtSurface(centerX - 2*i,centerY - i);
    else
        break
    end
end


%%%%%%%%%求8个RZ
%%%%% 1 /5
RmMin_1(1,1) = min(Bscan_1(1:100));
RmMin_1(1,2) = min(Bscan_1(101:200));
RmMin_1(1,3) = min(Bscan_1(201:300));
RmMin_1(1,4) = min(Bscan_1(301:400));
RmMin_1(1,5) = min(Bscan_1(401:500));
RmMax_1(1,1) = max(Bscan_1(1:100));
RmMax_1(1,2) = max(Bscan_1(101:200));
RmMax_1(1,3) = max(Bscan_1(201:300));
RmMax_1(1,4) = max(Bscan_1(301:400));
RmMax_1(1,5) = max(Bscan_1(401:500));
Rz_1(1,1) = RmMax_1(1,1) - RmMin_1(1,1);
Rz_1(1,2) = RmMax_1(1,2) - RmMin_1(1,2);
Rz_1(1,3) = RmMax_1(1,3) - RmMin_1(1,3);
Rz_1(1,4) = RmMax_1(1,4) - RmMin_1(1,4);
Rz_1(1,5) = RmMax_1(1,5) - RmMin_1(1,5);
Rz1 = mean(Rz_1);

RmMin_5(1,1) = min(Bscan_5(1:100));
RmMin_5(1,2) = min(Bscan_5(101:200));
RmMin_5(1,3) = min(Bscan_5(201:300));
RmMin_5(1,4) = min(Bscan_5(301:400));
RmMin_5(1,5) = min(Bscan_5(401:500));
RmMax_5(1,1) = max(Bscan_5(1:100));
RmMax_5(1,2) = max(Bscan_5(101:200));
RmMax_5(1,3) = max(Bscan_5(201:300));
RmMax_5(1,4) = max(Bscan_5(301:400));
RmMax_5(1,5) = max(Bscan_5(401:500));
Rz_5(1,1) = RmMax_5(1,1) - RmMin_5(1,1);
Rz_5(1,2) = RmMax_5(1,2) - RmMin_5(1,2);
Rz_5(1,3) = RmMax_5(1,3) - RmMin_5(1,3);
Rz_5(1,4) = RmMax_5(1,4) - RmMin_5(1,4);
Rz_5(1,5) = RmMax_5(1,5) - RmMin_5(1,5);
Rz5 = mean(Rz_5);

%%%%% 2 /8
RmMin_2(1,1) = min(Bscan_2(123:173));
RmMin_2(1,2) = min(Bscan_2(174:224));
RmMin_2(1,3) = min(Bscan_2(225:275));
RmMin_2(1,4) = min(Bscan_2(276:326));
RmMin_2(1,5) = min(Bscan_2(327:377));
RmMax_2(1,1) = max(Bscan_2(123:173));
RmMax_2(1,2) = max(Bscan_2(174:224));
RmMax_2(1,3) = max(Bscan_2(225:275));
RmMax_2(1,4) = max(Bscan_2(276:326));
RmMax_2(1,5) = max(Bscan_2(327:377));
Rz_2(1,1) = RmMax_2(1,1) - RmMin_2(1,1);
Rz_2(1,2) = RmMax_2(1,2) - RmMin_2(1,2);
Rz_2(1,3) = RmMax_2(1,3) - RmMin_2(1,3);
Rz_2(1,4) = RmMax_2(1,4) - RmMin_2(1,4);
Rz_2(1,5) = RmMax_2(1,5) - RmMin_2(1,5);
Rz2 = mean(Rz_2);

RmMin_8(1,1) = min(Bscan_8(123:173));
RmMin_8(1,2) = min(Bscan_8(174:224));
RmMin_8(1,3) = min(Bscan_8(225:275));
RmMin_8(1,4) = min(Bscan_8(276:326));
RmMin_8(1,5) = min(Bscan_8(327:377));
RmMax_8(1,1) = max(Bscan_8(123:173));
RmMax_8(1,2) = max(Bscan_8(174:224));
RmMax_8(1,3) = max(Bscan_8(225:275));
RmMax_8(1,4) = max(Bscan_8(276:326));
RmMax_8(1,5) = max(Bscan_8(327:377));
Rz_8(1,1) = RmMax_8(1,1) - RmMin_8(1,1);
Rz_8(1,2) = RmMax_8(1,2) - RmMin_8(1,2);
Rz_8(1,3) = RmMax_8(1,3) - RmMin_8(1,3);
Rz_8(1,4) = RmMax_8(1,4) - RmMin_8(1,4);
Rz_8(1,5) = RmMax_8(1,5) - RmMin_8(1,5);
Rz8 = mean(Rz_8);

%%%%% 3/7
RmMin_3(1,1) = min(Bscan_3(51:130));
RmMin_3(1,2) = min(Bscan_3(131:210));
RmMin_3(1,3) = min(Bscan_3(211:290));
RmMin_3(1,4) = min(Bscan_3(291:370));
RmMin_3(1,5) = min(Bscan_3(371:450));
RmMax_3(1,1) = max(Bscan_3(51:130));
RmMax_3(1,2) = max(Bscan_3(131:210));
RmMax_3(1,3) = max(Bscan_3(211:290));
RmMax_3(1,4) = max(Bscan_3(291:370));
RmMax_3(1,5) = max(Bscan_3(371:450));
Rz_3(1,1) = RmMax_3(1,1) - RmMin_3(1,1);
Rz_3(1,2) = RmMax_3(1,2) - RmMin_3(1,2);
Rz_3(1,3) = RmMax_3(1,3) - RmMin_3(1,3);
Rz_3(1,4) = RmMax_3(1,4) - RmMin_3(1,4);
Rz_3(1,5) = RmMax_3(1,5) - RmMin_3(1,5);
Rz3 = mean(Rz_3);

RmMin_7(1,1) = min(Bscan_7(51:130));
RmMin_7(1,2) = min(Bscan_7(131:210));
RmMin_7(1,3) = min(Bscan_7(211:290));
RmMin_7(1,4) = min(Bscan_7(291:370));
RmMin_7(1,5) = min(Bscan_7(371:450));
RmMax_7(1,1) = max(Bscan_7(51:130));
RmMax_7(1,2) = max(Bscan_7(131:210));
RmMax_7(1,3) = max(Bscan_7(211:290));
RmMax_7(1,4) = max(Bscan_7(291:370));
RmMax_7(1,5) = max(Bscan_7(371:450));
Rz_7(1,1) = RmMax_7(1,1) - RmMin_7(1,1);
Rz_7(1,2) = RmMax_7(1,2) - RmMin_7(1,2);
Rz_7(1,3) = RmMax_7(1,3) - RmMin_7(1,3);
Rz_7(1,4) = RmMax_7(1,4) - RmMin_7(1,4);
Rz_7(1,5) = RmMax_7(1,5) - RmMin_7(1,5);
Rz7 = mean(Rz_7);

%%%%% 4/6
RmMin_4(1,1) = min(Bscan_4(123:173));
RmMin_4(1,2) = min(Bscan_4(174:224));
RmMin_4(1,3) = min(Bscan_4(225:275));
RmMin_4(1,4) = min(Bscan_4(276:326));
RmMin_4(1,5) = min(Bscan_4(327:377));
RmMax_4(1,1) = max(Bscan_4(123:173));
RmMax_4(1,2) = max(Bscan_4(174:224));
RmMax_4(1,3) = max(Bscan_4(225:275));
RmMax_4(1,4) = max(Bscan_4(276:326));
RmMax_4(1,5) = max(Bscan_4(327:377));
Rz_4(1,1) = RmMax_4(1,1) - RmMin_4(1,1);
Rz_4(1,2) = RmMax_4(1,2) - RmMin_4(1,2);
Rz_4(1,3) = RmMax_4(1,3) - RmMin_4(1,3);
Rz_4(1,4) = RmMax_4(1,4) - RmMin_4(1,4);
Rz_4(1,5) = RmMax_4(1,5) - RmMin_4(1,5);
Rz4 = mean(Rz_4);

RmMin_6(1,1) = min(Bscan_6(123:173));
RmMin_6(1,2) = min(Bscan_6(174:224));
RmMin_6(1,3) = min(Bscan_6(225:275));
RmMin_6(1,4) = min(Bscan_6(276:326));
RmMin_6(1,5) = min(Bscan_6(327:377));
RmMax_6(1,1) = max(Bscan_6(123:173));
RmMax_6(1,2) = max(Bscan_6(174:224));
RmMax_6(1,3) = max(Bscan_6(225:275));
RmMax_6(1,4) = max(Bscan_6(276:326));
RmMax_6(1,5) = max(Bscan_6(327:377));
Rz_6(1,1) = RmMax_6(1,1) - RmMin_6(1,1);
Rz_6(1,2) = RmMax_6(1,2) - RmMin_6(1,2);
Rz_6(1,3) = RmMax_6(1,3) - RmMin_6(1,3);
Rz_6(1,4) = RmMax_6(1,4) - RmMin_6(1,4);
Rz_6(1,5) = RmMax_6(1,5) - RmMin_6(1,5);
Rz6 = mean(Rz_6);

Rz = (Rz1+Rz2+Rz3+Rz4+Rz5+Rz6+Rz7+Rz8)/8 %%   Rz  8条平均

Bscan = [Bscan_1;Bscan_2;Bscan_3;Bscan_4;Bscan_5;Bscan_6;Bscan_7;Bscan_8];
[maxx,maxy]  = find(Bscan == max(max(Bscan)));
[minx,miny]  = find(Bscan == min(min(Bscan)));
Rm = max(max(Bscan)) - min(min(Bscan))


% % ＃＃＃＃＃＃＃＃＃＃＃2D厚度分布＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃　
% % 
% w2=fspecial('average',[3 3]); %% 定义一个滤波器 
% averageDepth=imfilter(averageDepth,w2,'replicate');

fig2 = figure(2);
mesh(x,y,averageDepth);
set(gca,'FontSize',32);%,'FontSize',16
set(gca,'tickdir','in');
set(gca,'GridAlpha',0.8);
set(gca,'LineWidth',0.8);
% zlim([30,80]);
axis([-2.5,2.5,-2.5,2.5,0,150])
% axis([-2.5,2.5,-2.5,2.5,50,200])  %%样本比较厚
colorbar;
colormap(jet);
caxis([0,150]);
% caxis([50,200]);%%样本比较厚
xlabel('mm','FontName','Times New Roman','FontSize',40,'color','k');%x轴坐标
ylabel('mm','FontName','Times New Roman','FontSize',40,'color','k');%y轴坐标
zlabel('Thickness\um','FontName','Times New Roman','FontSize',40,'color','k');%z轴坐标
title('SkinDepth','FontName','Times New Roman','FontSize',40,'color','k'); %标题

set(fig2,'position',[100 100 1000 800]);

print(fig2,'-dbitmap',strcat(savePath,'\',name,'_Depth.bmp'));
% print(fig2,'-dbitmap',strcat(savePath,'\',name,'_Depth.fig'));
savefig(fig2,strcat(savePath,'\',name,'Depth.fig'),'compact');
%＃＃＃＃＃＃＃＃＃＃＃2D皮肤表面形态＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃　
fig107 = figure(107);
x1=(-rr:rr)*0.01; %x轴的像素大小以及分布
y1=(-rr:rr)*0.01; %y轴的像素大小以及分布
h=imagesc([x1 x1],[y1 y1],filtSurface2);
set(h,'alphadata',~isnan(filtSurface2));
set(gca,'FontSize',32);%,'FontSize',16
set(gca,'LineWidth',1.5);
colorbar;
colormap(jet);
caxis([-5,5]);
set(gca,'xtick',-2:1:2);   
set(gca,'ytick',-2:1:2);  
axis([-2,2,-2,2]);
xlabel('mm','FontName','Times New Roman','FontSize',32,'color','k');%x轴坐标
ylabel('mm','FontName','Times New Roman','FontSize',32,'color','k');%y轴坐标
title('texture','FontName','Times New Roman','FontSize',32,'color','k'); %标题
set(fig107,'position',[350 350 1000 800]);
print(fig107,'-dbitmap',strcat(savePath,'\',name,'texture.bmp'));
savefig(fig107,strcat(savePath,'\',name,'texture.fig'),'compact');

fig108 = figure(108);
x2=(-rr:rr)*0.01; %x轴的像素大小以及分布
y2=(-rr:rr)*0.01; %y轴的像素大小以及分布
h=imagesc([x2 x2],[y2 y2],sRa);
set(h,'alphadata',~isnan(sRa));
set(gca,'FontSize',32);%,'FontSize',16
set(gca,'LineWidth',1.5);
colorbar;
colormap(jet);
caxis([0,5]);
set(gca,'xtick',-2:1:2);   
set(gca,'ytick',-2:1:2);  
axis([-2,2,-2,2]);
xlabel('mm','FontName','Times New Roman','FontSize',32,'color','k');%x轴坐标
ylabel('mm','FontName','Times New Roman','FontSize',32,'color','k');%y轴坐标
title('Ra','FontName','Times New Roman','FontSize',32,'color','k'); %标题
set(fig108,'position',[350 350 1000 800]);
print(fig108,'-dbitmap',strcat(savePath,'\',name,'Ra.bmp'));
savefig(fig108,strcat(savePath,'\',name,'Ra.fig'),'compact');

%＃＃＃＃＃＃＃＃＃＃＃粗糙度的频数统计＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃　

fig105 = figure(105);clf;
h = histogram(sRa,'Normalization','probability','binwidth',0.1);
set(gca,'FontSize',24);%,'FontSize',16
% set(gca,'LineWidth',1.5);
xlabel('um','FontName','Times New Roman','FontSize',32,'color','k');%x轴坐标
ylabel('Frequency','FontName','Times New Roman','FontSize',32,'color','k');%y轴坐标
title('Ra frequency distribution','FontName','Times New Roman','FontSize',32,'color','k'); %标题
xlim([0,10]);
ylim([0,0.1]);
set(fig105,'position',[350 350 1000 800]);
print(fig105,'-dbitmap',strcat(savePath,'\',name,'Ra frequency distribution.bmp'));
savefig(fig105,strcat(savePath,'\',name,'Ra frequency distribution.fig'),'compact');

fig106 = figure(106);
[~,~]=contourf(x,y,filtSurface);%等高线
set(gca,'FontSize',24);%,'FontSize',16 set(gca,'LineWidth',1.5);
colorbar; colormap(jet); caxis([-20,20]); axis([-2.85,2.85,-2.85,2.85])
xlabel('mm','FontName','Times New Roman','FontSize',32,'color','k');%x轴坐标
ylabel('mm','FontName','Times New Roman','FontSize',32,'color','k');%y轴坐标
title('Surface distribution','FontName','Times NewRoman','FontSize',32,'color','k'); %标题 set(fig6,'position',[350 350 1000 800]); % text(2.5,3,num2str(Ra),'FontSize',20,'color','r' );
set(fig106,'position',[350 350 1000 800]);
print(fig106,'-dbitmap',strcat(savePath,'\',name,'Surface.bmp'));
savefig(fig106,strcat(savePath,'\',name,'Surface.fig'),'compact');
% h=imagesc(filtSurface);
% set(h,'alphadata',~isnan(filtSurface));
% colorbar;
% colormap(jet);
% caxis([-5,5]);
close all

end

