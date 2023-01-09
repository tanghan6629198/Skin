clc
clear
close all
%path = 'I:\20190910\Intensity3D\SR_4_DAY4\';
path = uigetdir('Select processed skin data path');      %选择数据文件夹
depth3D=load([path '\' 'depth3D.mat']);
surfaceLine=load([path  '\' 'surfaceLine.mat']);
bottomLine=load([path  '\' 'bottomLine.mat']);

b=bottomLine.bottomLine;
x=(-500:499)*0.009;
y=(-500:499)*0.009;
d=depth3D.depth_Line;

w2=fspecial('average',[9 9]); 
medDepth=imfilter(d,w2,'replicate');
% [x,y]=meshgrid(x,y);

%＃＃＃＃＃＃＃＃＃＃＃三维厚度分布＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃　

fig1 = figure(1);
%subplot(121);
mesh(x,y,medDepth);
set(gca,'FontSize',32);%,'FontSize',16
set(gca,'tickdir','in');
set(gca,'GridAlpha',0.8);
set(gca,'LineWidth',0.8);
axis([-3.5,3.5,-3.5,3.5,0,150])
colorbar;
colormap(jet);
caxis([0,150]);
xlabel('mm','FontName','Times New Roman','FontSize',40,'color','k');%x轴坐标
ylabel('mm','FontName','Times New Roman','FontSize',40,'color','k');%y轴坐标
zlabel('Thickness/\mum','FontName','Times New Roman','FontSize',40,'color','k');%z轴坐标

set(fig1,'position',[100 100 1000 800]);

print(fig1,'-dbitmap',strcat(path,'Depth.bmp'));
%＃＃＃＃＃＃＃＃＃＃＃三维厚度频数分布＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃　

fig2 = figure(2),

%subplot(122);
h=histogram((medDepth),30);
% hold on; 
%显示柱状图数值
% hBin=h.BinEdges(1:end-1)+h.BinWidth/2;
% text(hBin,h.Values+max(h.Values)/25,num2cell(h.Values));
%计算百分比
% Hpercent=round(h.Values/sum(h.Values)*100);
%加入百分号
% Hpercent2=num2cell(Hpercent);
% for i=1: length(Hpercent)
%     Hpercent2(i)={[num2str(Hpercent(i)),'%']};
% end
% text(hBin,h.Values+max(h.Values)/15,Hpercent2);%显示百分比

axis([20,130,0,40000])
set(gca,'FontSize',24);%,'FontSize',16
xlabel('Thickness/\mum','FontName','Times New Roman','FontSize',32,'color','k');%x轴坐标
ylabel('Frequrency','FontName','Times New Roman','FontSize',32,'color','k');%y轴坐标
% title('Thickness frequency distribution','FontName','Times New Roman','FontSize',28,'color','k');%z轴坐标

set(fig2,'position',[150 150 1000 800]);

print(fig2,'-dbitmap',strcat(path,'DepthHistogram.bmp'));
%% Ra

s=surfaceLine.surfaceLine*0.003493/1.38*10^3;
[X,Y] = meshgrid(x,y);
[fitresult, gof] = createFit(X, Y, s);
fitSurface = fitresult(X,Y);
flatData=(s-fitSurface);

w2=fspecial('average',[5 5]); %% 定义一个滤波器 
filtSurface00=imfilter(flatData,w2,'replicate');
filtSurface0 = abs(filtSurface00-max(filtSurface00(:)));
% medSurface=imrotate(medSurface, -90);
 
Ra = 0 ; sum = 0;count = 0;thickness=0;st=0;
for i=1:1000
    for j=1:1000
        if ~isnan(medDepth(i,j))
            sum =sum + filtSurface0(i,j);
            thickness = thickness+medDepth(i,j);
            count = count+1;
        else
            filtSurface(i,j)=nan;  %绘图需要
        end
    end   
end
thickness = thickness/(count); %平均厚度

num = 0;
for i=1:1000
    for j=1:1000
        if ~isnan(medDepth(i,j))%(~isnan(filtSurface(i,j))&(filtSurface(i,j)~=0))%
            st = st+(d(i,j)-thickness).^2;
           if medDepth(i,j)>=thickness
               num = num +1;
           end
        end
    end   
end
rate = num/count;%频数分布中大于均值厚度的占比
st = sqrt(1/(count)*st);%厚度标准差

%＃＃＃＃＃＃＃＃＃＃＃皮肤表面形态＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃　

fig3 = figure(3);
[c,h]=contourf(x,y,filtSurface0,5);%等高线
set(gca,'FontSize',32);%,'FontSize',16
set(gca,'LineWidth',1.5);
% set(h,'LevelList',[0 0.15 0.2 0.3])%设定等高线的值  ,'ShowText','on'
colorbar;
colormap(jet);
axis([-3.5,3.5,-3.5,3.5])
xlabel('mm','FontName','Times New Roman','FontSize',42,'color','k');%x轴坐标
ylabel('mm','FontName','Times New Roman','FontSize',42,'color','k');%y轴坐标


set(fig3,'position',[200 200 1000 800]);

print(fig3,'-dbitmap',strcat(path,'skinSurface.bmp'));
%＃＃＃＃＃＃＃＃＃＃＃皮肤表面粗糙度＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃＃　

meanSueface=sum/(count);
 Rq=0;
for i=1:1000
    for j=1:1000
        if ~isnan(filtSurface0(i,j))
             Ra = Ra + abs(filtSurface0(i,j)-meanSueface);
%              Rq = Rq + (filtSurface0(i,j)-meanSueface);
              
        end
    end   
end
Ra = 1/(count)*Ra;%粗糙度
% Rq = sqrt(1/(count)*Rq);
% Rmax = max(filtSurface0(:))-min(filtSurface0(:));



