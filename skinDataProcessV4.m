    function [surfaceLocation,bottomLocation,depth,r] = skinDataProcessV4( filePath )
% function [depth_Line,surfaceLine,bottomLine] = skinDataProcessV2( filePath )
%*******************************************************
%��ͬƤ�����������޸ģ� A_scan  B_scan  line   Depth
%���ܣ������㷨�Ż���ıȽ�
%��ɶȣ��������㷨
%��ũ��Zhao Ruihang
%ʱ�䣺2020.05.07

%*******************************************************

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% ������������
% ���µĲ���
% A_scan=1024;  % ��������
% B_scan=1024;   % the times of A-scan in per B-scan   %A-Scan������
% line=1024;     % the times of B-scan    %B-Scan������S


A_scan=1024;  % ��������
B_scan=1000;   % the times of A-scan in per B-scan   %A-Scan������
line=1000;     % the times of B-scan    %B-Scan������
r_range = 285;
Depth = 18;      %%Ƥ�� �����ֵ Day��ͬ  Ҫ�޸�
%% ��ȡ����

shft = B_scan*A_scan;  %һ��bscan�ܵĲ�������
fid=fopen(filePath,'r','n');
wb = waitbar(0,'��ȡ3D������...');
disp('��ȡ3D������...');
data3D = zeros(A_scan,B_scan,line,'single');
A = zeros(A_scan,B_scan,'single')*NaN;
for i=1:line   
    A = fread(fid,[A_scan,B_scan],'float32','n'); 
    if isempty(A)==1    %��ʾ�߼���������ǡ���Ҳ����ȡ��
        break;  %��A�ǿվ���ʱ������
    end
    fseek(fid, 4*i*shft, 'bof');  %ָ�������j*shft��Ԫ��
    data3D(:,:,i)=A;
    waitbar(i/line,wb,['��ȡ3D������...' num2str(100*i/line) '%']);
end
clear A;
close(wb);
% for i = 1:line/2
%     for j = 1:B_scan/2
%         data3Dnew(:,j,i) = data3D(:,2*j-1,2*i-1);
%     end
% end
%% �ü���Ч����
wb = waitbar(0,'�ü���Ч����...');
disp('�ü���Ч����...');

% ȷ����Χ
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

range = rangeMin:rangeMax;  %Ƥ���źŷ�Χ

% Thrmo�ı�Ե���
% section90 = zeros(A_scan,B_scan,'single')*NaN;
% section180 = zeros(A_scan,B_scan,'single')*NaN;
section90 = data3D(:,:,line/2);
section180 = data3D(:,B_scan/2,:);
%% ��Եʶ��


centerX = B_scan/2;
centerY = line/2;

% r = 36; %�뾶0.25��ֱ�� 0.5
%  r = 71;  % �뾶��0.5mm��ֱ��Ϊ 1
% r = 398; %�뾶 2.5��ֱ�� 5
% r = 278;
r=r_range;     %��θ������������з�Χ��ȷ��  %2.5mm
% r = 200;

r2 = r^2;
waitbar(1,wb,['�ü���Ч����...' num2str(20)  '%']);
nanLine = zeros(A_scan,1)*NaN('single');

for j = 1:B_scan
    for k = 1:line
        if((j-centerX)^2+(k-centerY)^2)>r2
          	data3D(:,j,k)=nanLine;
        end
    end
    waitbar(j/B_scan,wb,['�ü���Ч����...' num2str(20 + j/B_scan*80) '%']);
end

close(wb);

p1 = [centerX-r,centerX+r];   %x�������
p2 = [rangeMax,rangeMin];   %y�������

% %% ԭʼ��ͼ
% figure(1),clf;
% imagesc(section2D);
% hold on
% p1 = [range_begin90+20,range_end90-20];   %x�������
% p2 = [rangeMax,rangeMin];   %y�������
% plot([p1(1),p1(1)],[p2(1),p2(2)],'color','r','linewidth',2);   %��ʼ������
% plot([p1(2),p1(2)],[p2(1),p2(2)],'color','r','linewidth',2);   %��ֹ������
% hold off

% xlabel('mm','FontName','Times New Roman','FontSize',40,'color','k');%x������
% ylabel('mm','FontName','Times New Roman','FontSize',40,'color','k');%y������
% zlabel('Thickness/\mum','FontName','Times New Roman','FontSize',40,'color','k');%z������

%% �ָ�Ƥ������

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
wbSkin3dProcess = waitbar(100,'����Ƥ�����...');
zeroline = zeros(rangeMax +1 - rangeMin,1);
IndexK = 0;
for k = centerY - r  : centerY  + r     %ÿ��B-Scan���ν��д���     %ѡ��B-scan��ѡ��������B-Scan
%     k;
    IndexK = IndexK + 1;
    Intensity = data3D(range,:,k); %��ȡ��ǿ���ź�
    
    %���ҵ�ÿ��B-Scan��λ�ú���ֵ��С
    %�ж��ǲ��Ǽ�ֵ�㣺���׵�����>0��ô���Ǽ�Сֵ����֮Ϊ����ֵ
    %�ж���������ĵ���ֵ������+��- Ϊ����ֵ�㡣��-��+ Ϊ��Сֵ�㡣
    i_big = 1;
    i_small = 1;
    Index = 0;
    for j = centerX - r : centerX  + r     %ÿ��A-Scan���δ���   %ѡ��B-Scan������A-scan
%         j;
        Index = Index + 1;
        Data = Intensity(:,j);
        DataDiff = diff(Intensity(:,j));
        
%         %A-Scan����
%         figure(k),clf;
%         plot(Data);
        
        if (~isnan(Intensity(1,j))&&all(Intensity(:,j)~=zeroline))
%             problem = 1;  %���ڳ�������Ų�ı�־λ
%ÿ��A-scan������0����ֹ����            
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
%�жϼ���ֵ�ͼ�Сֵ��λ�ú�ǿ��ֵ
            for rangeIndex = 1:rangeMax-rangeMin-1
                if DataDiff(rangeIndex)>0 && DataDiff(rangeIndex+1) <0   %����ֵ��
%                     IndexJudge = 1; %ǰ�ڳ�����д����Ų��־λ
                    JudgeBig(i_big) = rangeIndex + 1;          %����ֵ��λ��
                    IntensityBig(i_big) = Data(rangeIndex + 1);   %����ֵ��ǿ��ֵ
                    i_big = i_big + 1;
                else
                    if DataDiff(rangeIndex) <0 && DataDiff(rangeIndex+1) >0 %��Сֵ��
%                         IndexJudge = 2; %ǰ�ڳ�����д����Ų�ı�־λ
                        JudgeSmall(i_small) = rangeIndex + 1;  %��Сֵ��λ��
                        IntensitySmall(i_small) = Data(rangeIndex + 1);  %��Сֵ��ǿ��ֵ
                        i_small = i_small +1;
                    end
                end
            end
%��A-Scan���б��涨λ
%            [m,n] = size(JudgeBig);
%            SeekBigIndex = JudgeBig;         % ���¸�ֵ����ֵ��λ�ã���ֹԭʼ���ݶ�ʧ--�ô�����
           SeekBigIntensity = IntensityBig; % ���¸�ֵ����ֵ����ǿ��ֵ����ֹԭʼ���ݶ�ʧ
           [~,maxNumIndex] = max(IntensityBig); %������ֵ��ǿ��ֵ �Լ�  ���ֵ�ǵڼ�����ֵ
           maxIndex = JudgeBig(maxNumIndex);          %������ֵ���������λ��  ��ʵ������λ��Ӧ�ü���rangeMin
           SeekBigIntensity(maxNumIndex) = 0; %���ֵ��ǿ��ֵ��Ϊ0���Բ�����һ����ֵ

           
           flag = 1 ;  %������ѭ�������жϷ�ֵ��׼ȷ��
           while(flag)                
               %������һ����ֵ
               [~,max2NumIndex] = max(SeekBigIntensity); %��ôδ�ֵ��ǿ��ֵ  �Լ�  �δ�ֵ�ǵڼ�����ֵ
               max2Index = JudgeBig(max2NumIndex);          %��ôδ�ֵ���������λ��
               SeekBigIntensity(max2NumIndex) = 0;              %�δ�ֵ��ǿ��ֵ��Ϊ0
               
               
               Day2 = Depth;
               if abs(maxIndex - max2Index) <= Depth
                   Day2 = Depth - 5;   %Ƥ����Ȳ�����ʱʹ��
               end
               if abs(maxIndex - max2Index) > Day2  %10pixel = 25um.������ֵ����С��25um--Day01�ĺ��--9DayΪ15--13DayΪ20
                   
                  %����ϱ�����±�������λ�úͶ�Ӧ��ǿ��ֵ 
                  surfaceRelativePeak = min(maxNumIndex,max2NumIndex);  %����ϱ�����ǵڼ�����ֵ
                  bottomRelativePeak = max(maxNumIndex,max2NumIndex); %����±�����ǵڼ�����ֵ
                  surfaceRelativeLocation = JudgeBig(surfaceRelativePeak); %��ȡ�ϱ����ֵ�����λ��
                  bottomRelativeLocation = JudgeBig(bottomRelativePeak);   %��ȡ�±����ֵ�����λ��
                  surfaceIntensity = IntensityBig(surfaceRelativePeak); %��ȡ�ϱ����ǿ��ֵ
                  bottomIntensity = IntensityBig(bottomRelativePeak);   %��ȡ�±����ǿ��ֵ
                  
                  %�Է�ֵǰ���������Ӧ����            
                  InspectRange =round( abs(surfaceRelativeLocation - bottomRelativeLocation) / 2); %������ = ������ֵ��һ�����
%                     InspectRange = 5;
%�ϱ�����                 
                  %�ж��ڼ����������ļ�����ֵ��
                  surfaceRangeMin = round(surfaceRelativeLocation - InspectRange);
                  surfaceRangeMax = round(surfaceRelativeLocation + InspectRange);
                  [~,n] = size(JudgeBig);
                  flagChange = 1;
                  for surfaceRangeIndex = 1:n
                      if flagChange == 1
                          if JudgeBig(surfaceRangeIndex) >= surfaceRangeMin   %�����ֵ�ڼ�����ķ�Χ
                            surfaceInspectNumMin = surfaceRangeIndex  ;  %ȷ����鷶Χ�ڵ���ʼ��ֵ�ǵڼ�����ֵ
                            flagChange = 2;
                          end
                      else
                          if flagChange == 2
                              if JudgeBig(surfaceRangeIndex) > surfaceRangeMax   %�����ֵ�ڼ�����ķ�Χ�⣬�������
                                surfaceInspectNumMax = surfaceRangeIndex - 1; %ȷ����鷶Χ�ڵ���ֹ��ֵ�ǵڼ�����ֵ
                                break;
                              end
                          end
                      end                    
                  end
                  %�ڼ������ڵĽϴ�ֵ��ֵ��ǿ��ֵ
                  [surfaceInspectIntensity,surfaceInspectRelativePeak] = max(SeekBigIntensity(surfaceInspectNumMin:surfaceInspectNumMax));
%                   surfaceInspectIntensity   %�ϴ��ֵ��ǿ��ֵ
                  surfaceInspectRelativePeak = surfaceInspectRelativePeak + surfaceInspectNumMin - 1; %�ϴ��ֵ�ǵڼ�����ֵ
                  surfaceInspectRelativeLocation = JudgeBig(surfaceInspectRelativePeak);      %�ϴ��ֵ�����λ��  
                  
                  %���ϱ����ǿ�Ȳ�ֵΪ
                  surfaceDifference = surfaceIntensity - surfaceInspectIntensity;
                  %�ϴ��ֵ���ϱ����ֵ�Ĳ�ֵ  ��ռ���ϱ����ֵ�İٷֱ�
                  surfacePercentage =surfaceDifference / surfaceIntensity * 100;
%                   surfacePercentage =vpa(surfaceDifference / surfaceIntensity*100);
%                   disp(['surfacePercentage=',char(surfacePercentage),'%'])   %����ٷֱ�

                  %�ж��Ƿ���Ҫ�滻λ�� 
                  surfaceInspectBiggerNum = 1;
                  if (surfacePercentage <= 20)  && (surfaceRelativeLocation - surfaceInspectRelativeLocation >= 0)   %���С��20% �����ڷ�ֵ֮ǰ��������滻
                      surfaceRelativeLocation = surfaceInspectRelativeLocation ;
                      surfaceIntensity = surfaceInspectIntensity;
                      surfaceInspectBiggerNum = surfaceInspectBiggerNum + 1;  %�ж�����Ҫ�滻��A-Scan
                  end
% % % �±�����
% %                   �ж��ڼ����������ļ�����ֵ��
%                   bottomRangeMin = round(bottomRelativeLocation - InspectRange);
%                   bottomRangeMax = round(bottomRelativeLocation + InspectRange);
% %                   [~,n] = size(SeekBigIndex); %�ظ���
%                   flagChange = 1;
%                   for bottomRangeIndex = 1:n
%                       if flagChange == 1
%                           if JudgeBig(bottomRangeIndex) >= bottomRangeMin   %�����ֵ�ڼ�����ķ�Χ
%                             bottomInspectNumMin = bottomRangeIndex;    %ȷ����鷶Χ�ڵ���ʼ��ֵ�ǵڼ�����ֵ
%                             flagChange = 2;
%                           end
%                       else
%                           if flagChange == 2
%                               if JudgeBig(bottomRangeIndex) >= bottomRangeMax   %�����ֵ�ڼ�����ķ�Χ�⣬�������
%                                 bottomInspectNumMax = bottomRangeIndex - 1; %ȷ����鷶Χ�ڵ���ֹ��ֵ�ǵڼ�����ֵ
%                                 break;
%                               end
%                           end
%                       end                    
%                   end
%                   %�ڼ������ڵĽϴ�ֵ��ֵ��ǿ��ֵ
%                   [bottomInspectIntensity,bottomInspectNumIndex] = max(SeekBigIntensity(bottomInspectNumMin:bottomInspectNumMax));
% %                   bottomInspectIntensity   %�ϴ��ֵ��ǿ��ֵ
%                   bottomInspectRelativePeak = bottomInspectNumIndex + bottomInspectNumMin - 1; %�ϴ��ֵ�ǵڼ�����ֵ
%                   bottomInspectRelativeLocation = JudgeBig(bottomInspectRelativePeak);      %�ϴ��ֵ�����λ��  
%                   
%                   %���ϱ����ǿ�Ȳ�ֵΪ
%                   bottomDifference = bottomIntensity - bottomInspectIntensity;
%                   %�ϴ��ֵ���±����ֵ�Ĳ�ֵ  ��ռ���±����ֵ�İٷֱ�
%                   bottomPercentage =bottomDifference / bottomIntensity*100;
% %                   bottomPercentage =vpa(bottomDifference / bottomIntensity*100);
% %                   disp(['bottomPercentage=',char(bottomPercentage),'%'])              % ����ٷֱ�
% 
%                   %�ж��Ƿ���Ҫ�滻λ�� 
%                   bottomInspectBiggerNum = 1;
%                   if (bottomPercentage <= 20)  & (bottomInspectRelativeLocation >= bottomRelativeLocation)   %���С��20% �����ڷ�ֵ֮��������滻 %�·�ֵ�Ż��㷨
% %                   if (bottomPercentage <= 20)  & (bottomInspectRelativeLocation <= bottomRelativeLocation)   %���С��20% �����ڷ�ֵ֮ǰ��������滻  %�Ϸ�ֵ�Ż��㷨
%                       bottomRelativeLocation = bottomInspectRelativeLocation ;
%                       bottomIntensity = bottomInspectIntensity;
%                       bottomInspectBiggerNum = bottomInspectBiggerNum + 1;   %�ж�����Ҫ�滻��A-scan
%                   end

                  flag = 0;  %������ѭ��
               else
                   SeekBigIntensity(max2NumIndex) = 0;  %�����10pixel���дδ�ֵ���֣���Ϊ����0
               end             
           end
           
           surfaceRelativeLocation = surfaceRelativeLocation+rangeMin-1;     %�ϱ���ķ�ֵλ��
           bottomRelativeLocation = bottomRelativeLocation+rangeMin-1;      %�±���ķ�ֵλ��
           depthNum = (abs(surfaceRelativeLocation-bottomRelativeLocation))*0.003493/1.38*10^3;   %���ֵ��ȷ��
            %����A-Scan��ѭ�����������������ܽ�   
            %���ֵ�������j--����ǰ���Ϊ0
%             surfacePercentage(Index) = surfacePercentage;
%             bottomPercentage(Index) = bottomPercentage;

        else   %�������λ��Ϊ0��������������ΧΪ0
%             problem = 0; %ǰ�ڳ�������Ų�
            surfaceRelativeLocation = nan; 
            bottomRelativeLocation = nan;   
            surfaceIntensity = nan; 
            bottomIntensity = nan ; 
            depthNum = nan;
%             surfacePercentage = nan;
%             bottomPercentage = nan;
        end
            %�����ݽ��д���
            surfaceLocation(Index,IndexK) = surfaceRelativeLocation; %��ȡ�ϱ����ֵ��λ��
            bottomLocation(Index,IndexK) = bottomRelativeLocation;  %��ȡ�±����ֵ��λ��
            surfaceIntensityA(Index,IndexK) = surfaceIntensity; %��ȡ�ϱ����ǿ��ֵ
            bottomIntensityA(Index,IndexK) = bottomIntensity;   %��ȡ�±����ǿ��ֵ
            depth(Index,IndexK) = depthNum;  %��ȡƤ�����ֵ
%             surfacePercentageA(Index,IndexK) = surfacePercentage; %Ϊ�˻�ȡ�ϱ������ֵ��Χ
%             bottomPercentageA(Index,IndexK) =  bottomPercentage; %Ϊ�˻�ȡ�ϱ������ֵ��Χ
    end
        %һ��B-Scan������A-scanѭ������
        % ��A-Scanͼ����м���ֵ��Сֵ���б�ע
%         hold on
%         scatter(JudgeBig,IntensityBig,'*');  %��漫��ֵ��
%         scatter(JudgeSmall,IntensitySmall,'^'); %��漫Сֵ��
%         plot([surfaceRelativeLocation,surfaceRelativeLocation],[maxPeak,max2Peak],'color','r','linewidth',2);   %��ʼ������
%         plot([bottomRelativeLocation,bottomRelativeLocation],[maxPeak,max2Peak],'color','r','linewidth',2);   %��ֹ������     
%         hold off   
        
        %���±�������λ�ú�ǿ��
%         surfaceIntensity
%         bottomIntensity
%         surfaceRelativeLocation
%         bottomRelativeLocation
        
%           %������ֵ��Ƶ��ͳ��
%           tabulate(surfacePercentageA(:))
%           fig2 = figure(2),clf;
%           histfit(surfacePercentageA);
%           i = 1;
%           print(fig2,'-dbitmap',strcat('D:\ϵͳ�Դ�\����\ProgramChange\Percentage','i','SurfacePercentage.bmp'));
%           i = i +1;
   
    waitbar(k/line,wbSkin3dProcess,['����Ƥ�����...' num2str(100*k/line) '%']);
end
close(wbSkin3dProcess);

%             SuanZi=fspecial('average',[3 3]); %% ����һ���˲��� 
%             averagesurfaceLocation=imfilter(surfaceLocation,SuanZi,'replicate');
%             averagebottomLocation=imfilter(bottomLocation,SuanZi,'replicate');
% 

% % ���������������Ż�
% [x,y]=size(depth)
% for i = 15:x-15
%     for j =15:y-15
%         if depth(i,j)-depth(i-1,j-1)>15
%             depth(i,j) = round((depth(i,j-10)+depth(i,j+10))/2);
%         end
%     end
% end

            %�鿴B-Scan�����±����ʶ���
%             [m,n] = size(surfaceLocation); 
%             fig1 = figure(1);hold on
%             for juli = 1:m
%                 surfaceEnd(juli+centerX-r-1) = averagesurfaceLocation(juli,line/2-centerY+r+1 );
%                 bottomEnd(juli+centerX-r-1) = averagebottomLocation(juli,line/2-centerY+r+1 );   
%             end
%             plot(surfaceEnd,'r');
%             plot(bottomEnd,'r');
%             hold off


%��֤�����Ƿ���ɹ�
% A = surfaceLocation(400,400)
% B = surfaceIntensityA(400,400)
% C = bottomLocation(400,400)
% D = bottomIntensityA(400,400)
end
%����ٷֱȵ�ģ��
%a=vpa(0.45);
%disp(['a=',char(a),'%'])