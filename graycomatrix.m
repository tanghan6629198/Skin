function [ASM,ENT,COR] = graycomatrix(Offset,I,N) 
%     Offset = [0,3];
%     N = 8;
    dx = Offset(1,1);
    dy = Offset(1,2);
%     Imax = max(max(I));
%     Imin = min(min(I));
    Imax = max(max(I(1:20,:)));
    Imin = min(min(I(1:20,:)));
%     a = size(I,1);
%     b = size(I,2); 
    a = 20;
    b = 100; 
    glcms = zeros(a,b);
    GLCMS = zeros(N,N);
    for i = 1:a
        for j = 1:b
            if I(i,j) == Imax
                glcms(i,j) = N;
            else
            glcms(i,j) = floor((I(i,j) - Imin) / (Imax - Imin) * N) + 1 ;  %%N 灰度级数
            end
        end
    end
    figure();imagesc(glcms); axis image; colormap('gray'); hold on;
    
    for i = 1-dx:a
        for j = 1:b-dy
            Px = glcms(i,j);
            Py = glcms(i+dx,j+dy);
            GLCMS(Px,Py) = GLCMS(Px,Py) + 1;
        end
    end   
    Psum = (a+dx)*(b-dy);
    
    for i = 1:N
        for j = 1:N
            GLCMS(i,j) = GLCMS(i,j)/Psum;
        end
    end
%%%&&&&&&&&&&&&&&&&
    ASM = sum(sum(GLCMS.^2)) ;%能量 
%%%&&&&&&&&&&&&&&&&
    GLCMSlog =log(GLCMS); 
    ENT = 0;
    for i = 1:N
        for j = 1:N
            if GLCMS(i,j) ~= 0
            ENT = ENT + (-GLCMS(i,j)*GLCMSlog(i,j)); %熵
            end
        end
    end
    ENT;
%%%&&&&&&&&&&&&&&&&
    COR = 0;ux = 0;uy = 0;deltax = 0;deltay = 0;COR2 = 0;
    for i = 1:N
        for j = 1:N
            ux = i * GLCMS(i,j) + ux;   %相关性中的μx
            uy = j * GLCMS(i,j) + uy;   %μy
        end
    end
    
    for i = 1:N
        for j = 1:N
            deltax = (i - ux)^2 * GLCMS(i,j) + deltax;   %相关性中的μx
            deltay = (i - uy)^2 * GLCMS(i,j) + deltay;   %μy
            COR = (i - ux) * (j - uy) * GLCMS(i,j) + COR ;
%            COR2 = i*j*GLCMS(i,j) + COR2;
            
        end
    end
    COR = COR/deltax/deltay;
%    COR2 = (COR2 - (ux * uy))/deltax/deltay
    
end 