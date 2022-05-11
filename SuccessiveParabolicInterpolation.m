%--------------------------------------------------------------------------------
% auther:Zhu Binglong
% Date:2022-5-11 15:58:09
% Function: Successive Parabolic Interpolation + Axesion transformation
% Reference: https://mathsfromnothing.cf/successive-parabolic-interpolation/?i=1
%--------------------------------------------------------------------------------
clear all;
clc;
close all;
format long;

figure(1)

% anonymous function
% [x,y]=fminbnd(@(x)x.^3-13*x.^2+20*x+100,-5,20)
% [x,y]=fminbnd(@(x)x.^2/10-2*sin(x),-3,5)

%% function 1
%  f=@(x)x.^3-13*x.^2+20*x+100;% best_position= 7.813437226365976,best_fitness= -60.369882183024487
% fplot(f,[-5 20]);hold on;
% [x,fo]=successiveparabolicinterpolation(f,12);

%% function 2 :https://github.com/gopherclass/successive-parabolic-interpolation
f=@(x)x.^2/10-2*sin(x);%best_position= 1.427551777634208 ,best_fitness=-1.775725653147415,
fplot(f,[-15 15]);hold on;
[x,fo]=successiveparabolicinterpolation(f,3);
saveas( 1, 'Result.jpg')

figure(2)
plot(1:1:numel(fo),fo,'*')
saveas( 2, 'FimWithIter.jpg')

disp(['best_position:' num2str(x(end),'% 20.16g') ',best_fitness:' num2str(fo(end),'% 20.16g')] );
% x(end)
% fo(end)
%% sub function
function [x,fo]=successiveparabolicinterpolation(f,x)
% https://mathsfromnothing.cf/successive-parabolic-interpolation/?i=1
%Example input to the Successive parabolic interpolation optimization method
%[x,fc,normc]=successiveparabolicinterpolation(f,x)
%f=@(x)x.^3-13*x.^2+20*x+100;
%x=12000000;
%
%Input
%f is a function handle of the function to be minimized
%x is the current position
%
%Output
%x is the new location
%fc is the function evaluated at the new location

xtol=1e-10;
ytol=1e-8;
delt=1;% trust domain radius
beta=1;% Maximum axial transformation radius
xc=[x-delt x x+delt];% trust domain
fc=f(xc);% solve Error value
plot(xc,fc,'o','markerfacecolor','w');% plot the first point
for k=1:20
    k % let's see if  k is properly
    %     if(xc(1)==xc(2) || xc(1)==xc(3) || xc(2)==xc(3))
    %         break;
    %     end
    xnew=(fc(1)*(xc(2)^2-xc(3)^2)-fc(2)*(xc(1)^2-xc(3)^2)+fc(3)*(xc(1)^2-xc(2)^2))/(2*(fc(1)*(xc(2)-xc(3))-fc(2)*(xc(1)-xc(3))+fc(3)*(xc(1)-xc(2))));
    if isnan(xnew)||isinf(xnew)
        xnew=xc(1);
    end
    ynew=f(xnew);
    plot(xnew,ynew,'o','markerfacecolor','w');hold on;%plot the vertex of the parabola
    if any(abs(xnew-xc)<xtol)||all(abs(ynew-fc)<ytol)
        break;% terminal condition
    end
    if ynew<=min(fc)% if the  fmin value gets smaller
        [maxfc,maxfcI]=sort(fc);% ascending order;
        xc(maxfcI(end))=xnew;
        fc(maxfcI(end))=ynew;% the biggest is instead
        xbest=xnew;
        ybest=ynew;
    else%如果找到的是最大值xnew
        
        %         %用新找的误差最大的和已有误差最小的点进行状态转移 不好，会有反复
        %         [maxfc,maxfcI]=sort(fc);
        %         xk=xc(maxfcI(1));%误差最小的x
        %         xkm1=xnew;%x(k-1) %误差最大的x
        
        %用已有的三个点中误差最大的和误差最小的进行转移操作
        [maxfc,maxfcI]=sort(fc);
        xk=xc(maxfcI(1));%已有的误差最小的x
        xkm1=xc(maxfcI(end));%x(k-1) %误差最大的x
        xkp1=xk+beta*rand*(xk-xkm1)/norm(xk-xkm1);%x(k+1) 转移算子
        
        xnew=xkp1;%找到的应该是更好的解
        ynew=f(xnew);
        
        
        if ynew<=min(fc)
            plot(xnew,ynew,'o','markerfacecolor','w');hold on;
            %             drawarrow([xkm1, xnew], [fc(maxfcI(end)), ynew])
            drawarrow([xkm1, fc(maxfcI(end))], [xnew, ynew])
            %             PlotLineArrow(gca, [xkm1, xnew], [fc(maxfcI(end)), ynew],'b', 'r');
            
            [maxfc,maxfcI]=sort(fc);%升序
            xc(maxfcI(end))=xnew;
            fc(maxfcI(end))=ynew;% the biggest is instead
            xbest=xnew;
            ybest=ynew;
        end
        
    end
    %     [xc,xcI]=sort(xc);
    %     fc=fc(xcI);
    plot(xc,fc,'o','markerfacecolor','w');hold on;
    [a,b,c]=parabola(xc,fc);%当前抛物线绘图
    x(k)=xbest;
    fo(k)=ybest;
end

plot(xbest(end),ybest(end),'*');hold on;% the best point with minimum f

end

% 绘制抛物线，并给出系数a b c
function [a,b,c]=parabola(xc,fc)
x1=xc(1);x2=xc(2);x3=xc(3);
y1=fc(1);y2=fc(2);y3=fc(3);
% x1=p1(1);y1=p1(2);
% x2=p2(1);y2=p2(2);
% x3=p3(1);y3=p3(2);
if x1~=x2 && x2~=x3 && x1~=x3
    a=-((-x2* y1+x3 *y1+x1* y2-x3* y2-x1* y3+x2* y3)/((-x1+x2) *(x2-x3) *(-x1+x3)));
    b=-((x2^2 *y1 - x3^2 *y1 - x1^2 *y2 + x3^2 *y2 + x1^2 *y3 - x2^2 *y3)/((x1 - x2) *(x1 - x3) *(x2 - x3)));
    c=-((-x2^2 *x3 *y1+x2 *x3^2 *y1+x1^2 *x3 *y2-x1 *x3^2 *y2-x1^2 *x2 *y3+x1 *x2^2 *y3)/((x1-x2)* (x1-x3) *(x2-x3)));
else
    a=0;b=0;c=0;
end

% P=[a,b,c];
% x=min(px):0.1:max(px);
px=xc;
x=linspace(min(px)-(max(px)-min(px))/0.5,max(px)+(max(px)-min(px))/0.5);
y=a*x.^2+b*x+c;
plot(x,y,'--');
hold on;
end

function drawarrow(x,y,lineType,ax)
switch nargin
    case 2
        lineType='arrow';
        ax=gca;
    case 3
        ax=gca;
end
% 调整坐标大小以适应箭头长度
xlim=ax.XLim;
ylim=ax.YLim;
xlimmin=xlim(1);xlimmax=xlim(2);
ylimmin=ylim(1);ylimmax=ylim(2);
if xlimmin>min(x(1),y(1)), xlimmin=min(x(1),y(1));end
if xlimmax<max(x(1),y(1)), xlimmax=max(x(1),y(1));end
if ylimmin>min(x(2),y(2)), ylimmin=min(x(2),y(2));end
if ylimmax<max(x(2),y(2)), ylimmax=max(x(2),y(2));end
ax.XLim = [xlimmin,xlimmax];
ax.YLim = [ylimmin,ylimmax];
xlim=ax.XLim;
ylim=ax.YLim;
pos=ax.Position;
x_ratio = pos(3)/(xlim(2)-xlim(1));
y_ratio = pos(4)/(ylim(2)-ylim(1)); % 缩放比例
orig_pos=[-xlim(1)*x_ratio+pos(1),-ylim(1)*y_ratio+pos(2)]; % figure坐标系中的原点坐标
x=x.*[x_ratio,y_ratio];y=y.*[x_ratio,y_ratio];
x=x+orig_pos;y=y+orig_pos;
% annotation(lineType,[x(1),y(1)],[x(2),y(2)])
ar=annotation(lineType,[x(1),y(1)],[x(2),y(2)]);
ar.Color='Red';
ar.LineStyle='--';
ar.LineWidth=1;
end
% ————————————————
% 版权声明：本文为CSDN博主「qq_41093957」的原创文章，遵循CC 4.0 BY-SA版权协议，转载请附上原文出处链接及本声明。
% 原文链接：https://blog.csdn.net/qq_41093957/article/details/115256114