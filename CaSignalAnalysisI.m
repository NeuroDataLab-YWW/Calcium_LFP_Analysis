  clc;
 clearvars;
 close all;
 %%
clearvars
[filename, pathname] = uigetfile({'*.mat';'*.dat';'*.m';'*.*'},'File Selector');
a1=pathname;
a2=filename;
c=[a1,a2];
a = importdata(c);
y = a.values;
rs=round(1/a.interval);
% y=data2;
rs=100;
x=1:length(y);
x= x/rs;
plot(x,y);
%%
Base_start1=input('Please input base: ');  
Base_end1=input('Please input base: '); 
Time_start1=input('Please input : ');   
Time_start2=input('Please input : ');  
Time_end2=input('Please input : ');
Time_inject=input('Please input : ');
Time_inject=Time_inject/60;
x1=x(Time_start1*rs:Time_end1*rs); 
x2=x(Time_start2*rs:Time_end2*rs); 
y1=y(Time_start1*rs:Time_end1*rs)';
y2=y(Time_start2*rs:Time_end2*rs)';
XX=x(Time_start1*rs:Time_end2*rs);
YY=y(Time_start1*rs:Time_end2*rs);
Baseline=mean(y(Base_start1*rs:Base_end1*rs));
X=[x1 x2];
Y=[y1 y2];
figure(2)
plot(X,Y);
%%

Style=1;
color=[1 0 0;0 1 0;0 0 1;0 1 1];
n=2;
[p,s,mu] = polyfit(X,Y,n); 
fit_y = polyval(p,XX,[],mu); 
Base_F=mean(fit_y);
dt_y= YY-fit_y';
F=dt_y*100/(Base_F-Baseline);
Min=floor(min(F)/10)*10;
Max=ceil(max(F)/10)*10;
figure(3);
subplot(2,1,1);plot(XX,YY);
hold on;
plot(XX,fit_y);
subplot(2,1,2);plot(XX,F);
figure(4)
switch Style
    case 1
plot(XX/60,smooth(F,100),'red');
    case 2
        point=(Time_inject*60-XX(1))*rs+1;
        f=F(1:point(1));       
        plot(XX(1:point(1))/60,smooth(f,100),'color',color(1,:));
        FF(1).f=f;
        hold on
        for i=1:length(Time_inject)
            if i==length(Time_inject)
                f=F(point(i):end);
                plot(XX(point(i):end)/60,smooth(f,100),'color',color(end,:))              
            else
                 f=F(point(i):point(i+1));
                 plot(XX(point(i):point(i+1))/60,smooth(f,100),'color',color(i+1,:));
            end
            FF(i+1).f=f;
        end        
end
hold on
for i=1:length(Time_inject)
line([Time_inject(i) Time_inject(i)],[Min Max],'color','k','linewidth',1.5,'linestyle','--');
end
xlabel('Time/min','FontSize',12,'FontWeight','Bold');
ylabel('\Delta F/F(%)','FontSize',12,'FontWeight','Bold');
set(gca,'box','off','tickdir','out','linewidth',1.5,'fontsize',12);
ylim([Min Max]);
xlim([min(XX/60) max(XX)/60])

%%
time_in = Time_inject*60;
data3_4=F((time_in-600-Time_start1)*100:(time_in+60*60-Time_start1)*100);
smooth(F,200);

data3_4=data3_4';

figure
plot(data3_4)
before_3_4=mean(data3_4(100*100:550*100));
after_3_4=mean(data3_4(13*60*100:17*60*100));

