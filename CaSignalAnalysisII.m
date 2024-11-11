% %data built up
 datapath = '/Users/FiberPhotometry/';
 lists = dir([datapath '*.mat']);
 k=length(lists);
 data2=[];
%%
for i=1:k
     name=lists(i).name;
a = importdata([datapath name]);
 data3=a;
 qushu= [1:6000:360000]
 data3=data3(:,qushu);
 data3=smooth(data3,1);
 data2=[data2;data3'];
end

%% draw errorline
m=size(data2,1);   
tie=length(data2);
time = [1:tie];
datamean = mean(data2);
datamean=smooth(datamean,100);
sem = std(data2)/sqrt(m);
sem = smooth(sem,100);
% figure;
figure;

drawErrorLine(time,datamean,sem,[0.98 0.47 0.74],[0.04 0.85 0.71],0.2); 


 set(gca, 'FontSize', 25)
xlabel('Time(min)','FontSize',14);
xticks([0 5 10 15 20 25 30 35 40 45 50 55 60 65 70])
 ylabel('\Delta F/F(%)','FontSize',14);
 set(gca, 'FontSize', 14)
 xlabel('Time(min)','FontSize',14);
 ylabel('\Delta F/F(%)','FontSize',14);
 ylim([-20 20]);
box off;
hold on;
%% 
% %data built up
 datapath = '/Users/FiberPhotometry/';
 lists = dir([datapath '*.mat']);
 k=length(lists);
 data2=[];
%%
for i=1:k
     name=lists(i).name;
a = importdata([datapath name]);
 data3=a;
 qushu= [1:6000:24000];
 data3=data3(:,qushu);
 data3=smooth(data3,1);
 data2=[data2;data3'];
end

%% draw errorline
m=size(data2,1);   
tie=length(data2);
time = [1:tie];
datamean = mean(data2);
datamean=smooth(datamean,100);
sem = std(data2)/sqrt(m);
sem = smooth(sem,100);

drawErrorLine(time,datamean,sem,'k','k',0.1);  %[0.98 0.47 0.74],[0.04 0.85 0.71],0.2


 set(gca, 'FontSize', 25)
xlabel('Time(min)','FontSize',14);
xticks([0 5 10 15 20 25 30 35 40 45 50 55 60 65 70])
 ylabel('\Delta F/F(%)','FontSize',14);
 set(gca, 'FontSize', 14)
 xlabel('Time(min)','FontSize',14);
 ylabel('\Delta F/F(%)','FontSize',14);
 ylim([-20 20]);
box off;
hold off;