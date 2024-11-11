
fs = 30000;
for m =1:8
    data=[];XY=[];power=[];
    data = Bandpass_Filter_CH8_Ch1.values((200+(m-1)*10)*fs:(205+(m-1)*10)*fs); 
    [~,data_wc1] = buttord([48.5 51.5]*2/fs,[49.9 50.1]*2/fs,0.1,100,'s');
    [data_b1,data_a1] = butter(1,data_wc1,'stop');
    data = filtfilt(data_b1,data_a1,data);


    for j = 1:size(data,2)
        figure(m);
        subplot(3,1,1)
        plotHz(fs, data(:,j));
        subplot(3,1,2)
        params = struct();
        params.Fs = fs;
        params.tapers = [3, 5]; 
        params.fpass = [1 100];
        fpass = [1 100];
        movingwin = [4, 2];
        [S,t,f]=mtspecgramc(data(:,j),movingwin,params);
        plot_matrix(S,t,f, 'l');
        subplot(3,1,3)
        power(:,j) = mean(S,1)';
        plot(f, power(:,j));
    end

    f=f'; 
    XY=[f,power];
    allgamma_power_mean(m)=mean(XY(106:390,2));
    lowgamma_power_mean(m)=mean(XY(106:214,2));
    highgamma_power_mean(m)=mean(XY(215:390,2));
    theta_power_mean(m)=mean(XY(14:35,2));
    alpha_power_mean(m)=mean(XY(36:52,2));
    Beta_power_mean(m)=mean(XY(62:105,2));
    Delta_power_mean(m)=mean(XY(1:13,2));
    gammamean_2535(m)=mean(XY(106:148,2));
end
allgamma_power_mean(:,all(allgamma_power_mean==0,1))=[];lowgamma_power_mean(:,all(lowgamma_power_mean==0,1))=[];highgamma_power_mean(:,all(highgamma_power_mean==0,1))=[];theta_power_mean(:,all(theta_power_mean==0,1))=[];alpha_power_mean(:,all(alpha_power_mean==0,1))=[];Beta_power_mean(:,all(Beta_power_mean==0,1))=[];Delta_power_mean(:,all(Delta_power_mean==0,1))=[];gammamean_2535(:,all(gammamean_2535==0,1))=[];
allgammamean=mean(allgamma_power_mean);lowmean=mean(lowgamma_power_mean);highmean=mean(highgamma_power_mean);
thetamean=mean(theta_power_mean);alphamean=mean(alpha_power_mean);betamean=mean(Beta_power_mean);gammamean_2535=mean(gammamean_2535);
Deltapowermean=mean(Delta_power_mean);
alldata={'allgammamean','lowmean','highmean','thetamean','alphamean','betamean','2535gammamean','Delta_power_mean';allgammamean,lowmean,highmean,thetamean,alphamean,betamean,gammamean_2535,Deltapowermean};


%% 
fs = 30000;
for m = 1:5
    data=[];XY=[];power=[];
    data = Bandpass_Filter_CH11_Ch1.values((845+m*20)*fs:(835+m*20)*fs); 
    [~,data_wc1] = buttord([48.5 51.5]*2/fs,[49.9 50.1]*2/fs,0.1,100,'s');
    [data_b1,data_a1] = butter(1,data_wc1,'stop');
    data = filtfilt(data_b1,data_a1,data);

    for j = 1:size(data,2) 
        figure(m);
        subplot(3,1,1)
        plotHz(fs, data(:,j));
        subplot(3,1,2)
        params = struct();
        params.Fs = fs;
        params.tapers = [3, 5]; 
        params.fpass = [1 100];
        fpass = [1 100];
        movingwin = [4, 2];
        [S,t,f]=mtspecgramc(data(:,j),movingwin,params);
        plot_matrix(S,t,f, 'l');
        subplot(3,1,3)
        power(:,j) = mean(S,1)';
        plot(f, power(:,j));
     end
    f=f'; 
    XY=[f,power];
    allgamma_power_mean(m)=mean(XY(106:390,2));
    lowgamma_power_mean(m)=mean(XY(106:214,2));
    highgamma_power_mean(m)=mean(XY(215:390,2));
    theta_power_mean(m)=mean(XY(14:35,2));
    alpha_power_mean(m)=mean(XY(36:52,2));
    Beta_power_mean(m)=mean(XY(62:105,2));
    gammamean_2535(m)=mean(XY(106:148,2));
end
allmean=mean(allgamma_power_mean);lowmean=mean(lowgamma_power_mean);highmean=mean(highgamma_power_mean);
thetamean=mean(theta_power_mean);alphamean=mean(alpha_power_mean);betamean=mean(Beta_power_mean);gammamean_2535=mean(gammamean_2535);
alldata={'allgammamean','lowmean','highmean','thetamean','alphamean','betamean','2535gammamean';allmean,lowmean,highmean,thetamean,alphamean,betamean,gammamean_2535};