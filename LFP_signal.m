%%
list_path = '/Users/.xlsx';
res_path = '/Users/';
fs = 30000;   
ds_fs = 500;   
smooth_win = 0.01;      
hp_freq = 0.6;         
bsl_range = [-600, 0];  
stats_range = [300, 600];  
plot_range = [300, 600];  
band_list = {[1, 3], [1, 100]};   
band_name = {'1-3 Hz', 'total'}; 
fft_win = 4;    
fft_step = 2;   
grp_clr = [ ... 
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980; ...
    0.9290, 0.6940, 0.1250; ...
    0.4940, 0.1840, 0.5560; ...
    0.6350, 0.0780, 0.1840; ...
    ];

ds = LFP_signal_normalization_bandpower(list_path, res_path, fs, ...
    ds_fs, smooth_win, hp_freq, bsl_range, stats_range, plot_range, ...
    band_list, band_name, fft_win, fft_step, grp_clr);