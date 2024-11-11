function ds = LFP_signal_normalization_bandpower(list_path, res_path, fs, ...
    ds_fs, smooth_win, hp_freq, bsl_range, stats_range, plot_range, ...
    band_list, band_name, fft_win, fft_step, grp_clr)

mkdir(res_path);
ds = table2struct(readtable(list_path));
if isempty(ds_fs); ds_fs = fs; 
fs0 = fs;
fs = ds_fs;
tx = (floor(min(plot_range)):1/fs:ceil(max(plot_range)))';
tx2 = (floor(min(plot_range)/300):1:ceil(max(plot_range))/300)'*300;
if ~isempty(hp_freq)
    try
        hhigh = design(fdesign.highpass('N,F3db', 6, hp_freq, fs), 'butter');
    catch
        hp_freq = [];
    end
end
hnotch = design(fdesign.notch(10, 50, 20, fs), 'butter');
grp_name = unique(arrayfun(@(x)(x.group), ds, 'UniformOutput', false));
ngrp = length(grp_name);
nband = length(band_name);

for nn = 1:length(ds)
    [~, fname] = fileparts(ds(nn).path);
    fname = strrep(fname, '.', '_');
    oname = sprintf('%.3d_%s_%s', nn, ds(nn).group, fname);
    fprintf('[%.3d/%.3d] %s >>>> ', nn, length(ds), ds(nn).path);

    try
        sig = load(ds(nn).path);
        vls = fields(sig);
        sig = sig.(vls{1});
        if isfield(sig, 'values')
      sig = sig.values;
        else
            sig = reshape(sig, [], 1);
        end
    catch
        fprintf('file error!\n');
        continue;
    end
    ds_fac = fs0 / fs;
    sig = double(sig(1:ds_fac:end));
    if fs > 100
        sig = filtfilt(hnotch.sosMatrix, hnotch.ScaleValues, sig);
    end
    nt = length(sig);
    ts = (0:1:nt-1)' / fs - ds(nn).time;

    is_bsl = ts > min(bsl_range) & ts < max(bsl_range);
    sig = (sig - mean(sig(is_bsl))) / std(sig(is_bsl));
    if ~isempty(hp_freq)
        sig = filtfilt(hhigh.sosMatrix, hhigh.ScaleValues, sig);
    end
    sig_sm = smooth(sig, smooth_win * fs);
    
    psd = cell(nband, 1);
    for kk = 1:nband
        [pxx, ~, txx] = mspectrum(sig, fft_win, fft_step, fs, band_list{kk});
        pxx = mean(pxx, 1) * range(band_list{kk});
        psd{kk} = interp1(txx - ds(nn).time, pxx, tx, 'linear', nan);
    end
    psd = cat(2, psd{:});

    is_stats = ts > min(stats_range) & ts < max(stats_range);
    [pxx, fxx] = pwelch(sig(is_stats), [], [], 2^nextpow2(nt), fs);
    bp = nan(1, nband);
    for kk = 1:nband
        in_freq = (fxx > min(band_list{kk})) & (fxx <= max(band_list{kk}));
        bp(kk) = mean(pxx(in_freq)) * range(band_list{kk});
    end  

    ds(nn).name = oname;
    ds(nn).sig = interp1(ts, sig_sm, tx, 'nearest', nan);
    ds(nn).psd = psd;
    ds(nn).bp = bp;

    figure; set(gcf, 'Position', [50, 50, 800, 300], 'Visible', 'off', ...
        'PaperOrientation', 'landscape'); 
    plot(ts, sig_sm, '-k', 'LineWidth', 1);
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12, ...
        'XLim', plot_range, 'XTick', tx2);
    xlabel('Time to treatment (s)', 'FontSize', 12);
    ylabel('Normalized LFP', 'FontSize', 12);
    title(ds(nn).path, 'Interpreter', 'none', 'FontSize', 12);    
    print(gcf, fullfile(res_path, oname), '-djpeg', '-r300');
    print(gcf, fullfile(res_path, oname), '-dpdf', '-r300', '-vector');
    close;

    fprintf('%s finished.\n', oname);
end

var_names = { ...
    arrayfun(@(x)(x.name), ds, 'UniformOutput', false), band_name, {'bandpower'}, ...
    };
write_excel_matrix(cat(1, ds(:).bp), var_names, 'name', 1, 2, 3, ...
    fullfile(res_path, 'normalized_bandpower.xlsx'));

figure; set(gcf, 'Position', [50, 50, 800, 300], ...
    'PaperOrientation', 'landscape'); hold on;
for ii = 1:ngrp
    in_ses = arrayfun(@(x)(strcmp(x.group, grp_name{ii})), ds);
    temp = arrayfun(@(x)(x.sig), ds(in_ses), 'UniformOutput', false);
    temp = cat(2, temp{:});
    c = grp_clr(ii, :) + (1 - grp_clr(ii, :)) * (1 - 0.3);
    plot(tx, temp, 'Color', c, 'LineWidth', 0.3);
end
hs = zeros(ngrp, 1);
for ii = 1:ngrp
    in_ses = arrayfun(@(x)(strcmp(x.group, grp_name{ii})), ds);
    temp = arrayfun(@(x)(x.sig), ds(in_ses), 'UniformOutput', false);
    temp = mean(cat(2, temp{:}), 2, 'omitnan');
    hs(ii) = plot(tx, temp, 'Color', grp_clr(ii, :), 'LineWidth', 1);
end
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12, ...
        'XLim', plot_range, 'XTick', tx2);
legend(hs, grp_name, 'Box', 'off');
ylim([-4 4]);
xlabel('Time to treatment (s)', 'FontSize', 12);
ylabel('Normalized LFP', 'FontSize', 12);
title('Normalized LFP traces', 'Interpreter', 'none', 'FontSize', 12);
print(gcf, fullfile(res_path, 'LFP_traces_individual'), '-djpeg', '-r300');
print(gcf, fullfile(res_path, 'LFP_traces_individual'), '-dpdf', '-r300', '-vector');
close;

[b, a] = butter(4, [1 3] / (fs / 2), 'bandpass');


for nn = 1:length(ds)
    try
        sig = ds(nn).sig;
        ds(nn).filtered_sig = filtfilt(b, a, sig); 
    catch
        fprintf('Error in filtering data for dataset %d\n', nn);
    end
end

figure; set(gcf, 'Position', [50, 50, 800, 300], 'PaperOrientation', 'landscape'); hold on;

for ii = 1:ngrp
    in_ses = arrayfun(@(x)(strcmp(x.group, grp_name{ii})), ds);
    temp = arrayfun(@(x)(x.filtered_sig), ds(in_ses), 'UniformOutput', false);
    temp = cat(2, temp{:});
    c = grp_clr(ii, :) + (1 - grp_clr(ii, :)) * (1 - 0.3);
    plot(tx, temp, 'Color', c, 'LineWidth', 0.3);
end

hs = zeros(ngrp, 1);
for ii = 1:ngrp
    in_ses = arrayfun(@(x)(strcmp(x.group, grp_name{ii})), ds);
    temp = arrayfun(@(x)(x.filtered_sig), ds(in_ses), 'UniformOutput', false);
    temp = mean(cat(2, temp{:}), 2, 'omitnan');
    hs(ii) = plot(tx, temp, 'Color', grp_clr(ii, :), 'LineWidth', 1);
end

set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12, 'XLim', plot_range, 'XTick', tx2);
legend(hs, grp_name, 'Box', 'off');
ylim([-1.5 1.5]);
xlabel('Time to treatment (s)', 'FontSize', 12);
ylabel('Normalized LFP (1-3 Hz)', 'FontSize', 12);
title('1-3 Hz LFP traces', 'Interpreter', 'none', 'FontSize', 12);

print(gcf, fullfile(res_path, 'LFP_traces_1-3Hz'), '-djpeg', '-r300');
print(gcf, fullfile(res_path, 'LFP_traces_1-3Hz'), '-dpdf', '-r300', '-vector');
close;


%%

[b, a] = butter(4, [30 80] / (fs / 2), 'bandpass');

for nn = 1:length(ds)
    try
        sig = ds(nn).sig; 
        ds(nn).filtered_sig = filtfilt(b, a, sig); 
    catch
        fprintf('Error in filtering data for dataset %d\n', nn);
    end
end

figure; set(gcf, 'Position', [50, 50, 800, 300], 'PaperOrientation', 'landscape'); hold on;

for ii = 1:ngrp
    in_ses = arrayfun(@(x)(strcmp(x.group, grp_name{ii})), ds);
    temp = arrayfun(@(x)(x.filtered_sig), ds(in_ses), 'UniformOutput', false);
    temp = cat(2, temp{:});
    c = grp_clr(ii, :) + (1 - grp_clr(ii, :)) * (1 - 0.3);
    plot(tx, temp, 'Color', c, 'LineWidth', 0.3);
end

hs = zeros(ngrp, 1);
for ii = 1:ngrp
    in_ses = arrayfun(@(x)(strcmp(x.group, grp_name{ii})), ds);
    temp = arrayfun(@(x)(x.filtered_sig), ds(in_ses), 'UniformOutput', false);
    temp = mean(cat(2, temp{:}), 2, 'omitnan');
    hs(ii) = plot(tx, temp, 'Color', grp_clr(ii, :), 'LineWidth', 1);
end

set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12, 'XLim', plot_range, 'XTick', tx2);
legend(hs, grp_name, 'Box', 'off');
ylim([-1 1]);
xlabel('Time to treatment (s)', 'FontSize', 12);
ylabel('Normalized LFP (30-90 Hz)', 'FontSize', 12);
title('30-90 Hz LFP traces', 'Interpreter', 'none', 'FontSize', 12);

print(gcf, fullfile(res_path, 'LFP_traces_30-90Hz'), '-djpeg', '-r300');
print(gcf, fullfile(res_path, 'LFP_traces_30-90Hz'), '-dpdf', '-r300', '-vector');
close;


%%

figure; set(gcf, 'Position', [50, 50, 800, 300], 'Visible', 'off', ...
    'PaperOrientation', 'landscape'); hold on;
for ii = 1:ngrp
    in_ses = arrayfun(@(x)(strcmp(x.group, grp_name{ii})), ds);
    nses = sum(in_ses);
    temp = arrayfun(@(x)(x.sig), ds(in_ses), 'UniformOutput', false);
    temp = cat(2, temp{:});
    plot_shaded_err(tx, mean(temp, 2, 'omitnan'), std(temp, 1, 2, 'omitnan') / sqrt(nses), ...
        {'Color', grp_clr(ii, :), 'LineWidth', 1});
end
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12, ...
        'XLim', plot_range, 'XTick', tx2);
legend(grp_name, 'Box', 'off');
xlabel('Time to treatment (s)', 'FontSize', 12);
ylabel('Normalized LFP', 'FontSize', 12);
title('Normalized LFP traces', 'Interpreter', 'none', 'FontSize', 12);
print(gcf, fullfile(res_path, 'LFP_traces_group'), '-djpeg', '-r300');
print(gcf, fullfile(res_path, 'LFP_traces_group'), '-dpdf', '-r300', '-vector');
close;

for kk = 1:nband
    figure; set(gcf, 'Position', [50, 50, 800, 300], ...
        'PaperOrientation', 'landscape'); hold on;
    for ii = 1:ngrp
        in_ses = arrayfun(@(x)(strcmp(x.group, grp_name{ii})), ds);
        temp = arrayfun(@(x)(x.psd(:, kk)), ds(in_ses), 'UniformOutput', false);
        temp = cat(2, temp{:});
        c = grp_clr(ii, :) + (1 - grp_clr(ii, :)) * (1 - 0.3);
        plot(tx, temp, 'Color', c, 'LineWidth', 0.3);
    end
    hs = zeros(ngrp, 1);
    for ii = 1:ngrp
        in_ses = arrayfun(@(x)(strcmp(x.group, grp_name{ii})), ds);
        temp = arrayfun(@(x)(x.psd(:, kk)), ds(in_ses), 'UniformOutput', false);
        temp = mean(cat(2, temp{:}), 2, 'omitnan');
        hs(ii) = plot(tx, temp, 'Color', grp_clr(ii, :), 'LineWidth', 1);
    end
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12, ...
        'XLim', plot_range, 'XTick', tx2);
    legend(hs, grp_name, 'Box', 'off');
    xlabel('Time to treatment (s)', 'FontSize', 12);
    ylabel('Power', 'FontSize', 12);
    title(sprintf('%s power traces', band_name{kk}), ...
        'Interpreter', 'none', 'FontSize', 12);
    figname = sprintf('%s_traces_individual', band_name{kk});
    print(gcf, fullfile(res_path, figname), '-djpeg', '-r300');
    print(gcf, fullfile(res_path, figname), '-dpdf', '-r300', '-vector');

    figure; set(gcf, 'Position', [50, 50, 800, 300], 'Visible', 'off', ...
        'PaperOrientation', 'landscape'); hold on;
    for ii = 1:ngrp
        in_ses = arrayfun(@(x)(strcmp(x.group, grp_name{ii})), ds);
        temp = arrayfun(@(x)(x.psd(:, kk)), ds(in_ses), 'UniformOutput', false);
        temp = cat(2, temp{:});
        plot_shaded_err(tx, mean(temp, 2, 'omitnan'), std(temp, 1, 2, 'omitnan') / sqrt(nses), ...
            {'Color', grp_clr(ii, :), 'LineWidth', 1});
    end
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12, ...
        'XLim', plot_range, 'XTick', tx2);
    legend(grp_name, 'Box', 'off');
    xlabel('Time to treatment (s)', 'FontSize', 12);
    ylabel('Power', 'FontSize', 12);
    title(sprintf('%s power traces', band_name{kk}), ...
        'Interpreter', 'none', 'FontSize', 12);
    figname = sprintf('%s_traces_group', band_name{kk});
    print(gcf, fullfile(res_path, figname), '-djpeg', '-r300');
    print(gcf, fullfile(res_path, figname), '-dpdf', '-r300', '-vector');
    close;
end

fprintf('saved in %s.\n', res_path);

end

function [pxxs, fxx, txx] = mspectrum(sig, twin, tstep, fs, freqr)
nt = fix(length(sig)/fs);
ntx = fix(nt/twin);
nts = ntx*fs*twin;
nstep = round(twin/tstep);
tstep = twin/nstep;
sig = cat(1, sig, zeros(nts+(nstep-1)*tstep*fs-length(sig), 1));
txx = (1:ntx*nstep)*tstep+twin/2-tstep;

nfft = 2^nextpow2(twin*fs);
fxx = linspace(0, fs/2, nfft/2+1);
infreq = fxx>freqr(1) & fxx<freqr(2);
fxx = fxx(infreq);
nfreq = sum(infreq);

tapers = permute(dpss(twin*fs, 3, 5)*sqrt(fs), [1, 3, 2]);
pxxs = zeros(nfreq, nstep, ntx);
for kk = 1:nstep
    idx = (1:nts)+round((kk-1)*tstep*fs);
    pxx = fft(reshape(sig(idx), [], ntx).*tapers, nfft)/fs;
    pxx = pxx(infreq,:,:);
    pxxs(:,kk,:) = mean(conj(pxx).*pxx, 3);
end
pxxs = pxxs(:,:);
end

function write_excel_matrix(outmat, labels, rname, rdim, cdim, sdim, filepath)
rnames = labels{rdim};
cnames = labels{cdim};
snames = labels{sdim};
rtbl = table(rnames, 'VariableNames', {rname});
outmat = reshape(permute(outmat, [rdim, cdim, sdim]), ...
    length(rnames), length(cnames), length(snames));
for kk = 1:length(snames)
    outtbl = array2table(outmat(:, :, kk), 'VariableNames', cnames);
    writetable(horzcat(rtbl, outtbl), filepath, 'Sheet', snames{kk});
end
end

function H = plot_shaded_err(x, y, errBar, lineProps)
patchSaturation = 0.3;

y=y(:)';
x=x(:)';
errBar = repmat(errBar(:)', 2, 1);

initialHoldStatus=ishold;
if ~initialHoldStatus, hold on,  end

H.mainLine=plot(x,y,lineProps{:});

mainLineColor=get(H.mainLine,'color');
patchColor=mainLineColor+(1-mainLineColor)*patchSaturation;

uE=y+errBar(1,:);
lE=y-errBar(2,:);

yP=[lE,fliplr(uE),lE(1)]';
xP=[x,fliplr(x),x(1)]';

H.edge=plot(xP,yP,'-');

set(H.edge, 'color',patchColor, ...
    'HandleVisibility','off', ...
    'Tag', 'shadedErrorBar_edge')

if strcmp(get(gca,'YAxisLocation'),'left')
    uistack(H.mainLine,'top')
end

if ~initialHoldStatus, hold off, end

end
