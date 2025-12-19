clear
close all
clc

addpath('functions')
addpath('fieldtrip-20211020'); % path to your own fieldtrip toolbox
%%
srate = 1200; % sampling rate, Hz
J = customcolormap_preset('red-white-blue');
target = 'stim-wostim';
on_save = 'true';
mkdir(sprintf('figures/%s',target));
save_dir = sprintf('figures/%s',target);
%% load data
onstim_1 = load('../patient_result/result_01102024_2.mat');
onstim_2 = load('../patient_result/result_01112024_1.mat');
offstim_1 = load('../patient_result/result_01052024_1.mat');
offstim_2 = load('../patient_result/result_01052024_2.mat');
offstim_3 = load('../patient_result/result_01102024_1.mat');
on_1 = load('../matfiles/01102024_2.mat');
on_2 = load('../matfiles/01112024_1.mat');
off_merge = load('../matfiles//merge_off.mat');
%% barplot for behavioral data
offstim_trial = [offstim_1.trial_img; offstim_2.trial_img; offstim_3.trial_img];
offstim_rating = [offstim_1.trial_rating; offstim_2.trial_rating; offstim_3.trial_rating];
tmp1 = find(offstim_trial==2);
bar_off = offstim_rating(tmp1);
bar_off(16:18) = [];

onstim_trial = [onstim_1.trial_img; onstim_2.trial_img];
onstim_rating = [onstim_1.trial_rating; onstim_2.trial_rating];
tmp1 = find(onstim_trial==2);
tmp1(tmp1==8) = [];
bar_on = onstim_rating(tmp1);

% bar graph

arr_1 = [mean(bar_off)];
arr_2 =  [mean(bar_on)];

arr_sem_1 = [sem(bar_off)];
arr_sem_2 = [sem(bar_on)];
figure('Renderer', 'painters', 'Position', [100 100 300 300]);

bar(1,arr_1,'r');
hold on;
bar(2,arr_2,'b');
errorbar(1,arr_1,arr_sem_1,'k');
errorbar(2,arr_2,arr_sem_2,'k');

legend({sprintf('w/o stim (n=%d)',size(bar_off,1)),sprintf('with stim (n=%d)',size(bar_on,1))},'Location','Southwest');
ylabel('Desire to use drugs (Craving)');
[h,p] = ttest2(bar_off,bar_on,'tail','both')
ylim([1 8]);
xticks([1 2]);  xticklabels({'',''});

%%
pow_on = [on_1.freq_avg_bsl_w.powspctrm; on_2.freq_avg_bsl_w.powspctrm];
pow_off = [off_merge.freq_avg_bsl_w.powspctrm ];

%% Line plot
for chn = 1 :5
    freq = on_1.freq_avg_bsl_w.freq;
    tmp = round(on_1.freq_avg_bsl_w.time,2);
    a_on = find(tmp==-1);
    b_on = find(tmp==6);

    tmp = round(off_merge.freq_avg_bsl_w.time,2);
    a_off = find(tmp==-1);
    b_off = find(tmp==6);

    time = on_1.freq_avg_bsl_w.time(1,a_on:b_on);

    data_on=squeeze(pow_on(:,chn,:,a_on:b_on)); %Trials x 1 x Freq x Time
    data_off=squeeze(pow_off(:,chn,:,a_off:b_off)); %Trials x 1 x Freq x Time

    for ff= 1:7
        figure('Renderer', 'painters', 'Position', [100 100 500 230]);

        if ff == 1; findx = [1 : 12];
        elseif ff == 2; findx =  [13 : 17];
        elseif ff == 3; findx =  [18 : 23];
        elseif ff == 4; findx =  [24 : 40];
        elseif ff == 5; findx =  [41: 58];
        elseif ff == 6; findx =  [62 :80];
        elseif ff == 7; findx =  [1: 17];
        end

        %Plot
        logplot=data_off(:,findx,:);
        STE=squeeze(nanstd(logplot,[],1)/sqrt((size(logplot,1))));
        ck_shadedErrorBar(time,(nanmean(nanmean(logplot,2),1)),nanmean(STE,1),'r',1);
        hold on
        logplot=data_on(:,findx,:);
        STE=squeeze(nanstd(logplot,[],1)/sqrt((size(logplot,1))));
        ck_shadedErrorBar(time,(nanmean(nanmean(logplot,2),1)),nanmean(STE,1),'b',1);
        vline(0,'k--');

        % cluster permutation
        data1= squeeze(nanmean(data_on(:,findx,:),2));
        data2= squeeze(nanmean(data_off(:,findx,:),2));
        permu = cluster_permutation_main_1d(data1,data2);

        psig=find(permu~=0);
        psig
        for pp=1:length(psig)
            h=line([time(psig(pp))-.3 time(psig(pp))+.3],[0 0]);set(h,'LineWidth',2);set(h,'Color','k');
        end

        title(['Average Power ' num2str(round(freq(findx(1)))) '-' num2str(round(freq(findx(end)))) ' Hz']);
        legend('',sprintf('w/o stim (n=%d)',size(data_off,1)),'',sprintf('with stim (n=%d)',size(data_on,1)),'Location','Northwest');
        xlabel('Time (s)');
        xlim([-1 6]);
        ylabel('Baseline corrected power (dB)');

        if on_save
            if isempty(psig)
                saveas(gcf,[sprintf('figures\\%s\\',target) sprintf('line_%s_%d.pdf',on_1.freq_avg_bsl.label{chn},ff)]);
                close
            else
                saveas(gcf,[sprintf('figures\\%s\\',target) sprintf('line_%s_%d_sig.pdf',on_1.freq_avg_bsl.label{chn},ff)]);
                close
            end
        end

        % bar graph
        tmp = round(on_1.freq_avg_bsl_w.time,2);
        a_on = find(tmp==4.5);
        b_on = find(tmp==6);

        tmp = round(off_merge.freq_avg_bsl_w.time,2);
        a_off = find(tmp==4.5);
        b_off = find(tmp==6);

        bar_off = squeeze(mean(mean(pow_off(:,chn,findx,a_off:b_off),3),4));
        bar_on = squeeze(mean(mean(pow_on(:,chn,findx,a_on:b_on),3),4));

        arr_1 = [mean(bar_off)];
        arr_2 =  [mean(bar_on)];

        arr_sem_1 = [sem(bar_off)];
        arr_sem_2 = [sem(bar_on)];
        figure('Renderer', 'painters', 'Position', [100 100 300 300]);

        bar(1,arr_1,'r');
        hold on;
        bar(2,arr_2,'b');
        errorbar(1,arr_1,arr_sem_1,'k');
        errorbar(2,arr_2,arr_sem_2,'k');

        legend({sprintf('w/o stim (n=%d)',size(bar_off,1)),sprintf('with stim (n=%d)',size(bar_on,1))},'Location','Northeast');
        ylabel('baseline corrected power (dB)');
        [h,p] = ttest2(bar_off,bar_on,'tail','right')
        title([' Average Power ' num2str(round(freq(findx(1)))) '-' num2str(round(freq(findx(end)))) ' Hz']);
        xticks([1 2]);  xticklabels({'',''});
        if on_save
            if p>=0.05
                saveas(gcf,[sprintf('figures\\%s\\',target) sprintf('bar_%s_%d_%.3f.pdf',on_1.freq_avg_bsl.label{chn},ff,p)]);
                close
            else
                saveas(gcf,[sprintf('figures\\%s\\',target) sprintf('bar_%s_%d_%.3f_sig.pdf',on_1.freq_avg_bsl.label{chn},ff,p)]);
                close
            end
        end

    end
end


%% FFT
cfg = [];
cfg.trials = on_1.watch_trial;
proc_ftdata_on_1 = ft_selectdata(cfg,on_1.proc_ftdata);
cfg.trials = on_2.watch_trial;
proc_ftdata_on_2 = ft_selectdata(cfg,on_2.proc_ftdata);
cfg = [];
proc_ftdata_on = ft_appenddata(cfg,proc_ftdata_on_1, proc_ftdata_on_2);

cfg = [];
cfg.trials = off_1.watch_trial;
proc_ftdata_off_1 = ft_selectdata(cfg,off_1.proc_ftdata);
cfg.trials = off_2.watch_trial;
proc_ftdata_off_2 = ft_selectdata(cfg,off_2.proc_ftdata);
cfg.trials = off_3.watch_trial(1:3);
proc_ftdata_off_3 = ft_selectdata(cfg,off_3.proc_ftdata);
cfg = [];
proc_ftdata_off = ft_appenddata(cfg,proc_ftdata_off_1, proc_ftdata_off_2,proc_ftdata_off_3);


end_idx = dsearchn(proc_ftdata_on.time{1}',[2;3]);
for ii = 1 : length(proc_ftdata_on.trial)
    fft_ftdata_on.trial{1,ii} = proc_ftdata_on.trial{1,ii}(:,end_idx(1):end_idx(2));
    fft_ftdata_on.time{1,ii} = proc_ftdata_on.time{1,ii}(:,end_idx(1):end_idx(2));
end
fft_ftdata_on.fsample =proc_ftdata_on.fsample;
fft_ftdata_on.sampleinfo = proc_ftdata_on.sampleinfo;
fft_ftdata_on.label = proc_ftdata_on.label;

for ii = 1 : length(proc_ftdata_off.trial)
    fft_ftdata_off.trial{1,ii} = proc_ftdata_off.trial{1,ii}(:,end_idx(1):end_idx(2));
    fft_ftdata_off.time{1,ii} = proc_ftdata_off.time{1,ii}(:,end_idx(1):end_idx(2));
end
fft_ftdata_off.fsample =proc_ftdata_off.fsample;
fft_ftdata_off.sampleinfo = proc_ftdata_off.sampleinfo;
fft_ftdata_off.label = proc_ftdata_off.label;

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.foi = [1:0.25:120]; 
cfg.tapsmofrq = 1; 
cfg.output = 'pow';
cfg.keeptrials  = 'yes';
cfg.pad= 'nextpow2';
freq = cfg.foi;
TFRS_on = ft_freqanalysis(cfg, fft_ftdata_on);
TFRS_off = ft_freqanalysis(cfg, fft_ftdata_off);

chn = 1;
figure('Renderer', 'painters', 'Position', [100 100 350 300],'visible','on');
logplot1=squeeze(TFRS_off.powspctrm(:,chn,:));
STE=squeeze(nanstd(logplot1,[],1)/sqrt((size(logplot1,1))));
ck_shadedErrorBar(TFRS_off.freq,nanmean(logplot1,1),nanmean(STE,1),'r',1); hold on;


logplot2=squeeze(TFRS_on.powspctrm(:,chn,:));
STE=squeeze(nanstd(logplot2,[],1)/sqrt((size(logplot2,1))));
ck_shadedErrorBar(TFRS_on.freq,nanmean(logplot2,1),nanmean(STE,1),'k',1);
legend('','w/o stim','','with stim');
xlabel('Frequency (Hz)');
ylabel('Power (\muV)')
xlim([1 55]);
set(gca,'YScale','log')
set(gca,'XScale','log')
grid on;

clear p psig h
for i = 1 : length(freq)
    p(i) = perm_mdiff(logplot1(:,i)',logplot2(:,i)',1000);
end
psig = find(p<0.05);
for pp=1:length(psig)
    h=line([freq(psig(pp))-.25 freq(psig(pp))+.25],[5 5]);set(h,'LineWidth',2);set(h,'Color','k');
end

legend('','w/o stim','','with stim');
