clear, close all
restoredefaultpath
addpath(genpath('Z:\Toolboxes\eeglab12_0_2_3b'));
addpath(genpath('Z:\Toolboxes\colormaps'));
addpath('Z:\Toolboxes\Statstuff\')

% IO definitions
readdir = 'Z:\TemplateSwitch\EEG\results\prestim_baseline\';
cd(readdir)
filz = dir('*eegeye_tfdecomp.mat');

%% load data

load(filz(1).name)
tf_all_pow  = zeros([length(filz) size(tf_pow)]);
tf_all_pow(1,:,:,:,:)   = tf_pow;

tf_all_sync  = zeros([length(filz) size(tf_sync)]);
tf_all_sync(1,:,:,:,:,:)   = tf_sync;

nsubjects = length(filz);

for subno=2:nsubjects
    load(filz(subno).name,'tf_pow','tf_sync','dim')
    tf_all_pow(subno,:,:,:,:)    = tf_pow;
    tf_all_sync(subno,:,:,:,:,:) = tf_sync;
end

%%

elecs2plot_lab = {'o1','oz','o2'}; 
elecs2plot_idx = zeros(1,length(elecs2plot_lab));
for chani=1:length(elecs2plot_lab)
    elecs2plot_idx(chani) = find(strcmpi({dim.chans.labels}, elecs2plot_lab{chani}));
end

figure
for sp=1:4
    subplot(4,2,sp)
    contourf(dim.times,dim.freqs,squeeze(mean(mean(tf_all_pow(:,sp,elecs2plot_idx,:,:),1),3)),50,'linecolor','none')
    set(gca,'yscale','log','ytick',[2 4 8 16 32],'clim',[-3 3])
end

elecs2plot_lab = {'fcz'}; 
elecs2plot_idx = zeros(1,length(elecs2plot_lab));
for chani=1:length(elecs2plot_lab)
    elecs2plot_idx(chani) = find(strcmpi({dim.chans.labels}, elecs2plot_lab{chani}));
end

for sp=1:4
    subplot(4,2,sp+4)
    contourf(dim.times,dim.freqs,squeeze(mean(mean(tf_all_pow(:,sp,elecs2plot_idx,:,:),1),3)),50,'linecolor','none')
    set(gca,'yscale','log','ytick',[2 4 8 16 32],'clim',[-1 1])
end

%%
elecs2plot_lab = {'po3','po4'};
elecs2plot_idx = zeros(1,length(elecs2plot_lab));
for chani=1:length(elecs2plot_lab)
    elecs2plot_idx(chani) = find(strcmpi({dim.chans.labels}, elecs2plot_lab{chani}));
end

tfwin = [-500 -250; 14 20];
time=dsearchn(dim.times',tfwin(1,:)')'; time=time(1):time(2);
freq=dsearchn(dim.freqs',tfwin(2,:)')'; freq=freq(1):freq(2);

figure
subplot(211)
contourf(dim.times,dim.freqs,squeeze(mean(mean(mean(tf_all_pow(:,:,elecs2plot_idx,:,:),1),2),3)),50,'linecolor','none')
set(gca,'yscale','log','ytick',[2 4 8 16 32],'clim',[-1 1],'xlim',[-500 500])
rectangle('position',[tfwin(1,1) tfwin(2,1) tfwin(1,2)-tfwin(1,1)  tfwin(2,2)-tfwin(2,1)])

subplot(212)
topoplot(squeeze(mean(mean(mean(mean( tf_all_pow(:,:,:,freq,time), 1),2),4),5)),dim.chans,'electrodes','off');


