clear;
addpath('Z:\Toolboxes\Git\tfdecomp\matlab\');

%% run tfdecomp

cfg = [];

cfg.eeglab_path = 'Z:\Toolboxes\eeglab12_0_2_3b';
cfg.readdir = 'Z:\TemplateSwitch\EEG\processed\laplacian\';
cfg.writdir = 'Z:\TemplateSwitch\EEG\results\';
cfg.filename = '*laplacian.mat';
cfg.projectname = 'templateswitch_eegeye_singletrial';
cfg.prefix = 'ppxx';

cfg.seeds = {};
cfg.channels = 1:64;
cfg.frequencies = [2 40 25]; % from min to max in nsteps
cfg.cycles = [3 12]; % min max number of cycles used for min max frequency
cfg.scale = 'log'; % whether the above frequency settings will be logarithmically (recommended) or linearly scaled
cfg.times2save = -1000:25:500;
cfg.keeptrials = true;
%cfg.trialrejection = true;

cfg.srate = 512;
cfg.epochtime = [-2 2.5];

cfg.report_progress = true;
cfg.save_output = true;
cfg.overwrite = false;
cfg.plot_output.show = 'no';

tfdecomp(cfg);

%%
readdir = 'Z:\TemplateSwitch\EEG\results\';
cd(readdir);
filz=dir('*singletrial_tfdecomp.mat');
filx=dir('Z:\TemplateSwitch\EEG\processed\laplacian\*.mat');

%%
for subno=3:length(filz)
    load(filz(subno).name);
    load(['Z:\TemplateSwitch\EEG\processed\laplacian\' filx(subno).name]);
    convresult{1} = permute(cat(4,tf_pow{1},tf_pow{3}), [4 1 2 3]);    
    regressors{1} = zeros(size(convresult{1},1),7);
    trialnums = zeros(size(convresult{1},1),1);
    k=1;
    for triali = 1:ALLEEG(1).trials
        regressors{1}(k,:) = ALLEEG(1).epoch(triali).regressors;
        trialnums(k) = ALLEEG(1).epoch(triali).trialnum;
        k=k+1;
    end
    for triali = 1:ALLEEG(3).trials
        regressors{1}(k,:) = ALLEEG(3).epoch(triali).regressors;
        trialnums(k) = ALLEEG(3).epoch(triali).trialnum;
        k=k+1;
    end
    
    % sort back into original trial sequence
    [~,I]=sort(trialnums);
    regressors{1} = regressors{1}(I,:);
    convresult{1} = convresult{1}(I,:,:,:);
    
    
    
    convresult{2} = permute(cat(4,tf_pow{2},tf_pow{4}), [4 1 2 3]);
    regressors{2} = zeros(size(convresult{2},1),7);
    trialnums = zeros(size(convresult{2},1),1);
    k=1;
    for triali = 1:ALLEEG(2).trials
        regressors{2}(k,:) = ALLEEG(2).epoch(triali).regressors;
        trialnums(k) = ALLEEG(2).epoch(triali).trialnum;
        k=k+1;
    end
    for triali = 1:ALLEEG(4).trials
        regressors{2}(k,:) = ALLEEG(4).epoch(triali).regressors;
        trialnums(k) = ALLEEG(4).epoch(triali).trialnum;
        k=k+1;
    end
    
    % sort back into original trial sequence
    [~,I]=sort(trialnums);
    regressors{2} = regressors{2}(I,:);
    convresult{2} = convresult{2}(I,:,:,:);
    
    save(filz(subno).name,'convresult','regressors');
    
end
