clear;

%% Laplacian

addpath(genpath('Z:\Toolboxes\eeglab12_0_2_3b\functions'))
addpath('Z:\Toolboxes\eeglab12_0_2_3b\customized_functions')
addpath('Z:\Toolboxes\Git\tfdecomp\matlab\')

% addpath('/Volumes/Repository2/CogPsy/J.van.Driel/Toolboxes/Git/tfdecomp/matlab')

% files to use in analysis
readdir = 'Z:\TemplateSwitch\EEG\processed\';
cd(readdir)
subs = dir('pp*');

%%
for subno=1:length(subs)
    
    cd([readdir subs(subno).name])
    filz = dir('*cleaned.mat');
    if isempty(filz), continue; end
    
    outputname = [readdir 'laplacian' filesep filz(1).name(1:4) '_templateswitch_eegeye_laplacian_stimlocked.mat'];
    
    if exist(outputname,'file'), continue; end
    
    fprintf('Preparing subject %i/%i for level-1 analysis\n',subno,length(subs))
    
    %% load data
        
    load(filz(1).name,'EEG');
    
    %% resort data to be saccade-locked --> done in tfdecomp_eegeye
    
%     [startpoint,baselinepoint,removetrial] = deal(zeros(length(EEG.epoch),1));
%     
%     for ei=1:EEG.trials-1
%         
%         try
%             cueonset = find(cell2mat(EEG.epoch(ei).eventlatency)==0); cueonset=cueonset(1);
%             % constraint of search disp onset always followed by a saccade
%             % sometimes there was a saccade just before onset followed by
%             % correct fixation, but these are often blinks or accidental
%             % saccades to correct targets
%             if strcmp(EEG.epoch(ei).eventtype(cueonset+1), 'saccade')
%                 removetrial(ei)=0;
%             else
%                 removetrial(ei)=1;
%                 continue
%             end
%             
%             eyeevent = cueonset+1;
% 
%             % find saccade position in time
%             [~,saccloc] = min(abs(EEG.times - EEG.epoch(ei).eventlatency{eyeevent})); % first event prior to fixated target is fixation onset (i.e. saccade offset); second before target trigger is saccade onset
%             [~,zeroloc] = min(abs(EEG.times - 0 ));
%             
%             % find start point and replace with saccade-aligned data
%             startpoint(ei)    = saccloc-zeroloc+1;
%             EEG.data(:,1:end-startpoint(ei)+1,ei) = EEG.data(:,startpoint(ei):end,ei);
%             
%         catch me, removetrial(ei)=1; % should be trials where there was a saccade during search display onset (stim onset is followed by a fixation event)
%         end
%     end
%     
%     EEG = pop_select(EEG,'notrial',find(removetrial));
    
    %% code stay/switch trials
    
    for ei=2:EEG.trials-1
        if EEG.epoch(ei).regressors(end)>0
            EEG.epoch(ei).trialtype = [EEG.epoch(ei).condition '_switch'];
        elseif EEG.epoch(ei).regressors(end)<0
            EEG.epoch(ei).trialtype = [EEG.epoch(ei).condition '_stay'];
        end
    end
        
    %% divide into conditions
    
    connames = {'forced_stay','free_stay','forced_switch','free_switch'};
    ALLEEG = EEG;
    for condi=1:length(connames)
        ALLEEG(condi) = EEG;
        ALLEEG(condi).data = ALLEEG(condi).data(:,:,(strcmpi({EEG.epoch.trialtype}, connames{condi})));
        ALLEEG(condi).epoch = ALLEEG(condi).epoch((strcmpi({EEG.epoch.trialtype}, connames{condi})));
        ALLEEG(condi).trials=size(ALLEEG(condi).data,3);
        ALLEEG(condi).setname=connames{condi};
    end
    
    %% apply Laplacian
    chanx=[ALLEEG(1).chanlocs(1:64).X];
    chany=[ALLEEG(1).chanlocs(1:64).Y];
    chanz=[ALLEEG(1).chanlocs(1:64).Z];
    
    for condi=1:length(ALLEEG)
        fprintf('Computing Laplacian operator for condition %i...\n', condi);
        [ ALLEEG(condi).data(1:64,:,:), ALLEEG(condi).G, ALLEEG(condi).H ] = laplacian_perrinX(ALLEEG(condi).data(1:64,:,:),chanx,chany,chanz);
    end

    %% save
    save(outputname,'ALLEEG');
    
end


%% run tfdecomp

cfg = [];

cfg.eeglab_path = 'Z:\Toolboxes\eeglab12_0_2_3b';
cfg.readdir = 'Z:\TemplateSwitch\EEG\processed\laplacian\';
cfg.writdir = 'Z:\TemplateSwitch\EEG\results\prestim_baseline\';

% cfg.eeglab_path = '/Volumes/Repository2/CogPsy/J.van.Driel/Toolboxes/eeglab12_0_2_3b';
% cfg.readdir = '/Volumes/Repository2/CogPsy/J.van.Driel/TemplateSwitch/EEG/processed/laplacian/';
% cfg.writdir = '/Volumes/Repository2/CogPsy/J.van.Driel/TemplateSwitch/EEG/results/prestim_baseline/';

cfg.filename = '*laplacian_stimlocked.mat';
cfg.projectname = 'templateswitch_eegeye';
cfg.relocking = true; % relocking to saccade
cfg.prefix = 'ppxx';

cfg.seeds = {'f5','f6','fcz','cz','oz','o1','o2','poz'};
cfg.channels = 1:64;
cfg.connectivity = 'pli'; % 'pli','iscp','both','none'
cfg.frequencies = [2 40 25]; % from min to max in nsteps
cfg.cycles = [3 12]; % min max number of cycles used for min max frequency
cfg.scale = 'log'; % whether the above frequency settings will be logarithmically (recommended) or linearly scaled
cfg.times2save = -1000:25:500;
cfg.basetime = [-400 -100]; % so whole epoch
cfg.baselinetype = 'conavg'; % 'conavg' or 'conspec'
cfg.erpsubtract = false; % if true, non-phase-locked (i.e. "induced") power will be computed by subtracting the ERP from each single trial
cfg.matchtrialn = true; % if true, conditions will be equated in terms of trial count (so SNR is comparable across conditions)

cfg.srate = 512;
cfg.epochtime = [-2 2.5];

cfg.report_progress = true;
cfg.save_output = true;
cfg.overwrite = false;
cfg.plot_output.show = 'no';
cfg.plot_output.chan = {'fcz'};
cfg.plot_output.freq = [4 8];
cfg.plot_output.time = [-500 -200];
cfg.plot_output.connames = {'forced_stay','free_stay','forced_switch','free_switch'};
cfg.plot_output.save = true;

[tf_pow, ~, tf_sync, dim] = tfdecomp_eegeye(cfg);

%% run tfmultiplot

cfg = [];

cfg.metric = {'pow'};
cfg.chan = {'po8','po4','o2','po7','po3','o1'};
cfg.freq = [8 14];
cfg.time = [-500 0];
cfg.scale = 'log';
cfg.connames = {'free switch > free stay'};
cfg.markevents = [0];
cfg.concomp = [4;2];

tfmultiplot(cfg,tf_pow,dim);

%%
cfg = [];

cfg.metric = {'pow'};
% cfg.chan = {'f5','f6'};
cfg.chan = {'po8','po4','o2','po7','po3','o1'};
cfg.freq = [8 14];
cfg.time = [-500 -100];
cfg.scale = 'log';
cfg.connames = {'forced stay','free stay','forced switch','free switch'};
cfg.markevents = [0];
%cfg.concomp = [3 4; 1 2];

tfmultiplot(cfg,tf_pow,dim);

%%
cfg = [];

cfg.metric = {'pow'};
cfg.chan = {'f5','f6'};
cfg.freq = [8 14];
cfg.time = [-900 -300];
cfg.scale = 'log';
cfg.connames = {'switch>stay'};
cfg.markevents = [0];
cfg.concomp = [3 4; 1 2];

tfmultiplot(cfg,tf_pow,dim);

%%
cfg = [];

cfg.metric = 'pli';
cfg.seed = {'fcz'};
cfg.chan = {'p1','p2','pz'};
cfg.freq = [8 14];
cfg.time = [-800 -400];
cfg.scale = 'log';
cfg.connames = {'switch>stay'};
cfg.concomp = [3 4; 1 2];
cfg.markevents = 0;

tfmultiplot(cfg,tf_sync,dim);

