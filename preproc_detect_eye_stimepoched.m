%% Preprocessing EEG I: steps prior to Independent Component Analysis
% Joram van Driel, VU Amsterdam, July 2016

%% setup preliminaries

% Paths
restoredefaultpath

%% Paths to packages

eeglab_path = 'Z:\Toolboxes\eeglab12_0_2_3b';
ft_path     = 'Z:\Toolboxes\fieldtrip-20150318';

path_list_cell = regexp(path,pathsep,'Split');

if (exist(ft_path,'dir') == 7) && ~any(ismember(ft_path,path_list_cell))
    addpath(ft_path,'-begin');
end

if (exist(eeglab_path,'dir') == 7) && ~any(ismember(eeglab_path,path_list_cell))
     addpath(eeglab_path);
     eeglab;
     
     path_list_cell = regexp(path,pathsep,'Split');
     % remove conflicting paths
     ft_on_eeglab = fullfile(eeglab_path,'external','fieldtrip-partial','utilities');
     if any(ismember(ft_on_eeglab,path_list_cell))
         rmpath(genpath(ft_on_eeglab));
     end
end
ft_defaults;

%%
% Data dirs
clear
eyeldir = 'Z:\TemplateSwitch\EYE\';
behvdir = 'Z:\TemplateSwitch\Behav\';
basedir = 'Z:\TemplateSwitch\EEG\';
readdir = [basedir 'raw\'];
writdir = [basedir 'processed\'];

% length of time for EEG epoch (in sec)
epochtime=[ -2 2.5 ];

cd(readdir)
sublist=dir('*.bdf');
sublist={sublist.name};

bhvlist=dir([behvdir '*.csv']);

%%
% #triggers:
% init = 100
% end_trigger = 200
% # practice
%     #forced
% PracForcSearchL = 21
% PracForcSearchR = 22
%     # free
% PracFreeSearchLL = 23
% PracFreeSearchLR = 24
% PracFreeSearchRL = 25
% PracFreeSearchRR = 26
% # no practice
%     #forced
% NPracForcSearchL = 31 # target left hemisphere
% NPracForcSearchR = 32 # target right hemisphere
%     #free
% NPracFreeSearchLL = 33 first target left hemisphere, second target left hemisphere
% NPracFreeSearchLR = 34
% NPracFreeSearchRL = 35
% NPracFreeSearchRR = 36
% # resposne triggers, target 1, target 2, a distractor, or nothing
% target1_fix = 1
% target2_fix = 2
% errorFix = 3
% notarget_fix = 4

%% Loop around subjects

for subno=1:length(sublist)
    
    %%
    if ~exist([writdir sublist{subno}(1:4)],'dir')
        mkdir([writdir sublist{subno}(1:4)])
    end
    
    % filenames assume a sensible basename of the raw .bdf file
    outfilename1=[ writdir sublist{subno}(1:4) filesep sublist{subno}(1:end-4) '_detectedEye.mat' ];
    outfilename2=[ writdir sublist{subno}(1:4) filesep sublist{subno}(1:end-4) '_epoched_rejectedTrials.mat' ];
    outfilename3=[ writdir sublist{subno}(1:4) filesep sublist{subno}(1:end-4) '_epoched_icaWeights.mat' ];
    outfilename4=[ writdir sublist{subno}(1:4) filesep sublist{subno}(1:end-4) '_epoched_cleaned.mat' ];

    % do all steps consecutively for one subject, or go to the step where
    % you left off
    
    reload = true;
    go_on = true;
    
    %%
    while go_on    
        if ~exist(outfilename1,'file')

            %% parse eye-tracking data to .mat file

            if ~exist([ eyeldir sublist{subno}(1:end-4) '.mat' ],'file')
                ET = parseeyelink([ eyeldir sublist{subno}(1:end-4) '.asc' ],[ eyeldir sublist{subno}(1:end-4) '.mat' ],'Sync_EEG_ET');
            end

            %% load data

            % read bdf file
            fprintf('Loading subject %i of %i\n',subno,length(sublist));
            EEG = pop_biosig(sublist{subno});

            %%
            origEEG = EEG; % backup

            % re-reference to linked mastoids/earlobes
            EEG = pop_reref(EEG,[find(strcmpi('EXG5',{EEG.chanlocs.labels})) find(strcmpi('EXG6',{EEG.chanlocs.labels}))],'refstate','averef');

            % re-reference EOG
            EEG.data(strcmpi({EEG.chanlocs.labels},'EXG1'),:) = squeeze(EEG.data(strcmpi({EEG.chanlocs.labels},'EXG1'),:)-EEG.data(strcmpi({EEG.chanlocs.labels},'EXG2'),:));
            EEG.chanlocs(strcmpi({EEG.chanlocs.labels},'EXG1')).labels = 'VEOG';

            EEG.data(strcmpi({EEG.chanlocs.labels},'EXG2'),:) = squeeze(EEG.data(strcmpi({EEG.chanlocs.labels},'EXG3'),:)-EEG.data(strcmpi({EEG.chanlocs.labels},'EXG4'),:));
            EEG.chanlocs(strcmpi({EEG.chanlocs.labels},'EXG2')).labels = 'HEOG';

            % remove unnecessary channels
            EEG = pop_select(EEG,'nochannel',67:length(EEG.chanlocs));

            % (high-pass) filter; eegfiltnew is quite fast
            EEG = pop_resample(EEG,512);
            EEG = pop_eegfiltnew(EEG,.5,0);

            EEG=pop_chanedit(EEG,'lookup','Z:\Toolboxes\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp');
            EEG.chanlocs(67:end)=[];

            %%
            EEGm = EEG;
            m=[EEGm.event.type];
            m=m-(64512);

            for i=1:length(m)
                EEGm.event(i).type=m(i);
            end

            mtable=tabulate(m);
            mtable(mtable(:,2)==0,:)=[];
            mtable=mtable(:,[1 2]);
            disp(mtable)
            keyboard
            
            %%
            EEG= EEGm;

            %% Eye-tracking synchronization

            EEG = eeg_checkset( EEG );

            % This function uses start/end triggers 100/200, imports only the 2nd
            % and 3d column of the eyelink data, corresponding to xgaze and ygaze,
            % and the booleans correspond to: import eye-events from eyelink (so
            % saccades and blinks are already imported as events); do regression to
            % check synchronization accuracy, do not filter eyetrack, and plot the
            % synchronization results

            EEG = pop_importeyetracker(EEG,[ eyeldir sublist{subno}(1:end-4) '.mat' ],[100 200] ,[2 3] ,{'gaze_x' 'gaze_y'},0,1,0,1);

            %-% Reject continuous data based on eye tracker

            % versionX needs argument rejectionmethod that can take value 1 or 2;
            % value 2 only marks the blinks with 'badeye' event code; this code is
            % then skipped during the saccade detectin routine

            EEG = rej_eyecontinX(EEG,[67 68] ,[1 1] ,[1680 1050] ,0,2);

            %-% Detect eye-movements and fixations with Eyetracker plugin

            % note: Eyelink may have also done this detection online; the option to
            % import these as events was set above;

            left_eye_xy     = [67 68]; % channels of left eye tracker
            right_eye_xy    = [];
            lambda_velocity = 7;       % for microsaccadas between 4 and 6
            mindur          = 4;       % minimum saccada duration in samples
            dva             = 0.0216;
            smooth          = 1;       % raw Eye data is smoothed to suppress noise
            globalthresh    = 1;       % one global threshold applied to all epochs (0: determined individually per epoch)
            clusterdist     = 26;      % when saccades within this value in samples are detected, they are treated as
            clustermode     = 2;       % option 4: one long saccade --> better to only take the first? comment by olaf %%%%%%%%%%%%%%%
            minvals         = [1 1];
            maxvals         = [1680 1050];

            % note: first set below to 1,0,0; this produces a figure of detection
            % quality, without storing the saccades/fixations as events; then if
            % satisfied, rerun with below set to 0,1,1; this skips the figure and
            % stores the saccades/fixations as events
            plotfig         = 1;
            writesac        = 0;
            writefix        = 0;

            EEG = detecteyemovementsX(EEG,left_eye_xy,right_eye_xy,lambda_velocity,mindur,dva,smooth,globalthresh,clusterdist,clustermode,minvals,maxvals,plotfig,writesac,writefix);
            keyboard
            
            plotfig         = 0;
            writesac        = 1;
            writefix        = 1;

            EEG = detecteyemovementsX(EEG,left_eye_xy,right_eye_xy,lambda_velocity,mindur,dva,smooth,globalthresh,clusterdist,clustermode,minvals,maxvals,plotfig,writesac,writefix);
           
            %%
            save(outfilename1,'EEG')
            reload = false;

        elseif ~exist(outfilename2,'file')
            
            if reload
                fprintf('Loading subject %i of %i for epoching and bad trial marking...\n',subno,length(sublist));
                load(outfilename1)
            else
                fprintf('Epoching and bad trial marking for subject %i of %i...\n',subno,length(sublist));
            end
            
            %% Epoch the data

            triggers={
                    'NPracForcSearchL' {'31'}; % target left hemisphere
                    'NPracForcSearchR' {'32'}; % target right hemisphere
                    'NPracFreeSearchLL' {'33'}; % first target left hemisphere, second target left hemisphere
                    'NPracFreeSearchLR' {'34'};
                    'NPracFreeSearchRL' {'35'};
                    'NPracFreeSearchRR' {'36'};
                };

            % now epoch on correct saccades after a search display
            EEG = pop_epoch( EEG,[triggers{:,2}], epochtime);
            EEG = applytochannels(EEG, [1:64] ,' pop_rmbase( EEG, []);');

            %% check for bad channels

            % divide data into data blocks of 10
            figure;
            nblocks = 10;
            ntrials=floor(EEG.trials/nblocks);
            for b=1:nblocks
                subplot(3,ceil(nblocks/3),b)
                newdata = reshape(squeeze(EEG.data(1:64,:,1+((b-1)*ntrials):b*ntrials)),64,size(EEG.data,2)*ntrials);
                zstd = std(zscore(newdata),[],2);
                topoplot(zstd,EEG.chanlocs(1:64),'electrodes','on','maplimits',[min(zstd) max(zstd)]);
                colorbar
            end
            colormap hot
            figure
            subplot(211)
            topoplot([],EEG.chanlocs(1:64),'style','empty','electrodes','labels');
            subplot(212)
            topoplot([],EEG.chanlocs(1:64),'style','empty','electrodes','numbers');
            keyboard

            %% custom-written artifact rejection based on Fieldtrip routines

            FT_EEG = eeglab2ft(EEG);

            cfg=[];
            cfg.channel = {EEG.chanlocs(1:64).labels}';
            dat = ft_selectdata(cfg,FT_EEG);
            dat = dat.trial;
            dat_filt = cell(size(dat));

            cfg=[];
            cfg.bpfreq = [110 140];
            cfg.bpfiltord = 6;
            cfg.bpfilttype = 'but';
            cfg.hilbert = 'yes';
            cfg.boxcar = 0.2;
            cfg.channels = 1:64;
            cfg.cutoff = 12;
            cfg.art_time = [-1 1.5]; % window within which artifacts are not allowed; note: also use this window for ICA!

            %-% loop over trials
            reverseStr='';
            numtrl=EEG.trials;
            tidx = dsearchn(FT_EEG.time{1}',cfg.art_time')';

            for ei=1:numtrl

                % display progress
                msg = sprintf('Filtering trial %i/%i...',  ei,numtrl);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));

                % filter in high-freq band
                tmpdat = ft_preproc_bandpassfilter(dat{ei},EEG.srate, cfg.bpfreq, cfg.bpfiltord, cfg.bpfilttype);
                tmpdat = ft_preproc_hilbert(tmpdat,cfg.hilbert);
                nsmp = round(cfg.boxcar*EEG.srate);
                if ~rem(nsmp,2)
                    % the kernel should have an odd number of samples
                    nsmp = nsmp+1;
                end
                tmpdat = ft_preproc_smooth(tmpdat, nsmp); % better edge behaviour
                dat_filt{ei} = double(tmpdat(:,tidx(1):tidx(2)));

                if ei==1
                    sumval = zeros(size(dat_filt{1},1), 1);
                    sumsqr = zeros(size(dat_filt{1},1), 1);
                    numsmp = zeros(size(dat_filt{1},1), 1);
                    numsgn = size(dat_filt{1},1);
                end

                % accumulate the sum and the sum-of-squares
                sumval = sumval + sum(dat_filt{ei},2);
                sumsqr = sumsqr + sum(dat_filt{ei}.^2,2);
                numsmp = numsmp + size(dat_filt{ei},2);
            end
            fprintf('\n')

            % avg and std
            datavg = sumval./numsmp;
            datstd = sqrt(sumsqr./numsmp - (sumval./numsmp).^2);

            zmax = cell(1, numtrl);
            zsum = cell(1, numtrl);
            zindx = cell(1, numtrl);

            indvec = ones(1,numtrl);
            for ei = 1:numtrl
                % initialize some matrices
                zmax{ei}  = -inf + zeros(1,size(dat_filt{ei},2));
                zsum{ei}  = zeros(1,size(dat_filt{ei},2));
                zindx{ei} = zeros(1,size(dat_filt{ei},2));

                nsmp          = size(dat_filt{ei},2);
                zdata         = (dat_filt{ei} - datavg(:,indvec(ei)*ones(1,nsmp)))./datstd(:,indvec(ei)*ones(1,nsmp));  % convert the filtered data to z-values
                zsum{ei}   = nansum(zdata,1);                   % accumulate the z-values over channels
                [zmax{ei},ind] = max(zdata,[],1);            % find the maximum z-value and remember it
                zindx{ei}      = cfg.channels(ind);                % also remember the channel number that has the largest z-value

                zsum{ei} = zsum{ei} ./ sqrt(numsgn);
            end
            cfg.threshold = median([zsum{:}]) + abs(min([zsum{:}])-median([zsum{:}])) + cfg.cutoff;

            figure
            plot([zsum{:}])
            hold on
            plot([1 length([zsum{:}])], [cfg.threshold cfg.threshold],'m')

            rejected_trials = zeros(1,numtrl);
            for ei=1:numtrl
                if sum(zsum{ei}>cfg.threshold)>0
                    rejected_trials(ei) = 1;
                end
            end

            EEG.reject.rejmanual = rejected_trials;
            EEG.reject.rejmanualE = zeros(EEG.nbchan,EEG.trials);
            pop_eegplot(EEG,1,1,0)
            keyboard % at this point, do following checks:
            % - check if rejected trials make sense, also check zsum and threshold figure
            % - turn back to topoplots of possible bad channels, and verify (add them to chans2interp txt file)

            if sum(EEG.reject.rejmanual)~=sum(rejected_trials)
                rejected_trials = EEG.reject.rejmanual;
            end

            save(outfilename2,'EEG','rejected_trials','cfg'); % the cfg variable contains the rejection settings, so these can be traced back (e.g. z-val cutoff)
            reload = false;
            
        elseif ~exist(outfilename3,'file')

            if reload
                fprintf('Loading subject %i of %i for ICA...\n',subno,length(sublist));
                load(outfilename2)
            else
                fprintf('Running ICA for subject %i of %i...\n',subno,length(sublist));
            end

            %% Independent Component Analysis

            % remove trials with artifacts detected in previous step
            EEG=pop_select(EEG,'notrial',find(rejected_trials));
            EEG=pop_select(EEG,'time', cfg.art_time);

            %% sets bad channels to zero if necessary

            fid=fopen([writdir 'chans2interp.txt'],'r');
            chans2interp={};
            while ~feof(fid)
                aline=regexp(fgetl(fid),'\t','split');
                chans2interp{size(chans2interp,1)+1,1}=aline{1};
                for i=2:length(aline)
                    chans2interp{size(chans2interp,1),i}=aline{i};
                end
            end

            chanind=1:64;
            subject_prefix = sublist{subno}(1:4);
            chans = chans2interp(strcmpi(chans2interp(:,1),subject_prefix),2:end);
            if ~isempty(chans{1})
                bad_chansidx = zeros(1,length(chans));
                for c=1:length(chans)
                    if ~isempty(chans{c}), bad_chansidx(c) = find(strcmpi({EEG.chanlocs.labels},chans(c))); end
                end
                bad_chansidx(bad_chansidx==0)=[];
                chanind(bad_chansidx)=[];
            end

            %% run ICA

            EEG=pop_runica(EEG,'icatype','runica','dataset',1,'options',{'extended',1},'chanind',chanind);
            ICAEEG.icawinv = EEG.icawinv;
            ICAEEG.icasphere = EEG.icasphere;
            ICAEEG.icaweights = EEG.icaweights;
            ICAEEG.icachansind = EEG.icachansind;

            save(outfilename3,'ICAEEG')
            reload = false;
            
        elseif ~exist(outfilename4,'file')
            
            if reload
                fprintf('Loading subject %i of %i for final cleaning...\n',subno,length(sublist));
                load(outfilename2)
                load(outfilename3)
            else
                fprintf('Final cleaning steps for subject %i of %i...\n',subno,length(sublist));
            end


            %% add search array position
            %- based on OpenSesame output
            %- also add trial number so we can keep track of trial selection in
            %- later stages

            %- get opensesame output data
            toggle=false;
            k=[];
            j=1;
            % colnames = {'previous_target','fixated_circle','df','s1_cat','s1_color','s2_cat','s2_color'};
            colnames = {'s1_x','s1_y','s2_x','s2_y','s3_x','s3_y','s4_x','s4_y','s5_x','s5_y','s6_x','s6_y'};
            searcharraypos=zeros(1600,6,2);
            fid = fopen([behvdir bhvlist(subno).name]);
            aline=regexp(fgetl(fid),',','split');
            for c=1:length(aline)
                if sum(strcmpi(aline(c),colnames))>0
                    k(j)=c;
                    j=j+1;
                end
            end
            colnames = aline(k);
            p=1;
            while ~feof(fid)
                aline=regexp(fgetl(fid),',','split');
                if str2double(aline(12))<1
                    continue
                end
                stringdat = aline(k);
                searcharraypos(p,1,1:2) = str2double(stringdat(1:2));
                searcharraypos(p,2,1:2) = str2double(stringdat(3:4));
                searcharraypos(p,3,1:2) = str2double(stringdat(5:6));
                searcharraypos(p,4,1:2) = str2double(stringdat(7:8));
                searcharraypos(p,5,1:2) = str2double(stringdat(9:10));
                searcharraypos(p,6,1:2) = str2double(stringdat(11:12));
                p=p+1;
            end
            fclose(fid);


            %% add regressors and array positions to epoch structure

            removetrial = zeros(1,EEG.trials);
            targetval = zeros(EEG.trials,2);

            for ei=1:EEG.trials-1

                try
                    tmp = find(cell2mat(EEG.epoch(ei).eventlatency)>0 & str2double(EEG.epoch(ei).eventtype)<3);
                    response = tmp(1);
                    cueonset = find(cell2mat(EEG.epoch(ei).eventlatency)==0); cueonset=cueonset(1);

                    if strcmp(EEG.epoch(ei).eventtype(response-1), 'fixation') && strcmp(EEG.epoch(ei).eventtype(response-2), 'saccade')
                        eyeevent = response-2;
                    elseif strcmp(EEG.epoch(ei).eventtype(response-1), 'saccade') && strcmp(EEG.epoch(ei).eventtype(response+1), 'fixation')
                        eyeevent = response-1;
                    else
                        removetrial(ei)=1;
                        continue
                    end

                    if str2double(EEG.epoch(ei).eventtype(cueonset))<33
                        targetval(ei,1)=1;
                        EEG.epoch(ei).condition = 'forced';
                    elseif str2double(EEG.epoch(ei).eventtype(cueonset))>32
                        targetval(ei,1)=2;
                        EEG.epoch(ei).condition = 'free';
                    end

                    targetval(ei,2) = str2double(EEG.epoch(ei).eventtype(response(1)));

                    % get saccade/fixation properties for GLM
                    tmpreg(1) = cell2mat(EEG.epoch(ei).eventlatency(eyeevent)); % when was the saccade with respect to search display onset
                    tmpreg(2) = cell2mat(EEG.epoch(ei).eventduration(eyeevent)); % how long did it take to fixate
                    tmpreg(3) = cell2mat(EEG.epoch(ei).eventsac_amplitude(eyeevent)); % what was the saccade amplitude
                    tmpreg(4) = cell2mat(EEG.epoch(ei).eventsac_vmax(eyeevent)); % what was the saccade velocity
                    tmpreg(5) = cell2mat(EEG.epoch(ei).eventsac_endpos_x(eyeevent)) - cell2mat(EEG.epoch(ei).eventsac_startpos_x(eyeevent)); % x of end of saccade relative to start (i.e. retinotopic displacement x-coordinate)
                    tmpreg(6) = cell2mat(EEG.epoch(ei).eventsac_endpos_y(eyeevent)) - cell2mat(EEG.epoch(ei).eventsac_startpos_y(eyeevent));

                    % get saccade trajectory from cue onset until correction
                    % fixation
                    saccpositions = zeros(10,2);
                    k=1;
                    toggle=true;
                    for ev=cueonset+1:length(EEG.epoch(ei).eventtype)
                        if strcmp(EEG.epoch(ei).eventtype(ev),'saccade')
                            if toggle
                                saccpositions(k,1)=cell2mat(EEG.epoch(ei).eventsac_startpos_x(ev));
                                saccpositions(k,2)=cell2mat(EEG.epoch(ei).eventsac_startpos_y(ev));
                                saccpositions(k+1,1)=cell2mat(EEG.epoch(ei).eventsac_endpos_x(ev));
                                saccpositions(k+1,2)=cell2mat(EEG.epoch(ei).eventsac_endpos_y(ev));
                                k=k+2;
                            end
                        elseif str2double(EEG.epoch(ei).eventtype(ev))<3
                            toggle=false;
                        end
                    end
                    saccpositions(saccpositions(:,1)==0,:)=[];
                    saccpositions(:,1)=saccpositions(:,1)-840;
                    saccpositions(:,2)=saccpositions(:,2)-525;

                    EEG.epoch(ei).regressors = tmpreg;
                    EEG.epoch(ei).searcharraypos = squeeze(searcharraypos(ei,:,:));
                    EEG.epoch(ei).saccade_trajectory = saccpositions;

                catch me, removetrial(ei)=1;
                end
            end

            %% check if this went OK

            ei=randi([1 EEG.trials],1,1);
            %ei=1210;
            disp(EEG.epoch(ei))
            figure
            hold on
            try
                plot(EEG.epoch(ei).searcharraypos(:,1),EEG.epoch(ei).searcharraypos(:,2),'o','markersize',20)
                if strcmp(EEG.epoch(ei).condition,'forced')
                    plot(EEG.epoch(ei).searcharraypos(1,1),EEG.epoch(ei).searcharraypos(1,2),'ro','markerfacecolor','r','markersize',20)
                elseif strcmp(EEG.epoch(ei).condition,'free')
                    plot(EEG.epoch(ei).searcharraypos(1:2,1),EEG.epoch(ei).searcharraypos(1:2,2),'ro','markerfacecolor','r','markersize',20)
                end
                plot(EEG.epoch(ei).saccade_trajectory(1,1),EEG.epoch(ei).saccade_trajectory(1,2),'-*k')
                plot(EEG.epoch(ei).saccade_trajectory(:,1),EEG.epoch(ei).saccade_trajectory(:,2),'-^','markersize',2)
            catch me
            end

            %% reject trials with irregular saccade trajectory

            badtrajectory = zeros(1,EEG.trials);
            for ei=1:EEG.trials

                try
                    sacc_angles = zeros(1,size(EEG.epoch(ei).saccade_trajectory,1)-1);
                    for si=1:size(EEG.epoch(ei).saccade_trajectory,1)-1
                        sacc_angles(si) = rad2deg(atan((EEG.epoch(ei).saccade_trajectory(si+1,2)-EEG.epoch(ei).saccade_trajectory(1,2))/(EEG.epoch(ei).saccade_trajectory(si+1,1)-EEG.epoch(ei).saccade_trajectory(1,1))));
                    end

                    targ_angle=[];
                    if strcmp(EEG.epoch(ei).condition,'forced')
                        targ_angle = rad2deg(atan((EEG.epoch(ei).searcharraypos(1,2)-EEG.epoch(ei).saccade_trajectory(1,2))/(EEG.epoch(ei).searcharraypos(1,1)-EEG.epoch(ei).saccade_trajectory(1,1))));
                    elseif strcmp(EEG.epoch(ei).condition,'free')
                        targ_angle(1) = rad2deg(atan((EEG.epoch(ei).searcharraypos(1,2)-EEG.epoch(ei).saccade_trajectory(1,2))/(EEG.epoch(ei).searcharraypos(1,1)-EEG.epoch(ei).saccade_trajectory(1,1))));
                        targ_angle(2) = rad2deg(atan((EEG.epoch(ei).searcharraypos(2,2)-EEG.epoch(ei).saccade_trajectory(1,2))/(EEG.epoch(ei).searcharraypos(2,1)-EEG.epoch(ei).saccade_trajectory(1,1))));
                    end

                    marktrial=zeros(1,length(targ_angle));
                    for ti=1:length(targ_angle)
                        if abs(sacc_angles-targ_angle(ti))>(45/2)
                            marktrial(ti)=1;
                        end
                    end
                    if sum(marktrial)>1
                        badtrajectory(ei)=1;
                    end
                catch me
                end
            end


            %% get trial position along streaks of stay/switch trials

            streaks = zeros(EEG.trials,2);   

            %- do some juggling
            targetval(targetval(:,2)==2,2)=targetval(targetval(:,2)==2,2)-2;
            targetval(:,3)=targetval(:,2);
            targetval(2:end,2)=targetval(2:end,2)-targetval(1:end-1,2);
            targetval(1,2)=0; targetval(end,2)=0;
            targetval(:,2)=abs(targetval(:,2));
            targetval(targetval(:,2)==0,2) = targetval(targetval(:,2)==0,2)-1;

            % get streaks of stay strials and switches and determine position of
            % trial with respect to previous and next switch

            switchi=0; stayi=0;
            for ei=length(targetval):-1:2
                if targetval(ei,1)==targetval(ei-1,1)

                    if targetval(ei,2)==1
                        switchi=switchi+1;
                        stayi=0;
                    elseif targetval(ei,2)==-1
                        switchi=0;
                        stayi=stayi+1;
                    end
                    streaks(ei,1:2) = [switchi stayi];
                else
                    streaks(ei,:)=0;
                end
            end

            streaks(:,2)=streaks(:,2).*-1;
            streaks(find(removetrial),:)=0;

            for ei=1:EEG.trials
                EEG.epoch(ei).regressors = [EEG.epoch(ei).regressors sum(streaks(ei,:))];
                EEG.epoch(ei).trialnum = ei;
            end


            %% remove artifact and error trials

            EEG=pop_select(EEG,'notrial',find(rejected_trials+removetrial+badtrajectory));

            %% add ICA weights to EEG structure

            EEG.icachansind = ICAEEG.icachansind;
            EEG.icasphere = ICAEEG.icasphere;
            EEG.icaweights = ICAEEG.icaweights;
            EEG.icawinv = ICAEEG.icawinv;
            EEG = eeg_checkset(EEG);
            clear ICAEEG

            %% Temporarily remove bad channels

            fid=fopen([writdir 'chans2interp.txt'],'r');
            chans2interp={};
            while ~feof(fid)
                aline=regexp(fgetl(fid),'\t','split');
                chans2interp{size(chans2interp,1)+1,1}=aline{1};
                for i=2:length(aline)
                    chans2interp{size(chans2interp,1),i}=aline{i};
                end
            end

            chanind=1:64;
            subject_prefix = sublist{subno}(1:4);
            chans = chans2interp(strcmpi(chans2interp(:,1),subject_prefix),2:end);
            if ~isempty(chans{1})
                bad_chansidx = zeros(1,length(chans));
                for c=1:length(chans)
                    if ~isempty(chans{c}), bad_chansidx(c) = find(strcmpi({EEG.chanlocs.labels},chans(c))); end
                end
                bad_chansidx(bad_chansidx==0)=[];
                chanind(bad_chansidx)=[];
            else
                bad_chansidx = 0;
            end

            %%
            if sum(bad_chansidx)>0

                EEG2 = pop_select(EEG,'nochannel',[bad_chansidx 65 66]);
            else
                
                EEG2 = pop_select(EEG,'nochannel',[65 66]);
                %% detect oculomotor ICs informed by eye-tracker
            end
            
            [EEG2,~,bad] = pop_eyetrackerica(EEG2,'saccade','fixation',[5 5] ,1.1,3,1,1);
            prompt = 'Do you agree with selected components? y/n --> ';
            agree = input(prompt,'s');
            if strcmp(agree,'y')
                comps2remove = bad;
            elseif strcmp(agree,'n')
                prompt = 'Please provide alternative component selection: --> ';
                newbad = input(prompt);
                comps2remove = newbad;
            end
            EEG2 = pop_subcomp( EEG2, comps2remove, 0);
            
            %% put IC-cleaned channels back in original EEG structure and interpolate bad ones
            
            good_chans = 1:EEG.nbchan;
            if sum(bad_chansidx)>0
                good_chans([bad_chansidx 65 66])=[];
                EEG.data(good_chans,:,:) = EEG2.data;
                EEG = eeg_interp(EEG,bad_chansidx);
            else
                good_chans([65 66])=[];
                EEG.data(good_chans,:,:) = EEG2.data;
            end
            
            clear EEG2

            save(outfilename4,'EEG','comps2remove','bad_chansidx');
            go_on = false;
        else
            go_on = false;
        end
    end
end


