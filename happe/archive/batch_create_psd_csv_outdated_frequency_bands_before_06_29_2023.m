%% BATCH_CREATE_PSD_CSV
% Creates CSVs with average power in each frequency band from LEAPP data 
% processed via HAPPE and BEAPP.
% Author: Gerardo Parra, Last updated: 2023-01-25

%% TODO: Consider adding relative/log power
% clear workspace to avoid interference
clear

%% inputs - computer setting & data location
computer  = 'pc';    % 'mac' or 'pc'
rcfs_path = fullfile( 'Groups','LCN-Nelson-LEAP','Groups','BCH',        ...
                       'EEG','eeg_processed_data' );
data_path = rcfs_path;
subdir    = '';
id_prefix = alphanumericsPattern(2); % prefix to subject id numbers

%% INPUTS - age groups
% ages_to_run = {1,1,1,1,1,1,1,1,1,1,1,1};
ages_to_run = {0,0,0,0,0,0,0,0,0,1,0,0};
age_names   = {             ...
                '03_months',...
                '04_months',...
                '05_months',...
                '06_months',...
                '08_months',...
                '10_months',...
                '12_months',...
                '14_months',...
                '16_months',...
                '18_months',...
                '24_months',...
                '30_months',...
              };
age_paths = age_names;
% input one corresponding to each age group
beapp_tags = {
                '_2023-04-25',  ...
                '_2023-04-25',  ...
                '_2023-04-25',  ...
                '_2023-04-25',  ...
                '_2023-04-25',  ...
                '_2023-04-25',  ...
                '_2023-04-25',  ...
                '_2023-04-25',  ...
                '_2023-04-25',  ...
                '_2023-04-25',  ...
                '_2023-04-25',  ...
                '_2023-04-25'
             };

%% INPUTS - tasks
tasks_to_run = {1,0,0,0,0,0,0};
% run tags to load from
run_tags  = {
                '_BASELINE',        ...
                '_AUDITORY',        ...
                '_VEP',             ...
                '_DISENGAGEMENT',   ...
                 '_FLANKER',        ...
                '_WM_ISI',          ...
                '_WM_stim_duration' ...
            };
task_names = {'baseline','auditory','vep', ...
              'disengagement','flanker','wm_isi','wm_duration'};

%% INPUTS - frequency bands
bands_to_run = {0,1,1,1,1};
band_names   = {'delta','theta','alpha','beta','gamma'};
ranges       = {[0.5 3],[3 7],[7 10],[11 20],[20 40]};
% ranges       = {[0.5 4],[4 8],[8 12],[12 30],[30 47]};

%% INPUTS - rois
rois_to_run = {1,1,1,1,1};
roi_names   = {'all_1020','frontal','parietal','temporal','occipital'};
all_1020    = [ 22   9  33  36  58  62  96  70  83  45  52  92  11];
frontal     = [ 11  24 124];
parietal    = [ 52  62  92];
temporal    = [ 45 108];
occipital   = [ 70 83 75];
% frontal     = [  3   4  19  23  24 124];
% parietal    = [ 52  58  92  96  68];
chans       = {all_1020,frontal,parietal,temporal,occipital};

%% run happe+er scripts
% set up paths based on computer setting
switch computer
    case 'mac'
        path_prefix = '\\rc-fs.tch.harvard.edu\dmc-nelson'; 
    case 'pc'
        path_prefix = 'Z:';
end
data_path = fullfile(path_prefix,data_path);

% get structs for age groups, tasks, and freq bands to run
ages = struct( 'path',age_paths,'tag',beapp_tags,'name',age_names,'run',...
                ages_to_run );
tasks = struct('tag',run_tags,'name',task_names,'run',tasks_to_run);
bands = struct( 'name',band_names,'range',ranges,'run',bands_to_run,    ...
                'fvec_i',[] );
rois  = struct('name',roi_names,'chans',chans,'run',rois_to_run);
run_ages  = cellfind(ages_to_run);  run_tasks = cellfind(tasks_to_run);
run_bands = cellfind(bands_to_run); run_rois  = cellfind(rois_to_run);

% iterate through age groups
for i_a = run_ages
    age = ages(i_a);
    % prepare output folder and struct
    out_path = fullfile(data_path,age.path,subdir,'8 - PSD');
    if ~exist(out_path,'dir'), mkdir(out_path); end

    % iterate through tasks
    for i_t = run_tasks
        task = tasks(i_t);
        % initialize output struct
        PSD = struct('mean',[],'log',[]);
        methods = fieldnames(PSD); n_m = length(methods);
        % get files with matching tag and extension
        ftag  = ['*' task.tag '*.mat'];
        fpath = fullfile(age.path,subdir,'7 - BEAPP',['psd' age.tag],ftag);
        fpath = fullfile(data_path,fpath);
        files = dir(fpath); n_files = length(files);
        if ~n_files, disp(['No files found for ' age.name]), continue; end
        % iterate through task files
        for i_f = 1:n_files
            % load file data
            file = files(i_f);
            file.path = fullfile(file.folder,file.name);
            load(file.path,'file_proc_info','eeg_wfp','f')
            % add subject id to tables
            id = get_subj_id(file.name,3,id_prefix);
            for i_m = 1:n_m, PSD.(methods{i_m})(i_f).id = id; end

            % only need to do once: get freq band & chan indices
            if i_f == 1
                % get band indeces for fvec
                for i_b = run_bands
                    band = bands(i_b);
                    bands(i_b).fvec_i = get_indeces(f{1},band.range,1);
                end
                % get chan indeces for happe data (since rows from channels 
                % removed  during processing are removed from final data 
                % cells, unlike in beapp) 
                for i_r = run_rois
                    roi = rois(i_r);
                    chan_T = file_proc_info.chanlocs;
                    rois(i_r).chans_i = get_chan_indeces(roi.chans,chan_T);
                end
            end

            %% calculate power in each frequency band
            % iterate through conditions
            n_cond = size(eeg_wfp,2);
            for i_c = 1:n_cond
                trial_mean = mean(eeg_wfp{i_c},3,'omitnan');
                % iterate through rois
                for i_r = run_rois
                    roi = rois(i_r);
                    roi_mean = mean(trial_mean(roi.chans_i,:),1,'omitnan');
                    % iterate through bands
                    for i_b = run_bands
                        band = bands(i_b);
                        band_mean = mean( roi_mean(:,band.fvec_i),2,    ...
                                          'omitnan' );
                        T_col = [roi.name '_' band.name '_' int2str(i_c)];
                        PSD(i_c).mean(i_f).(T_col) = band_mean;
                        PSD(i_c).log(i_f).(T_col)  = log(band_mean);
                    end
                end
            end
        end

        % convert structs to tables and save
        for i_c = 1:n_cond
            for i_m = 1:n_m
                method = methods{i_m};
                % get table name
                if n_cond == 1, cond = ''; 
                else, cond = ['_' int2str(i_c)]; end
                out_fname = [ 'bch_' age.name '_' task.name cond      ...
                              '_psd_' method '_' get_date_str() '.csv' ];
                out_fname = fullfile(out_path,out_fname);
                % convert struct to table and save
                T = struct2table(PSD.(method));
                writetable(T,out_fname,'WriteVariableNames',true)
            end
        end
    end
end

%% helper functions
function i = cellfind(c)
    i = find(cell2mat(c));
end

function id = get_subj_id(fname,id_len,pfx)
    id_pat   = digitsPattern(id_len);
    nums     = extract(fname,id_pat);
    file_pfx = extractBefore(fname,id_pat);
    file_pfx = extract(file_pfx,pfx);
    if isempty(nums)
        error(['Subject id not found for: ' fname])
    elseif ~contains(file_pfx,pfx)
        disp(['Check filename for: ' fname])
    end
    id = [file_pfx{1} '_' nums{1}];
end