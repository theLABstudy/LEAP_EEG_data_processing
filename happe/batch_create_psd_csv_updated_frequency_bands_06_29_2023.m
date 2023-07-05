%% BATCH_CREATE_PSD_CSV
% Creates CSVs with average power in each frequency band from LEAPP data 
% processed via HAPPE and BEAPP.
% Author: Gerardo Parra, Last updated: 2023-01-25

%% TODO: Consider adding relative/log power
% clear workspace to avoid interference
clear

%% inputs - computer setting & data location
computer  = 'mac';    % 'mac' or 'pc'
rcfs_path = fullfile( 'Groups','LCN-Nelson-LEAP','Groups','BCH',        ...
                       'EEG','eeg_processed_data' );
% rcfs_path = fullfile( 'Groups','LCN-Nelson-LEAP','Groups','BCH',        ...
%                        'EEG','eeg_processed_data','mff_for_processing' );
data_path = rcfs_path;
subdir    = '';
id_prefix = alphanumericsPattern(2); % prefix to subject id numbers

%% INPUTS - age groups
ages_to_run = {1,1,1,1,1,1,1,1,1,1,1,1};
% ages_to_run = {0,1,1,0,1,1,0,1,1,1,1,1};
% ages_to_run = {0,0,0,0,0,0,0,0,0,0,0,0};
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
                '36_months',...
              };
age_paths = age_names;
% input one corresponding to each age group
beapp_tags = {
                '_2023-05-25',  ...
                '_2023-06-15',  ...
                '_2023-06-15',  ...
                '_2023-05-25',  ...
                '_2023-06-15',  ...
                '_2023-06-15',  ...
                '_2023-06-15',  ...
                '_2023-06-15',  ...
                '_2023-06-15',  ...
                '_2023-06-15',  ...
                '_2023-06-15',  ...
                '_2023-06-15',  ...
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
bands_to_run = {1,1,1,1,1,1,1};
band_names   = {'delta','theta','low_alpha','high_alpha','beta','low_gamma','high_gamma'}; %updated 6/29
ranges       = {[2 3.99],[4 5.99],[6 8.99],[9 11.99],[12 30],[30 45],[65 90]}; %updated 6/29

%% INPUTS - rois
rois_to_run = {1,1,1,1,1,1,1,1,1,1,1};
roi_names = {'occipital','L_Temporal','R_Temporal','temporal', ...
             'L_Parietal','R_Parietal','parietal','L_lat_frontal', ...
             'R_lat_frontal','Frontal_central','frontal'};
occipital       = [70, 75, 83, 71, 76];
L_Temporal      = [45, 40, 46, 41, 36];
R_Temporal      = [108, 109, 102, 103, 104];
temporal        = [45, 40, 46, 41, 36, 108, 109, 102, 103, 104];
L_Parietal      = [52, 53, 59, 60];
R_Parietal      = [92, 86, 85, 91];
parietal        = [52, 53, 59, 60, 92, 86, 85, 91];
L_lat_frontal   = [33, 27, 34, 28, 24];
R_lat_frontal   = [122, 123, 116, 117,124];
Frontal_central = [6, 13, 112, 12, 5];
frontal         = [33, 27, 34, 28, 24, 122, 123, 116, 117, 124, ...
                    6, 13, 112, 12, 5];
chans           = {occipital, L_Temporal, R_Temporal, temporal, ... 
                   L_Parietal, R_Parietal, parietal, L_lat_frontal, ...
                   R_lat_frontal, Frontal_central, frontal};
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
        PSD = struct('mean',[],'log',[], 'log10',[]);
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
                        PSD(i_c).log10(i_f).(T_col)  = log10(band_mean);
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