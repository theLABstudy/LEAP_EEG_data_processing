% happe_segment() - A helper script for HAPPE that segments the data
%                   according to user-specification and paradigm. Does not
%                   work on data that has already been segmented.
%
% Usage: 
%   >> EEG = happe_segment(EEG, params)
%
% Inputs:
%   EEG    - The EEG, in EEGLAB (Delorme & Makeig, 2004) format, to
%            segment.
%   params - The struct containing all the relevent information for 
%            segmentation, including segmentation bounds.
%
% Outputs:
%   EEG    - The EEG, in EEGLAB (Delorme & Makeig, 2004) format, after 
%            segmentation.
%
% Author: L.J. Gabard-Durnam, PINE Lab at Northeastern University, 2021
%         A.D. Monachino, PINE Lab at Northeastern University, 2022
%
% This file is part of HAPPE.
% Copyright 2018, 2021 Alexa Monachino, Kelsie Lopez, Laurel Gabard-Durnam
%
% HAPPE is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
% 
% HAPPE is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
% details.
% 
% You should have received a copy of the GNU General Public License along
% with HAPPE. If not, see <https://www.gnu.org/licenses/>.

function EEG = happe_segment_v4(EEG, params)
if EEG.trials == 1
    fprintf('Segmenting...\n') ;
    % TASK/ERP - SEGMENT USING TAGS
    if params.paradigm.task
        % Check to ensure all event types are char vectors
        for indx=1:size(EEG.event,2)
            if ~ischar(EEG.event(indx).type)
                EEG.event(indx).type = num2str(EEG.event(indx).type) ;
            end
        end

        % For ERPs, adjust the event latencies by the user-specified
        % offset.
        if params.paradigm.ERP.on
            sampOffset = EEG.srate*params.segment.offset/1000 ;
            for i = 1:size(EEG.event, 2); EEG.event(i).latency = ...
                    EEG.event(i).latency + sampOffset ; end
        end
        
        % Try segmenting the data using user-specified tags, catching any
        % errors that occur in the process. If unable to find any of the
        % tags, will throw a specific error to help the user troubleshoot.
        try
            %% TODO: GP 2022-10-04
            %   - adapt beapp_extract_condition_labels to update
            %     EEG.event.type with EEG.event.mffkey_cond from most
            %     recent event with EEG.event.code == 'TRSP'
            %      - should add setting to inputParameters to indicate
            %        whether this is necessary
            %   - potentially, add something like add_events_eeglab_struct
            %     to get event info from table into EEG struct
            %   - pop_epoch should then work. may need to tell EEGLAB to
            %     search for type field, rather than code
            if ~isfield(params.paradigm,'use_cond_labels')
                params.paradigm.use_cond_labels = 0;
            end
            if params.paradigm.use_cond_labels
                [evt_info,cba,~] = extract_condition_labels( EEG.filename, ...
                    'event-related', [], {EEG.event}, params.paradigm );
                EEG.event = evt_info{1};
                EEG = pop_epoch(EEG, cba.Condition_Name', ...
                    [params.segment.start, params.segment.end], 'verbose', ...
                    'no', 'epochinfo', 'yes') ;
            else
                EEG = pop_epoch(EEG, params.paradigm.onsetTags, ...
                    [params.segment.start, params.segment.end], 'verbose', ...
                    'no', 'epochinfo', 'yes') ;
            end
        catch ME
            if strcmp(ME.message, ['pop_epoch(): empty epoch range (no ' ...
                    'epochs were found).'])
                error('HAPPE:noTags', ['ERROR: No specified event/onset ' ...
                    'tags in the data.']) ;
            else; rethrow(ME) ;
            end
        end
        
    % BASELINE/RESTING - SEGMENT USING LENGTH: Use the user-specified
    % segment length to epoch the data.
    else
        % GP added 2023-03-06: fix baseline segmentation
        EEG = pop_selectevent(EEG,'type',{'boundary','a_on','a_of'},    ...
              'deleteevents','on');
        bl_evts = [find(strcmp({EEG.event.type},'a_of')) find(strcmp({EEG.event.type},'a_on'))];
        evt_start = max( min(bl_evts)-1, 1); 
        pnt_start = EEG.event(evt_start).latency - 1; % added 03.31.2023 by GG & GP
        EEG = eeg_eegrej(EEG,[1 pnt_start]);
        bl_evts = [find(strcmp({EEG.event.type},'a_of')) find(strcmp({EEG.event.type},'a_on'))];
        evt_end = min( max(bl_evts)+1, length(EEG.event));
        pnt_end = (EEG.event(evt_end).latency + 1); % added 03.31.2023 by GG & GP
        if pnt_end > 10000
            EEG = eeg_eegrej(EEG,[pnt_end EEG.pnts]);
        end

% if no event markers after last a_on, check if there is a lot of data after
% to decide whether to include as the end of baseline, or other data to cut
% off

% for edge case scenario, the following line is more general and should
% catch it as well as correct segmentation

        EEG = eeg_eegrej(EEG,[pnt_start + 180*EEG.srate, EEG.pnts]); % experimental - added 04/04/23 by GG

        EEG = eeg_regepochs(EEG, 'recurrence', params.segment.length, ...
            'limits', [0 params.segment.length], 'rmbase', [NaN]) ;
    end
else; fprintf('Cannot segment data that has already been segmented.\n') ;
end

end

% select events from baseline period from input param onsetTags
function EEG = pop_select_bl_events(EEG,params)
    n_evt = length(EEG.event);
    % get indeces of events by type
    bl_inds = find_evt_inds(EEG,params.paradigm.onsetTags);
    if isempty(bl_inds)
        % if empty, try older event tags
        bl_inds = find_evt_inds(EEG,{'a_of'});
    end
    if isempty(bl_inds), EEG = []; return; end
    X_inds  = find_evt_inds(EEG,{'X'});
    other_task_inds = setdiff(1:n_evt,[bl_inds X_inds]);
    % serpate rows into baseline and non-baseline events
    non_bl = other_task_inds > bl_inds(1);
    % if non_bl is 0s, baseline after other tasks
    if ~non_bl
        non_bl_rng = 1:(bl_inds(1)-1);
    else
        after_bl_inds = other_task_inds(non_bl);
        % get non/baseline ranges
        bl_rng = bl_inds(1):(after_bl_inds(1)-1);
        non_bl_rng = setdiff(1:n_evt,bl_rng);
    end
    % remove epochs with non baseline events
    rm_epochs = unique([EEG.event(non_bl_rng).epoch]);
    EEG = pop_select(EEG,'notrial',rm_epochs);
end

% find indices of events by input type
function evt_inds = find_evt_inds(EEG,evt_type)
    evt_inds = find( strcmp({EEG.event.type},evt_type) );
end

% function EEG = select_bl_events(EEG,params)
%     % get EEG.event row indeces for event tags of interest
%     event_ind = find(strcmp({EEG.event.code},params.paradigm.onsetTags));
%     % get indeces for events following events of interest
%     event_ind = unique([event_ind event_ind+1]);
%     % keep only these rows in event table
%     EEG.event = EEG.event(event_ind);
%     % subtract latency of first row from other rows s.t. t=0 -->
%     % first event
%     latencies = [EEG.event.latency]; latencies = latencies - latencies(1);
%     for i = 1:length(EEG.event), EEG.event(i).latency = latencies(i); end
%     % update xmax to reflect latency of last kept event
%     EEG.xmax = EEG.event(end).latency/EEG.srate;
% end