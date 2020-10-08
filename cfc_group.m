% CONFIDENCE FORCED-CHOICE TOOLBOX  v0.3
%
% cfc_group
%   This function groups trials from a confidence forced-choice experiment.
%   The idea is to replace the raw data with group of trials that have the
%   same properties (same stimuli in both intervals, same responses, and
%   same tasks).
%
% INPUT:
%   'raw_data': matrix of:
%       1st col: stimulus intensity interval 1
%       2nd col: stimulus intensity interval 2
%       3rd col: perceptual decision interval 1 (1 = 'A',  0 = 'B')
%       4th col: perceptual decision interval 2 (1 = 'A',  0 = 'B')
%       5th col: confidence choice for interval 1 (1 = chosen,  0 = declined)
%       6th col: confidence choice for interval 2 (1 = chosen,  0 = declined)
%       7th col: stimulus task for interval 1 (optional)
%       8th col: stimulus task for interval 2 (optional)
%
% OPTIONAL PARAMETERS:
%   'bins': ask to bin stimuli for each task in (value) bins
%                (useful when stimulus can take any value within a range)
%
%   'merge_2intervals': consider that the two intervals are interchangeable
%                default value is 'false'
%                compact the data even further, so faster modelling
%                dangerous use if there is an interval bias (e.g. preference for 1st interval)
%
% OUTPUT:
%   'grouped_data': grouped data per (s1, s2, r1, r2):
%       1st col: stimulus intensity interval 1
%       2nd col: stimulus intensity interval 2
%       3rd col: perceptual decision interval 1
%       4th col: perceptual decision interval 2
%       5th col: nb of confidence choices for interval 1
%       6th col: nb of confidence choices for interval 2
%       7th col: stimulus task for interval 1
%       8th col: stimulus task for interval 2
%
%   'conf_choice_prob': list of Type 2 choice probability to select interval 1
%       (rather than interval 2), ordered in the same order as 'grouped_data'
%
%   'intrvl_bias_flag': z-score of interval bias
%
%   'intrvl_bias_pairs': list of interval biases for pairs of symmetric intervals
%       (for debugging purposes only)
%
%
% EXAMPLES OF USE:
%   grouped_data = cfc_group(raw_data);
%   grouped_data = cfc_group(raw_data, 'bins', 5);
%   grouped_data = cfc_group(raw_data, 'merge_2intervals', true);
%   [grouped_data, conf_choice_prob] = cfc_group(raw_data);
%   
%
% 21-SEP-2020 - pascal mamassian
 
% 29-FEB-2020 - pascal mamassian
% 10-JUN-2020 - pm: cleaned up


function [grouped_data, conf_choice_prob, intrvl_bias_flag, intrvl_bias_pairs] = ...
    cfc_group(raw_data, varargin)

    % -> default optional arguments
    dflt_nb_bins = 0;          % 0 is a message to not bin input data
    dflt_merge_2intervals = false;  % assume intervals are not interchangeable
    
    % -> parse all arguments
    ip = inputParser;
    valid_nb_bins = @(xx) isnumeric(xx) && isscalar(xx) && (xx > 0);
    addRequired(ip, 'raw_data', @isnumeric);
    addParameter(ip, 'bins', dflt_nb_bins, valid_nb_bins);
    addParameter(ip, 'merge_2intervals', dflt_merge_2intervals, @islogical);
    parse(ip, raw_data, varargin{:});
    nb_bins = ip.Results.bins;
    merge_2intervals = ip.Results.merge_2intervals;


    % -> initialize variables
    intrvl_bias_flag = NaN;
    nb_trials = size(raw_data, 1);

    % -> make sure the task is explicit for each interval
    if (size(raw_data, 2) < 7)
        % -> by default, assume that there is only one sensory task
        tasks_nn = repmat([1, 1], nb_trials, 1);
        raw_data = [raw_data, tasks_nn];
    end
    
    % -> check how many tasks were run
    tsk_vals = unique(raw_data(:, 7:8));
    nb_tasks = length(tsk_vals);
    
    
    % -> re-organize the data in bins
    if (nb_bins > 0)
        stm_data = [raw_data(:, [1, 7]); raw_data(:, [2, 8])];
        stm_mean = NaN(nb_tasks, nb_bins);
        inds_mat = NaN(nb_trials, 2);   % bin indices for the two intervals
        for tt = 1:nb_tasks
            inds3 = (stm_data(:,2) == tt);
            inds4 = find(inds3);
            stm_tsk_vals = stm_data(inds3);
            qq = quantile(stm_tsk_vals, nb_bins - 1);
            [~, ~, stm_tsk_bin] = histcounts(stm_tsk_vals, [-inf, qq, inf]);
            for nn = 1:nb_bins
                inds5 = stm_tsk_bin == nn;
                stm_mean(tt, nn) = mean(stm_tsk_vals(inds5));
                inds6 = inds4(inds5);
                inds1 = (inds6 <= nb_trials);    % indices in 1st interval
                inds2 = (inds6 > nb_trials);     % indices in 2nd interval
                inds_mat(inds6(inds1), 1) = nn;
                inds_mat(inds6(inds2)-nb_trials, 2) = nn;
            end
        end
                
        % -> replace raw stimulus value by mean of bin
        for rr = 1:nb_trials
            tsk1 = raw_data(rr, 7);
            tsk2 = raw_data(rr, 8);
            raw_data(rr, 1) = stm_mean(tsk1, inds_mat(rr, 1));
            raw_data(rr, 2) = stm_mean(tsk2, inds_mat(rr, 2));
        end
    end
    
    % -> extract the 'choice probability', i.e. prob to choose 1st interval
    type2_choice_prob = raw_data(:, 5) ./ sum(raw_data(:, 5:6), 2);
    
    
    % -> build the output matrix 'grouped_data' and the vector 'conf_choice_prob'
    [knd_vals, ~, knd_ic] = unique(raw_data(:, [1:4, 7:8]), 'rows');
    nb_knds = size(knd_vals, 1);
    conf_choice_prob = NaN(nb_knds, 1);
    grouped_data = NaN(nb_knds, 8);
    grouped_data(:, [1:4, 7:8]) = knd_vals;

    for ww = 1:nb_knds
        inds = (knd_ic == ww);     % indices that have the same kind
        type2_choice_vals = type2_choice_prob(inds);

        nn1 = sum(type2_choice_vals);
        nn0 = sum(1 - type2_choice_vals);
        
        conf_choice_prob(ww) = nn1 / (nn1 + nn0);
        grouped_data(ww, 5:6) = [nn1, nn0];
    end
    
    
    % -> ******************************************************************
    % -> more work to do if intervals are interchangeable
    
    intrvl_bias_pairs = NaN(size(grouped_data,1), 4);   % intrvl1, intrvl2, size1, size2
    
    if (merge_2intervals)

        % -> combine symmetric pairs (s2, s1) -> (s1, s2)
        nb_knds2 = 0;
        grouped_data2 = NaN(size(grouped_data));
        knd_vals2 = NaN(size(knd_vals));
        for ww = 1:nb_knds
            uu = knd_vals(ww, :);
            vv = [uu(2), uu(1), uu(4), uu(3), uu(6), uu(5)];    % entries 5&6 are actually 7&8
            pp = grouped_data(ww, :);
            qq = [pp(2), pp(1), pp(4), pp(3), pp(6), pp(5), pp(8), pp(7)];
            [aa, bb] = ismember(vv, knd_vals2, 'rows');
            if (aa)
                new_vals = pp([6,5]);
                grouped_data2(bb, [5,6]) = grouped_data2(bb, [5,6]) + new_vals;
                intrvl_bias_pairs(bb, 4) = sum(new_vals);
                intrvl_bias_pairs(bb, 2) = new_vals(1) / intrvl_bias_pairs(bb, 4);

            else
                nb_knds2 = nb_knds2 + 1;
                knd_vals2(nb_knds2, :) = uu;
                
                grouped_data2(nb_knds2, :) = pp;
                        
                intrvl_bias_pairs(nb_knds2, 3) = sum(grouped_data2(nb_knds2, [5,6]));
                intrvl_bias_pairs(nb_knds2, 1) = grouped_data2(nb_knds2, 5) / intrvl_bias_pairs(nb_knds2, 3);
            end
        end
        
        % -> shrink matrices to actual values
        grouped_data = grouped_data2(1:nb_knds2, :);
        intrvl_bias_pairs = intrvl_bias_pairs(1:nb_knds2, :);
        
        % -> clean the order
        for ww = 1:nb_knds2
            pp = grouped_data(ww, :);
            qq = [pp(2), pp(1), pp(4), pp(3), pp(6), pp(5), pp(8), pp(7)];
            
            rr = intrvl_bias_pairs(ww, :);
            ss = [rr(2), rr(1), rr(4), rr(3)];
            
            if (pp(8) < pp(7))
                grouped_data(ww, :) = qq;
                intrvl_bias_pairs(ww, :) = ss;
            elseif (pp(7) == pp(8))
                if (pp(2) < pp(1))
                    grouped_data(ww, :) = qq;
                    intrvl_bias_pairs(ww, :) = ss;
                end
            end

        end
        

        % -> check whether there is an interval bias and issue a warning if yes
        prob1 = intrvl_bias_pairs(:, 1);
        prob2 = intrvl_bias_pairs(:, 2);
        
        diffe = norminv(prob1) - norminv(prob2);
        zscore = nanmean(diffe) / nanstd(diffe);
        intrvl_bias_flag = zscore;
        if (abs(zscore) > 1.0)
            fprintf('\nWARNING: the data contain an interval bias (z-score = %7.3f)\n', zscore);
            fprintf('You may consider keeping the two intervals separated, i.e. use:\n');
            fprintf('grouped_data = cfc_group(raw_data, ''merge_2intervals'', false);\n\n');
        end
        
        conf_choice_prob = grouped_data(:, 5) ./ sum(grouped_data(:, 5:6), 2);
    end
    
    
end

% -> THE END

