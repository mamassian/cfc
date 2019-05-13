% -> confidence example for simulations and fit
%
%
% 23-MAY-2016 - Pascal Mamassian
% 05-NOV-2017 - pm: added fit
% 24-APR-2018 - pm: split integral space into subspaces to place discontinuities at boundaries
% 20-JAN-2019 - pm: adapted for cfc format
% 15-FEB-2019 - pm: generalize for different sensory tasks
% 13-MAR-2019 - pm: normalize confidence evidence by sensory sensitivity


close all;
rng_s = rng('shuffle');


% ----------------------
% -> sensory parameters
% ----------------------

% -> number of sensory tasks
% sens_tasks = 1;
sens_tasks = [1, 2];
nb_tasks = length(sens_tasks);

% -> noise of measurement (transduction): stdev
sens_noise1 = 1.0;
% sens_noise1 = 1.3;
% sens_noise2 = 1.0;   
sens_noise2 = 1.3;   

% -> sensory criterion
sens_crit1 = 0.0;
% sens_crit1 = 0.4;
sens_crit2 = 0.0;
% sens_crit2 = 0.4;

% ----------------------
% -> confidence parameters
% ----------------------

% -> added noise to combined (original + new sample)
% conf_noise1 = 0.0;
% conf_noise1 = 0.1;
conf_noise1 = 0.5;  % <<
% conf_noise1 = 1.0;
% conf_noise1 = 1.5;
% conf_noise1 = 10.0;

% conf_noise2 = 0.0;
% conf_noise2 = 0.1;
% conf_noise2 = 0.5;
conf_noise2 = 1.0;  % <<
% conf_noise2 = 1.5;
% conf_noise2 = 2.0;


% -> boost towards super-ideal: 0 = ideal;  1 = super-ideal
% conf_boost1 = 0.0;
% conf_boost1 = 0.1;
conf_boost1 = 0.2;  % <<
% conf_boost1 = 0.4;
% conf_boost1 = 1.0;

% conf_boost2 = 0.0;
conf_boost2 = 0.1;  % <<
% conf_boost2 = 0.2;
% conf_boost2 = 1.0;


% -> confidence criterion
conf_crit1 = sens_crit1;
conf_crit2 = sens_crit2;


% -> confidence bias: overconfidence for task 2 (relative to task 1)
conf_bias = 0.0;    % 0 = no bias
% conf_bias = 0.5;    % log(over_cfd2 / over_cfd1)


% -> bias to choose 1st interval in confidence pair as a log-odds in [0, 1]
intrvl_bias = 0.5;  % 0.5: no bias
% intrvl_bias = 0.65;  % bias > 0.5: in favour of interval 1

intrvl_log_odds = log(intrvl_bias / (1 - intrvl_bias));


% ----------------------
% -> add noise to the parameters

sens_noise1 = abs(sens_noise1 + randn*0.4);
sens_noise2 = abs(sens_noise2 + randn*0.4);
sens_crit1 = sens_crit1 + randn*0.4;
sens_crit2 = sens_crit2 + randn*0.4;
conf_noise1 = abs(conf_noise1 + randn*0.2);
conf_noise2 = abs(conf_noise2 + randn*0.2);
conf_boost1 = conf_boost1 + randn*0.1;
conf_boost2 = conf_boost2 + randn*0.1;
conf_boost1 = max(conf_boost1, 0.0);  conf_boost1 = min(conf_boost1, 1.0);
conf_boost2 = max(conf_boost2, 0.0);  conf_boost2 = min(conf_boost2, 1.0);
conf_crit1 = sens_crit1;
conf_crit2 = sens_crit2;
conf_bias = conf_bias + randn*0.4;
intrvl_bias = intrvl_bias + randn*0.1;
intrvl_bias = max(intrvl_bias, 0.0);  intrvl_bias = min(intrvl_bias, 1.0);


% ----------------------
% -> stimuli
% ----------------------

method_cst_stm = 0;    % simulate stimuli taken from a range (e.g. staircase)
% method_cst_stm = 1;    % simulate method of constant stimuli

if (method_cst_stm)
    % -> list of stimuli with different difficulty levels
    sens_intens1 = -1.5:0.75:1.5;
    sens_intens2 = -2.0:1.0:2.0;

    nb_sens_intens1 = length(sens_intens1);
    nb_sens_intens2 = length(sens_intens2);
else
    % -> ranges of difficulty levels
    sens_intens1_min = -1.5;
    sens_intens1_max =  1.5;
    sens_intens2_min = -1.5;
    sens_intens2_max =  1.5;
end


% -> these are actually (confidence) pairs of trials
% nb_trials = 5;
% nb_trials = 100;
% nb_trials = 200;
nb_trials = 1000;       % -> equivalent to an experiment of 1h
% nb_trials = 10000;    % <<
% nb_trials = 1000000;


% ----------------------
% -> fit
% ----------------------

if (~method_cst_stm)
    nb_bins = 4;
%     nb_bins = 5;
end

compute_efficiency = true;
% compute_efficiency = false;

split_conf_noise_boost = true;
% split_conf_noise_boost = false;

% add_conf_crit = true;
add_conf_crit = false;

% add_conf_bias = true;
add_conf_bias = false;

% add_intrvl_bias = true;
add_intrvl_bias = false;


% ----------------------
% -> other parameters
% ----------------------

% -> range to plot fit
xmin = -2.5;
xmax =  2.5;
fit_x = xmin:0.05:xmax;


% -----------------------------------------------------------------------
% -> you probably don't want to mess up this file beyond this point...
% -----------------------------------------------------------------------

tasks_nn = NaN(nb_trials, 2);   % stimulus tasks for intervals 1, 2
stims_nn = NaN(nb_trials, 2);   % stimulus intensities for intervals 1, 2
percs_nn = NaN(nb_trials, 2);   % percepts for intervals 1, 2
choic_nn = NaN(nb_trials, 2);   % confidence choice: [1, 0] if interval 1, [0, 1] if interval 2
sens_smpl = NaN(1, 2);  % the 2 sensory samples in a confidence pair used for Type 1
conf_smpl = NaN(1, 2);  % the 2 stimuli in a confidence pair used for Type 2
resp1 = NaN(1, 2);  % decision based on sensory evidence (Type 1)
resp2 = NaN(1, 2);  % decision based on confidence evidence (Type 2)
conf_evd = NaN(1, 2);       % confidence evidence

for tt = 1:nb_trials
    
    % -> interval kk of confidence pair
    for intrv = 1:2
        
        if (intrv == 1)
            % -> choose randomly one of the tasks
            ind_task = randi(nb_tasks);  % index of task for interval kk
        else
            % -> impose the other task in the other interval 
            ind_task = nb_tasks + 1 - tasks_nn(tt, 1);
        end
        tasks_nn(tt, intrv) = sens_tasks(ind_task); % which task

        % -> choose stimulus and make response
        switch tasks_nn(tt, intrv)
            case 1
                % -> choose randomly a stimulus in the set of possible stimuli
                if (method_cst_stm)
                    ind_intens = randi(nb_sens_intens1);   % index of stimulus for interval 'intrv'
                    stims_nn(tt, intrv) = sens_intens1(ind_intens);  % (noise-less) stimulus intensity
                else
                    stim_val = rand*(sens_intens1_max - sens_intens1_min) + sens_intens1_min;
                    stims_nn(tt, intrv) = stim_val;  % (noise-less) stimulus intensity
                end

                % -> take independent noisy samples of the stimuli
                sens_smpl(intrv) = stims_nn(tt, intrv) + randn * sens_noise1;

                % -> sensory decision based on side of sensory criterion
                resp1(intrv) = (sens_smpl(intrv) > sens_crit1);
                percs_nn(tt, intrv) = resp1(intrv);
        
                % -> add fraction super-ideal
                boosted_sens_smpl = (1 - conf_boost1)*sens_smpl(intrv) + conf_boost1*stims_nn(tt, intrv);

                % -> scale to sensory sensitivity
                scaled_smpl = boosted_sens_smpl / sens_noise1;

                % -> corrupt by Type 2 noise
                conf_smpl(intrv) = scaled_smpl + randn * conf_noise1;

                % -> confidence evidence
                conf_evd(intrv) = abs(conf_smpl(intrv) - conf_crit1);
                
                % -> sensory decision based on confidence sample
                resp2(intrv) = (conf_smpl(intrv) > sens_crit1);

             case 2
                % -> choose randomly a stimulus in the set of possible stimuli
                if (method_cst_stm)
                    ind_intens = randi(nb_sens_intens2);   % index of stimulus for interval 'intrv'
                    stims_nn(tt, intrv) = sens_intens2(ind_intens);  % (noise-less) stimulus intensity
                else
                    stim_val = rand*(sens_intens2_max - sens_intens2_min) + sens_intens2_min;
                    stims_nn(tt, intrv) = stim_val;  % (noise-less) stimulus intensity
                end

                % -> take independent noisy samples of the stimuli
                sens_smpl(intrv) = stims_nn(tt, intrv) + randn * sens_noise2;

                % -> sensory decision based on side of sensory criterion
                resp1(intrv) = (sens_smpl(intrv) > sens_crit2);
                percs_nn(tt, intrv) = resp1(intrv);
                
                % -> add fraction super-ideal
                boosted_sens_smpl = (1 - conf_boost2)*sens_smpl(intrv) + conf_boost2*stims_nn(tt, intrv);
                
                % -> scale to sensory sensitivity & confidence bias
                scaled_smpl = boosted_sens_smpl / sens_noise2 * exp(conf_bias);

                % -> corrupt by Type 2 noise
                conf_smpl(intrv) = scaled_smpl + randn * conf_noise2;

                % -> confidence evidence
                conf_evd(intrv) = abs(conf_smpl(intrv) - conf_crit2);
                
                % -> sensory decision based on confidence sample
                resp2(intrv) = (conf_smpl(intrv) > sens_crit2);
        end
    end

    % -> add interval bias as a log-odds
    pref_intrvl1 = log(conf_evd(1) / conf_evd(2)) + intrvl_log_odds;
    
    % -> are you more confident on 1st interval?
    choice = (pref_intrvl1 > 0);   % 1 if intrv1, 0 if intrv2
    choice_intrvl = 2 - choice;    % 1 if intrv1, 2 if intrv2
    
    % -> if sensory decision based on confidence evidence is inconsistent
    %    with perceptual decision, that perceptual decision might well be wrong
    if (resp2(choice_intrvl) ~= resp1(choice_intrvl))
        choice_intrvl = 3 - choice_intrvl;    % use the other interval as more confident
    end
    
    choic_nn(tt, choice_intrvl) = 1;
    choic_nn(tt, 3 - choice_intrvl) = 0;
end


% ------------------------------------------------------------------------
% -> group data by identical stimuli and responses
raw_data = [stims_nn, percs_nn, choic_nn, tasks_nn];
if (method_cst_stm)
    if (intrvl_bias == 0.5)
        [grouped_data, resp2_data] = cfc_group_03(raw_data);
    else
        [grouped_data, resp2_data] = cfc_group_03(raw_data, 'intrvl_sym', false);
    end
else
    if (intrvl_bias == 0.5)
        [grouped_data, resp2_data] = cfc_group_03(raw_data, 'bins', nb_bins);
    else
        [grouped_data, resp2_data] = cfc_group_03(raw_data, ...
            'bins', nb_bins, 'intrvl_sym', false);
    end
end

% -> plot simulated choices by reponses
cfc_plot_02(grouped_data, 'choi_by_resp', true);


% ------------------------------------------------------------------------
% -> quality of simulations

% -> order of parameters
params_model(01) = sens_noise1;
params_model(02) = sens_crit1;
params_model(03) = sens_noise2;
params_model(04) = sens_crit2;
params_model(05) = conf_noise1;
params_model(06) = conf_boost1;
params_model(07) = conf_noise2;
params_model(08) = conf_boost2;
params_model(09) = conf_crit1;
params_model(10) = conf_crit2;
params_model(11) = conf_bias;
params_model(12) = intrvl_bias;

[params_chosen, params_nn1] = cfc_core_07(grouped_data, params_model);
log_like = sum(log(params_chosen) .* grouped_data(:,5) + ...
    log(1 - params_chosen) .* grouped_data(:,6));
fprintf('Log-Likelihood = %7.3f\n', log_like);

cfc_plot_02(grouped_data, 'human_model', params_chosen);
xlabel('Parameterized Choice');     % change x-label
ylabel('Simulated Choice');         % change y-label
fig = findobj('type','figure');
fig2 = fig(1);


% ------------------------------------------------------------------------
% -> get fit

cfc_struct = cfc_fit_07(grouped_data, ...
    'efficiency', compute_efficiency, ...
    'conf_noise_boost', split_conf_noise_boost, ...
    'conf_crit', add_conf_crit, ...
    'conf_bias', add_conf_bias, ...
    'intrvl_bias', add_intrvl_bias);

fit_noise1 = cfc_struct.param_noise1;
fit_crit1 =  cfc_struct.param_crit1;

if (split_conf_noise_boost || add_conf_crit || add_conf_bias || add_intrvl_bias)
    fit_noise2 = cfc_struct.param_noise2;         % sdtev of noise for Type 2 decision (0 = ideal)
    fit_crit2 = cfc_struct.param_crit2;           % criterion for Type 2 decision (criterion)
    fit_boost2 = cfc_struct.param_boost2;         % fraction super-ideal (1 - fraction ideal)
    fit_nn1_full = cfc_struct.nn1_full;
    nn1_plot_full = fit_nn1_full;
end

% -> add best fit to check quality of fit
figure(fig2);
plot(cfc_struct.chosen_full, resp2_data, '+', 'MarkerSize', 12, ...
    'Color', [0.0, 0.8, 0.5], 'LineWidth', 2);


if (compute_efficiency)
    fit_efficiency = cfc_struct.param_efficiency;
    fprintf('Efficiency = %7.3f\n', fit_efficiency);

    % -> plot data against ideal
    cfc_plot_02(grouped_data, 'human_model', cfc_struct.chosen_ideal);
    
    fig = findobj('type','figure');
    set(fig(1), 'Position', [1215 600 430 340]);
    xlabel('Ideal Choice');         % change x-label
    ylabel('Simulated Choice');     % change y-label
end


% -> plot psychometric functions
cfc_plot_02(grouped_data, 'psychometric', cfc_struct.chosen_full);

% -> add cumulative Gaussian fit to unsorted data
for tt = 1:nb_tasks
    subplot(1, nb_tasks, tt);

    fit_y = normcdf(fit_x, fit_crit1(tt), fit_noise1(tt));
    plot(fit_x, fit_y, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3);
    xlim([-2.5, 2.5]);
end



% **** THE END ****
