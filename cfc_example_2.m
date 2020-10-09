% CONFIDENCE FORCED-CHOICE TOOLBOX  v0.3
%
% cfc_example_2
%   Confidence example for simulations and fit.
%   Two tasks, with a confidence bias.
%
% 09-OCT-2020 - Pascal Mamassian


% ----------------------
% -> prepare to simulate an experiment
% ----------------------

% -> define some parameters for simulated experiment
simul_params = struct;

% -> number of sensory tasks
sens_tasks = [1, 2];
nb_tasks = length(sens_tasks);

for task_no = 1:nb_tasks
    
    method_cst_stm = 1;    % simulate method of constant stimuli

    % -> list of stimuli with different difficulty levels
    if (task_no == 1)
        sens_intens = -2.25:0.75:2.25;
    else
        sens_intens = -1.5:0.75:1.5;
    end

    simul_params(task_no).sens_intens = sens_intens;

    % -> these are actually (confidence) pairs of trials
    nb_trials = 10000;

    simul_params(task_no).method = method_cst_stm;
    simul_params(task_no).nb_trials = nb_trials;

end


% ----------------------
% -> sensory parameters
% ----------------------

% -> noise of measurement (transduction): stdev
sens_noise1 = 2.0;
sens_noise2 = 1.0;

% -> sensory criterion
sens_crit1 = 0.75;
sens_crit2 = 0.25;


% ----------------------
% -> confidence parameters
% ----------------------

% -> added noise to combined (original + new sample)
conf_noise1 = 0.5;
conf_noise2 = 1.0;

% -> boost towards super-ideal: 0 = ideal;  1 = super-ideal
conf_boost1 = 0.2;
conf_boost2 = 0.1;

% -> confidence criterion
conf_crit1 = 0.0;
conf_crit2 = 0.0;

% -> confidence bias: overconfidence for task 2 (relative to task 1)
% conf_bias = 1.0;    % 1.0 = no bias
conf_bias = 1.5;    % over_cfd2 / over_cfd1

% -> bias to choose 1st interval in confidence pair
intrvl_bias = 0.0;  % 0.0: no bias


% ----------------------
% -> simulate the experiment
% ----------------------

% -> create structure of parameters for simulated experiment
params_val = struct;
params_val.tasks_list = 1:nb_tasks;
params_val.sens_noise = [sens_noise1, sens_noise2];
params_val.sens_crit  = [sens_crit1, sens_crit2];
params_val.conf_noise = [conf_noise1, conf_noise2];
params_val.conf_boost = [conf_boost1, conf_boost2];
params_val.conf_crit  = [conf_crit1, conf_crit2];
params_val.intrvl_bias = intrvl_bias;
params_val.conf_bias  = [1.0, conf_bias];


% -> simulate the experiment and store data in 'raw_data' matrix
raw_data = cfc_simul_discrim(simul_params, params_val);

% -> prepare the data for plots and fit
[grouped_data, conf_choice_prob] = cfc_group(raw_data);

% -> plot Type 1 psychometric functions
cfc_plot(grouped_data, 'type1_psychometric', true);



% ----------------------
% -> fit model to data
% ----------------------

% -> create structure of chosen parameters to fit for model
%    start with '1' and increment
params_set = struct;
params_set.sens_noise = [1, 2];
params_set.sens_crit  = [3, 4];
params_set.conf_noise = [5, 6];
params_set.conf_boost = [7, 8];
params_set.conf_crit  = [0, 0];
params_set.intrvl_bias = 0;
params_set.conf_bias  = [0, 9];

cfc_struct = cfc_fit(grouped_data, 'model_parameters', params_set);


% -> collect estimated parameters
fit_noise1 = cfc_struct.sens_noise;
fit_crit1 =  cfc_struct.sens_crit;

fit_noise2 = cfc_struct.conf_noise;         % sdtev of noise for Type 2 decision (0 = ideal)
fit_crit2 = cfc_struct.conf_crit;           % criterion for Type 2 decision (criterion)
fit_boost2 = cfc_struct.conf_boost;         % fraction super-ideal (1 - fraction ideal)


% -> plot choices from simulated observer against best fit of model
cfc_plot(grouped_data, 'human_model', cfc_struct.choice_prob_model);

% -> plot Type 2 psychometric functions
cfc_plot(grouped_data, 'type2_psychometric', cfc_struct.choice_prob_model);


% -> THE END
