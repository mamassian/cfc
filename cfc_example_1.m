% CONFIDENCE FORCED-CHOICE TOOLBOX  v0.3
%
% cfc_example_1
%   Confidence example for simulations and fit.
%   Only one task, using the method of constant stimuli.
%
% 09-OCT-2020 - Pascal Mamassian


% ----------------------
% -> prepare to simulate an experiment
% ----------------------

% -> define some parameters for simulated experiment
simul_params = struct;

% -> simulate method of constant stimuli
simul_params.method = 1;    

% -> list of stimuli with different difficulty levels
simul_params.sens_intens = -1.5:0.5:1.5;

% -> these are actually (confidence) pairs of trials
simul_params.nb_trials = 10000;


% ----------------------
% -> sensory parameters
% ----------------------

% -> noise of measurement (transduction): stdev
sens_noise1 = 1.0;

% -> sensory criterion
sens_crit1 = 0.25;


% ----------------------
% -> confidence parameters
% ----------------------

% -> added noise to combined (original + new sample)
conf_noise1 = 0.5;

% -> boost towards super-ideal: 0 = ideal;  1 = super-ideal
conf_boost1 = 0.2;

% -> confidence criterion
conf_crit1 = 0.0;

% -> confidence bias: overconfidence for task 2 (relative to task 1)
conf_bias = 1.0;    % 1.0 = no bias

% -> bias to choose 1st interval in confidence pair
intrvl_bias = 0.0;  % 0.0: no bias


% ----------------------
% -> simulate the experiment
% ----------------------

% -> create structure of parameters for simulated experiment
params_val = struct;
params_val.tasks_list = 1;
params_val.sens_noise = sens_noise1;
params_val.sens_crit  = sens_crit1;
params_val.conf_noise = conf_noise1;
params_val.conf_boost = conf_boost1;
params_val.conf_crit  = conf_crit1;
params_val.intrvl_bias = intrvl_bias;


% -> simulate the experiment and store data in 'raw_data' matrix
raw_data = cfc_simul_discrim(simul_params, params_val);

% -> prepare the data for plots and fit
grouped_data = cfc_group(raw_data);

% -> plot Type 1 psychometric functions
cfc_plot(grouped_data, 'type1_psychometric', true);

% -> plot simulated choices by reponses
cfc_plot(grouped_data, 'choi_by_resp', true);


% ----------------------
% -> fit model to data
% ----------------------

% -> create structure of chosen parameters to fit for model
%    start with '1' and increment
params_set = struct;
params_set.sens_noise = 1;
params_set.sens_crit  = 2;
params_set.conf_noise = 3;
params_set.conf_boost = 4;
params_set.conf_crit  = 0;
params_set.intrvl_bias = 0;

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
