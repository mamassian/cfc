% core function for confidence forced-choice prediction
%
%   computes the probability of choosing interval 1 for a given trial
%   of a confidence forced-choice experiment
%
% INPUT:
%   wrap_data: grouped data per (s1, s2, r1, r2):
%       1st col: stimulus intensity interval 1
%       2nd col: stimulus intensity interval 2
%       3rd col: perceptual decision interval 1
%       4th col: perceptual decision interval 2
%       5th col: nb of confidence choice for interval 1
%       6th col: nb of confidence choice for interval 2
%       7th col: stimulus task for interval 1 (optional)
%       8th col: stimulus task for interval 2 (optional)
%   model_params: vector of model parameters:
%       1. sens_noise1: sensory (Type 1) sdtev of noise for task 1 (0 = perfectly sensitive)
%       2. sens_crit1:  sensory (Type 1) criterion for task 1
%       3. sens_noise2: sensory (Type 1) sdtev of noise for task 2
%       4. sens_crit2:  sensory (Type 1) criterion for task 2
%       5. conf_noise1: confidence (Type 2) sdtev of noise for task 1 (0 = ideal)
%       6. conf_boost1: fraction super-ideal for task 1 (0 = ideal, 1 = super-ideal)
%       7. conf_noise2: confidence (Type 2) sdtev of noise for task 2 (0 = ideal)
%       8. conf_boost2: fraction super-ideal for task 2 (0 = ideal, 1 = super-ideal)
%       9. conf_crit1:  confidence (Type 2) criterion for task 1 (criterion)
%       10. conf_crit2:  confidence (Type 2) criterion for task 2 (criterion)
%       11. conf_bias:   overconfidence for task 2 (relative to task 1)
%       12. intrvl_bias: bias in favour of interval 2
%
% OUTPUT:
%   type2_resp: Type 2 choice probability to select interval 1 (rather than 2)
%   nn1_lst: nb of confidence choices for intervals 1 and 2
%
%
% 05-NOV-2017 - pascal mamassian
% 24-APR-2018 - pm: split integral space into subspaces to place discontinuities at boundaries
% 29-DEC-2018 - pm: take into account type 2 criterion
% 06-JAN-2019 - pm: computes all 4 type 2 choice probabilities
% 07-FEB-2019 - pm: expand to all stimuli of an experiment
% 10-FEB-2019 - pm: introduce 'small triangle' for type 2 criterion
% 15-FEB-2019 - pm: generalize for different sensory tasks
% 13-MAR-2019 - pm: normalize confidence evidence by sensory sensitivity
% 21-MAR-2019 - pm: do not fit second noise2 parameter (undetermined)
% 29-MAR-2019 - pm: take care of extreme probability values


function [type2_resp, nn1_lst] = cfc_core_08(wrap_data, model_params)

% -> nb. of cases of (s1, s2, r1, r2, t1, t2)
nb_knds = size(wrap_data, 1);
type2_resp = NaN(nb_knds, 1);

if (size(wrap_data, 2) < 7)
    % by default, assume that there is only one sensory task
    tasks_nn = repmat([1, 1], nb_knds, 1);
    wrap_data = [wrap_data, tasks_nn];
end


% -> tolerances for Matlab integral functions (smaller is better)
atol = 1e-6;
% rtol = 1e-3;  % <<
rtol = 1e-2;

% -> display on-going information
% verbose_on = 1;
verbose_on = 0;

% -> order of parameters
ind_1noise1 = 1;
ind_1crit1  = 2;
ind_2noise1 = 3;
ind_2crit1  = 4;
ind_1noise2 = 5;
ind_1boost  = 6;
ind_2noise2 = 7;
ind_2boost  = 8;
ind_1crit2  = 9;
ind_2crit2 = 10;
ind_bias2  = 11;
ind_intrvl = 12;

% -> Type 1 parameters
sens_noise_task = NaN(1, 2);  % sensory noise, per task
sens_crit_task = NaN(1, 2);   % sensory criterion, per task

sens_noise_task(1) = model_params(ind_1noise1);
sens_crit_task(1)  = model_params(ind_1crit1);
sens_noise_task(2) = model_params(ind_2noise1);
sens_crit_task(2)  = model_params(ind_2crit1);


% -> Type 2 parameters
conf_noise_task = NaN(1, 2);  % confidence noise, per task
conf_boost_task = NaN(1, 2);  % confidence boost, per task
conf_crit_task = NaN(1, 2);  % confidence criterion, per task

conf_noise_task(1) = model_params(ind_1noise2);
conf_boost_task(1) = model_params(ind_1boost);
conf_noise_task(2) = model_params(ind_2noise2);
conf_boost_task(2) = model_params(ind_2boost);
conf_crit_task(1)  = model_params(ind_1crit2);
conf_crit_task(2)  = model_params(ind_2crit2);

conf_bias_task2    = model_params(ind_bias2);
intrvl_bias        = model_params(ind_intrvl);

% -> rescale by confidence bias
conf_scale = NaN(1, 2); 
conf_scale(1) = 1 / sens_noise_task(1);  % rescale by sensory sensitivity
conf_scale(2) = exp(conf_bias_task2) / sens_noise_task(2);

% -> transform interval bias into log-odds
intrvl_log_odds = log(intrvl_bias / (1 - intrvl_bias));

% -> constrain (2noise2 = 1noise2) if there are 2 tasks but 2noise2 is NaN
if (~isnan(sens_noise_task(2)) && isnan(conf_noise_task(2)))
    conf_noise_task(2) = conf_noise_task(1);
end


% -> bounds for the integrals, per task
% intsz = Inf .* ones(1, nb_tasks);     % in theory what should be used, but takes too long
intsz_task = 10.0 .* sens_noise_task;



for kk = 1:nb_knds

    % -> Type 1 tasks (1 or 2)
    tsk1 = wrap_data(kk, 7);    % task for interval 1
    tsk2 = wrap_data(kk, 8);    % task for interval 2
    
    % -> Stimulus intensities
    mu1 = wrap_data(kk, 1);
    x1_min = mu1 - intsz_task(tsk1);
    x1_max = mu1 + intsz_task(tsk1);

    mu2 = wrap_data(kk, 2);
    x2_min = mu2 - intsz_task(tsk2);
    x2_max = mu2 + intsz_task(tsk2);

    % -> Type 1 responses (1 = 'A', 0 = 'B')
    resp1 = wrap_data(kk, 3);
    resp2 = wrap_data(kk, 4);

    sens_noise1 = sens_noise_task(tsk1);
    sens_crit1 = sens_crit_task(tsk1);
    conf_noise1 = conf_noise_task(tsk1);
    conf_crit1 = conf_crit_task(tsk1);
    conf_boost1 = conf_boost_task(tsk1);
    
    sens_noise2 = sens_noise_task(tsk2);
    sens_crit2 = sens_crit_task(tsk2);
    conf_noise2 = conf_noise_task(tsk2);
    conf_crit2 = conf_crit_task(tsk2);
    conf_boost2 = conf_boost_task(tsk2);
    
    conf_scale1 = conf_scale(tsk1);
    conf_scale2 = conf_scale(tsk2);
    
    
    % -> probability density function of (xx, yy):
    %    P(s1, s2)
    noisy1_hdl = @(xx, yy) ...
        noisy1_fun(xx, yy, mu1, mu2, sens_noise1, sens_noise2);


    % -> probability to choose interval 1 for a pair of (smpl1, smpl2):
    %    P(ch=1 | s1, s2)
    choice2_hdl = @(xx, yy) ...
        choice2_fun(xx, yy, mu1, mu2, ...
                sens_crit1, conf_noise1, conf_crit1, conf_boost1, conf_scale1, ...
                sens_crit2, conf_noise2, conf_crit2, conf_boost2, conf_scale2);

    % -> probability to choose interval 1 for a pair of (smpl1, smpl2):
    %    P(ch=1)
    choice1_hdl = @(xx, yy) ...
        noisy1_hdl(xx, yy) .* choice2_hdl(xx, yy);



    % -> go through all 4 possible scenarii
    if (resp1 == 1)
        if (resp2 == 1)
            % -> Type 1 responses were 'A' and 'A'
            % -> P(r1=1, r2=1)
            joint_prob = (1 - normcdf(sens_crit1, mu1, sens_noise1)) * ...
                (1 - normcdf(sens_crit2, mu2, sens_noise2));
            % -> P(ch=1 | r1=1, r2=1)
            choose_intrvl1 = integral2(choice1_hdl, sens_crit1, x1_max, ...
                sens_crit2, x2_max,'AbsTol',atol,'RelTol',rtol);
        else
            % -> Type 1 responses were 'A' and 'B'
            joint_prob = (1 - normcdf(sens_crit1, mu1, sens_noise1)) * ...
                normcdf(sens_crit2, mu2, sens_noise2);
            choose_intrvl1 = integral2(choice1_hdl, sens_crit1, x1_max, ...
                x2_min, sens_crit2,'AbsTol',atol,'RelTol',rtol);
        end
    else
        if (resp2 == 1)
            % -> Type 1 responses were 'B' and 'A'
            joint_prob = normcdf(sens_crit1, mu1, sens_noise1) * ...
                (1 - normcdf(sens_crit2, mu2, sens_noise2));
            choose_intrvl1 = integral2(choice1_hdl, x1_min, sens_crit1, ...
                sens_crit2, x2_max,'AbsTol',atol,'RelTol',rtol);
        else
            % -> Type 1 responses were 'B' and 'B'
            joint_prob = normcdf(sens_crit1, mu1, sens_noise1) * ...
                normcdf(sens_crit2, mu2, sens_noise2);
            choose_intrvl1 = integral2(choice1_hdl, x1_min, sens_crit1, ...
                x2_min, sens_crit2,'AbsTol',atol,'RelTol',rtol);
        end
    end

    prob_intrvl1 = choose_intrvl1 / joint_prob;
    
    if (prob_intrvl1 < 1e-3)
        new_prob_intrvl1 = 0.0;
    elseif (prob_intrvl1 > (1.0 - 1e-3))
        new_prob_intrvl1 = 1.0;
    else
        evid_intrvl1 = log(prob_intrvl1 ./(1 - prob_intrvl1)) + intrvl_log_odds;
        new_prob_intrvl1 = 1 / (1/exp(evid_intrvl1) + 1);
    end
    type2_resp(kk) = new_prob_intrvl1;

    if (verbose_on)
        if (mod(kk, 10) == 0)
            fprintf('. ');
        else
            fprintf('.');
        end
    end

end
if (verbose_on)
    fprintf('\n');
end

nn1_choices = wrap_data(:, 5:6);
nn1_tmp1 = sum(nn1_choices, 2);
nn1_tmp2 = round(nn1_tmp1 .* type2_resp);
nn1_lst = [nn1_tmp2, nn1_tmp1 - nn1_tmp2];


% ------------------------------------------------------------------------
% -> internal functions for model

% -> probability density function of (smpl1, smpl2): 
%    P(s1, s2)
function prob = noisy1_fun(smpl1, smpl2, mu1, mu2, sd1, sd2) 
    prob = normpdf(smpl1, mu1, sd1) .* normpdf(smpl2, mu2, sd2);
end


% -> probability to choose interval 1 for a pair of (smpl1, smpl2):
%    P(ch=1 | s1, s2)
function prob = choice2_fun(smpl1, smpl2, mu1, mu2, ...
        s_crit1, c_noise1, c_crit1, c_boost1, c_scale1, ...
        s_crit2, c_noise2, c_crit2, c_boost2, c_scale2) 
    
    % -> this function is often called within a 2D integral, i.e. with
    %    a matrix of samples
    [smplx_nb, smply_nb] = size(smpl1);
    prob = NaN(smplx_nb, smply_nb);
    
    if (c_noise1 == 0)
        c_noise1 = 1e-6;
%         c_noise1 = 1e-3;
    end
    if (c_noise2 == 0)
        c_noise2 = 1e-6;
%         c_noise2 = 1e-3;
    end
        
    for smplx_ii=1:smplx_nb
        for smply_ii=1:smply_nb
            smpl_val1 = smpl1(smplx_ii, smply_ii);
            smpl_val2 = smpl2(smplx_ii, smply_ii);
            
            % -> perceptual decision (1 = 'A', 0 = 'B') for intervals 1 & 2
            intrvl1_dec = (smpl_val1 > s_crit1);
            intrvl2_dec = (smpl_val2 > s_crit2);

            % -> sliced means for intervals 1 & 2
            c_Ox = c_scale1 * (smpl_val1 + (mu1 - smpl_val1) * c_boost1);
            c_Oy = c_scale2 * (smpl_val2 + (mu2 - smpl_val2) * c_boost2);

            
            % -> place relevant points
            o_Ox = s_crit1;
            o_Oy = s_crit2;

            p_Ox = - s_crit2 + c_crit1 + c_crit2;
            p_Oy =   s_crit2;

            q_Ox =   s_crit1;
            q_Oy = - s_crit1 + c_crit1 + c_crit2;

            r_Ox = c_crit1;
            r_Oy = c_crit2;

            s_Ox = s_crit1;
            s_Oy = s_crit2 - c_crit1 + c_crit2;

            t_Ox = s_crit1 + c_crit1 - c_crit2;
            t_Oy = s_crit2;

            % -> rescale the axes so that the confidence noise distribution
            %    is bivariate Normal with unit variance (and no covariance)
            c_Oxr = c_Ox / c_noise1;     c_Oyr = c_Oy / c_noise2;
            o_Oxr = o_Ox / c_noise1;     o_Oyr = o_Oy / c_noise2;
            p_Oxr = p_Ox / c_noise1;     p_Oyr = p_Oy / c_noise2;
            q_Oxr = q_Ox / c_noise1;     q_Oyr = q_Oy / c_noise2;
            r_Oxr = r_Ox / c_noise1;     r_Oyr = r_Oy / c_noise2;
            s_Oxr = s_Ox / c_noise1;     s_Oyr = s_Oy / c_noise2;
            t_Oxr = t_Ox / c_noise1;     t_Oyr = t_Oy / c_noise2;


            % -> rotate coordinates Oxy clockwise to Ouv so that Ou
            %    is parallel to the separating line between AB and BA
            slope_yx1 = - c_noise1 / c_noise2;     % slope of that line
            [~, c_Ov1] = cfc_rotate(c_Oxr, c_Oyr, slope_yx1);
            [~, r_Ov1] = cfc_rotate(r_Oxr, r_Oyr, slope_yx1);
            
            % -> rotate coordinates Oxy clockwise to Ou'v' so that Ov'
            %    is parallel to the separating line between AA and BB
            slope_yx2 = - c_noise2 / c_noise1;     % slope of that line
            [c_Ou2, ~] = cfc_rotate(c_Oxr, c_Oyr, slope_yx2);
            [r_Ou2, ~] = cfc_rotate(r_Oxr, r_Oyr, slope_yx2);


            % -> compute "small triangle" if there is a conf_crit
            if ((c_crit1 ~= s_crit1) || (c_crit2 ~= s_crit2))
                % -> WARNING: not fully tested !!!
                
                % -> treat separately {AB, BA} and {AA, BB} 
                if (intrvl1_dec ~= intrvl2_dec)

                    add_triangle = sign(c_crit1 - s_crit1) * ...
                        sign(c_crit2 - s_crit2);

                    triangle_vertices = [o_Oxr, o_Oyr; p_Oxr, p_Oyr; q_Oxr, q_Oyr];

                    % -> shift relative to center of bivariate Gaussian                
                    triangle_vertices = triangle_vertices - [c_Oxr, c_Oyr];
                    
                    triangle_area = bivariate_any_triangle(triangle_vertices);
                    triangle_area = triangle_area * add_triangle;
                else
                    add_triangle1 = sign(c_crit2 - s_crit2);
                    add_triangle2 = sign(c_crit1 - s_crit1);

                    triangle_vertices1 = [t_Oxr, t_Oyr; p_Oxr, p_Oyr; r_Oxr, r_Oyr];
                    triangle_vertices1 = triangle_vertices1 - [c_Oxr, c_Oyr];
                    triangle_area1 = bivariate_any_triangle(triangle_vertices1);
                    triangle_area1 = triangle_area1 * add_triangle1;

                    triangle_vertices2 = [s_Oxr, s_Oyr; r_Oxr, r_Oyr; q_Oxr, q_Oyr];
                    triangle_vertices2 = triangle_vertices2 - [c_Oxr, c_Oyr];
                    triangle_area2 = bivariate_any_triangle(triangle_vertices2);
                    triangle_area2 = triangle_area2 * add_triangle2;
                end
                
            else
                triangle_area = 0.0;
                triangle_area1 = 0.0;
                triangle_area2 = 0.0;
            end     % end of "small triangle"
    
            
            % -> go through all 4 possible scenarii
            if (intrvl1_dec == 1)
                if (intrvl2_dec == 1)
                    % -> Type 1 responses were 'A' and 'A'
                    prob(smplx_ii, smply_ii) = 1 - normcdf(r_Ou2, c_Ou2, 1.0);
                    prob(smplx_ii, smply_ii) = prob(smplx_ii, smply_ii) - triangle_area1 + triangle_area2;
                else
                    % -> Type 1 responses were 'A' and 'B'
                    prob(smplx_ii, smply_ii) = 1 - normcdf(r_Ov1, c_Ov1, 1.0);
                    prob(smplx_ii, smply_ii) = prob(smplx_ii, smply_ii) + triangle_area;
                end
            else
                if (intrvl2_dec == 1)
                    % -> Type 1 responses were 'B' and 'A'
                    prob(smplx_ii, smply_ii) = normcdf(r_Ov1, c_Ov1, 1.0);
                    prob(smplx_ii, smply_ii) = prob(smplx_ii, smply_ii) - triangle_area;
                else
                    % -> Type 1 responses were 'B' and 'B'
                    prob(smplx_ii, smply_ii) = normcdf(r_Ou2, c_Ou2, 1.0);
                    prob(smplx_ii, smply_ii) = prob(smplx_ii, smply_ii) + triangle_area1 - triangle_area2;
                end
            end
            
        end
    end
    

    % -> rotate coordinates Oxy clockwise to Ouv so that Ou
    %    is parallel to the line that has slope 'aa'
    function [mu, mv] = cfc_rotate(mx, my, aa)
        cc = sqrt(1 + aa^2);
        mu = (     mx + aa* my) / cc;
        mv = (-aa* mx +     my) / cc;
    end

end

end