% obs diag (state space RMSE at obs locations) RMSE, ensemble spread
% 2023/04/17

%% load .mat data 
clear all; 
load('obs_seq_tmp.mat')
[n_obs, n_ens, n_t] = size(PFF_prior);
%}

%% analysis

% transform back to the state space:
obs_ss       = log(obs);
truth_ss     = log(truth);
PFF_prior_ss = log(PFF_prior);
PFF_post_ss  = log(PFF_post );

prior_mean_PFF  = squeeze(mean( PFF_prior_ss, 2));
post_mean_PFF   = squeeze(mean( PFF_post_ss , 2));

prior_std_PFF  = squeeze( mean(std( PFF_prior_ss, 0, 2),1) );
post_std_PFF   = squeeze( mean(std( PFF_post_ss , 0, 2),1) );

prior_rmse_PFF  = zeros(n_t, 1);
post_rmse_PFF   = zeros(n_t, 1);

for i=1:n_t
    prior_rmse_PFF(i)  = sqrt( mean( (prior_mean_PFF(:,i) - truth_ss(:,i)).^2 ) );
    post_rmse_PFF(i)   = sqrt( mean( (post_mean_PFF(:,i) - truth_ss(:,i)).^2 ) );
end

% also do the rank histogram for the prior
rk_PFF = zeros(n_obs, n_t);

for i=1:n_t
    for j=1:n_obs
        [vals, rank]  = sort([truth_ss(j,i) squeeze(PFF_prior_ss(j,:,i))],'ascend');
        rk_PFF(j,i) = find(rank==1)-1; % rank is between 0 and n_ens
    end
end


exp_name='expslp5_1day_1yr';
da_config='PFFli_ker_08cap_infR3eakffg_ensavgHT';

fn_output = ['obs_loc_diag_',exp_name,'_',da_config,'.mat']
save(fn_output, 'prior_rmse_PFF', 'post_rmse_PFF', 'prior_std_PFF', 'post_std_PFF','rk_PFF')
%}


