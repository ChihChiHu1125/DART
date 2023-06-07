% obs diag (state space RMSE at obs locations) RMSE, ensemble spread
% 2023/04/17

%% load .mat data 
clear all; 
load('obs_seq_tmp.mat')
[n_obs, n_ens, n_t] = size(EAKF_prior);
%}

%% analysis

% transform back to the state space:
obs_ss       = log(obs);
truth_ss     = log(truth);
EAKF_prior_ss = log(EAKF_prior);
EAKF_post_ss  = log(EAKF_post );

prior_mean_EAKF  = squeeze(mean( EAKF_prior_ss, 2));
post_mean_EAKF   = squeeze(mean( EAKF_post_ss , 2));

prior_std_EAKF  = squeeze( mean(std( EAKF_prior_ss, 0, 2),1) );
post_std_EAKF   = squeeze( mean(std( EAKF_post_ss , 0, 2),1) );

prior_rmse_EAKF  = zeros(n_t, 1);
post_rmse_EAKF   = zeros(n_t, 1);

for i=1:n_t
    prior_rmse_EAKF(i)  = sqrt( mean( (prior_mean_EAKF(:,i) - truth_ss(:,i)).^2 ) );
    post_rmse_EAKF(i)   = sqrt( mean( (post_mean_EAKF(:,i) - truth_ss(:,i)).^2 ) );
end

% also do the rank histogram for the prior
rk_EAKF = zeros(n_obs, n_t);

for i=1:n_t
    for j=1:n_obs
        [vals, rank]  = sort([truth_ss(j,i) squeeze(EAKF_prior_ss(j,:,i))],'ascend');
        rk_EAKF(j,i) = find(rank==1)-1; % rank is between 0 and n_ens
    end
end


exp_name='slp2_1day';
da_config='EAKF_cutoff_020';

fn_output = ['obs_loc_diag_',exp_name,'_',da_config,'.mat']
save(fn_output, 'prior_rmse_EAKF', 'post_rmse_EAKF', 'prior_std_EAKF', 'post_std_EAKF','rk_EAKF')
%}


