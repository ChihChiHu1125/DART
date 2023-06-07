% state space diagnostics (on linux)

%% Read data:

clear all;

exp_name='expslp5_1day_1yr';
da_config='PFFli_ker_08cap_infR3eakffg_ensavgHT';

fn_prior = ['preassim_',exp_name,'_',da_config,'.nc'];
fn_post  = ['analysis_',exp_name,'_',da_config,'.nc'];
fn_true  = '/home/chihchi/scratch/Bgrid_project_NEW/share/observation_files/true_state_long.nc';

ps_prior = squeeze( ncread(fn_prior, 'ps') );
ps_post  = squeeze( ncread(fn_post , 'ps') );
ps_true  = ncread(fn_true , 'ps');

t_prior = squeeze( ncread(fn_prior, 't') );
t_post  = squeeze( ncread(fn_post , 't') );
t_true = ncread(fn_true, 't');

u_prior = squeeze( ncread(fn_prior, 'u') );
u_post  = squeeze( ncread(fn_post , 'u') );
u_true = ncread(fn_true, 'u');

v_prior = squeeze( ncread(fn_prior, 'v') );
v_post  = squeeze( ncread(fn_post , 'v') );
v_true = ncread(fn_true, 'v');
%}

%% Analysis of RMSE and spread
ns = 81; % the start time 
nf =  4; % obs frequency
[nx, ny, nz, ne, nt] = size(u_prior);

ps_prior_mean = mean(ps_prior,3);
ps_post_mean  = mean(ps_post ,3);

t_prior_mean = squeeze( mean(t_prior,4));
t_post_mean  = squeeze( mean(t_post ,4));

u_prior_mean = squeeze( mean(u_prior,4));
u_post_mean  = squeeze( mean(u_post ,4));

v_prior_mean = squeeze( mean(v_prior,4));
v_post_mean  = squeeze( mean(v_post ,4));

% mass range
range_lat = [1:30];
range_lon = [1:60];

ps_rmse_prior = squeeze( sqrt( mean(mean( (ps_prior_mean(range_lon, range_lat,:,:) - ps_true(range_lon, range_lat,:,ns:nf:ns+nf*(nt-1) )).^2 ,1),2) ) );
ps_rmse_post  = squeeze( sqrt( mean(mean( (ps_post_mean(range_lon, range_lat,:,:)  - ps_true(range_lon, range_lat,:,ns:nf:ns+nf*(nt-1) )).^2 ,1),2) ) );

ps_spread_prior = squeeze( mean(mean(squeeze( std(ps_prior,0,3) ),1),2) );
ps_spread_post  = squeeze( mean(mean(squeeze( std(ps_post ,0,3) ),1),2) );

t_rmse_prior = squeeze( sqrt( mean(mean(mean( (t_prior_mean(range_lon, range_lat,:,:,:) - t_true(range_lon, range_lat,:,ns:nf:ns+nf*(nt-1) )).^2 ,1),2),3) ) );
t_rmse_post  = squeeze( sqrt( mean(mean(mean( (t_post_mean(range_lon, range_lat,:,:,:)  - t_true(range_lon, range_lat,:,ns:nf:ns+nf*(nt-1) )).^2 ,1),2),3) ) );

t_spread_prior = squeeze( mean(mean(mean(squeeze( std(t_prior,0,4) ),1),2),3) );
t_spread_post  = squeeze( mean(mean(mean(squeeze( std(t_post ,0,4) ),1),2),3) );

% wind range
range_lat = [1:29];
range_lon = [1:60];

u_rmse_prior = squeeze( sqrt( mean(mean(mean( (u_prior_mean(range_lon, range_lat,:,:,:) - u_true(range_lon, range_lat,:,ns:nf:ns+nf*(nt-1) )).^2 ,1),2),3) ) );
u_rmse_post  = squeeze( sqrt( mean(mean(mean( (u_post_mean(range_lon, range_lat,:,:,:)  - u_true(range_lon, range_lat,:,ns:nf:ns+nf*(nt-1) )).^2 ,1),2),3) ) );

u_spread_prior = squeeze( mean(mean(mean(squeeze( std(u_prior,0,4) ),1),2),3) );
u_spread_post  = squeeze( mean(mean(mean(squeeze( std(u_post ,0,4) ),1),2),3) );

v_rmse_prior = squeeze( sqrt( mean(mean(mean( (v_prior_mean(range_lon, range_lat,:,:,:) - v_true(range_lon, range_lat,:,ns:nf:ns+nf*(nt-1) )).^2 ,1),2),3) ) );
v_rmse_post  = squeeze( sqrt( mean(mean(mean( (v_post_mean(range_lon, range_lat,:,:,:)  - v_true(range_lon, range_lat,:,ns:nf:ns+nf*(nt-1) )).^2 ,1),2),3) ) );

v_spread_prior = squeeze( mean(mean(mean(squeeze( std(v_prior,0,4) ),1),2),3) );
v_spread_post  = squeeze( mean(mean(mean(squeeze( std(v_post ,0,4) ),1),2),3) );

%% Analysis of individual member error:
ps_memerr_prior = zeros(nt,ne); 
ps_memerr_post  = zeros(nt,ne);

for t=1:nt
    for i=1:ne
        ps_memerr_prior(t,i) = mean(mean(abs(ps_prior(:,:,i,t) - ps_true(:,:,1,ns+nf*(t-1) ))));
        ps_memerr_post (t,i) = mean(mean(abs(ps_post (:,:,i,t) - ps_true(:,:,1,ns+nf*(t-1) ))));
    end
end

t_memerr_prior = zeros(nt,ne); 
t_memerr_post  = zeros(nt,ne);

for t=1:nt
    for i=1:ne
        t_memerr_prior(t,i) = mean(mean(mean(abs(t_prior(:,:,:,i,t) - t_true(:,:,:,1,ns+nf*(t-1) )))));
        t_memerr_post (t,i) = mean(mean(mean(abs(t_post (:,:,:,i,t) - t_true(:,:,:,1,ns+nf*(t-1) )))));
    end
end

u_memerr_prior = zeros(nt,ne); 
u_memerr_post  = zeros(nt,ne);

for t=1:nt
    for i=1:ne
        u_memerr_prior(t,i) = mean(mean(mean(abs(u_prior(:,:,:,i,t) - u_true(:,:,:,1,ns+nf*(t-1) )))));
        u_memerr_post (t,i) = mean(mean(mean(abs(u_post (:,:,:,i,t) - u_true(:,:,:,1,ns+nf*(t-1) )))));
    end
end

v_memerr_prior = zeros(nt,ne); 
v_memerr_post  = zeros(nt,ne);

for t=1:nt
    for i=1:ne
        v_memerr_prior(t,i) = mean(mean(mean(abs(v_prior(:,:,:,i,t) - v_true(:,:,:,1,ns+nf*(t-1) )))));
        v_memerr_post (t,i) = mean(mean(mean(abs(v_post (:,:,:,i,t) - v_true(:,:,:,1,ns+nf*(t-1) )))));
    end
end
%}

%% Save files:
output_name = ['state_diag_',exp_name,'_',da_config,'.mat'];
save(output_name, ...
     'ps_rmse_prior', 'ps_rmse_post', 'ps_spread_prior', 'ps_spread_post', 'ps_memerr_prior', 'ps_memerr_post', ...
     't_rmse_prior',  't_rmse_post',  't_spread_prior',  't_spread_post',  't_memerr_prior',  't_memerr_post',  ...
     'u_rmse_prior',  'u_rmse_post',  'u_spread_prior',  'u_spread_post',  'u_memerr_prior',  'u_memerr_post',  ...  
     'v_rmse_prior',  'v_rmse_post',  'v_spread_prior',  'v_spread_post',  'v_memerr_prior',  'v_memerr_post')
