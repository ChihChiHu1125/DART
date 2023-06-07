% state space diagnostics new (includes rank histogram) (on linux)
% 2023/05/23

%% Read data:

clear all;

exp_name='slp_1day_1yr';
da_config='PFFli_ker_01_ni_30';


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

ps_spread_prior = squeeze( sqrt( mean(mean(squeeze( var(ps_prior,0,3) ),1),2) ) );
ps_spread_post  = squeeze( sqrt( mean(mean(squeeze( var(ps_post ,0,3) ),1),2) ) );

t_rmse_prior = squeeze( sqrt( mean(mean(mean( (t_prior_mean(range_lon, range_lat,:,:,:) - t_true(range_lon, range_lat,:,ns:nf:ns+nf*(nt-1) )).^2 ,1),2),3) ) );
t_rmse_post  = squeeze( sqrt( mean(mean(mean( (t_post_mean(range_lon, range_lat,:,:,:)  - t_true(range_lon, range_lat,:,ns:nf:ns+nf*(nt-1) )).^2 ,1),2),3) ) );

t_spread_prior = squeeze( sqrt( mean(mean(mean(squeeze( var(t_prior,0,4) ),1),2),3) ) );
t_spread_post  = squeeze( sqrt( mean(mean(mean(squeeze( var(t_post ,0,4) ),1),2),3) ) );

% wind range
range_lat = [1:29];
range_lon = [1:60];

u_rmse_prior = squeeze( sqrt( mean(mean(mean( (u_prior_mean(range_lon, range_lat,:,:,:) - u_true(range_lon, range_lat,:,ns:nf:ns+nf*(nt-1) )).^2 ,1),2),3) ) );
u_rmse_post  = squeeze( sqrt( mean(mean(mean( (u_post_mean(range_lon, range_lat,:,:,:)  - u_true(range_lon, range_lat,:,ns:nf:ns+nf*(nt-1) )).^2 ,1),2),3) ) );

u_spread_prior = squeeze( sqrt( mean(mean(mean(squeeze( var(u_prior,0,4) ),1),2),3) ) );
u_spread_post  = squeeze( sqrt( mean(mean(mean(squeeze( var(u_post ,0,4) ),1),2),3) ) );

v_rmse_prior = squeeze( sqrt( mean(mean(mean( (v_prior_mean(range_lon, range_lat,:,:,:) - v_true(range_lon, range_lat,:,ns:nf:ns+nf*(nt-1) )).^2 ,1),2),3) ) );
v_rmse_post  = squeeze( sqrt( mean(mean(mean( (v_post_mean(range_lon, range_lat,:,:,:)  - v_true(range_lon, range_lat,:,ns:nf:ns+nf*(nt-1) )).^2 ,1),2),3) ) );

v_spread_prior = squeeze( sqrt( mean(mean(mean(squeeze( var(v_prior,0,4) ),1),2),3) ) );
v_spread_post  = squeeze( sqrt( mean(mean(mean(squeeze( var(v_post ,0,4) ),1),2),3) ) );
%}

%% Analysis of the rank histogram

ns = 81; % the start time 
nf =  4; % obs frequency

[nx, ny, nz, ne, nt] = size(t_prior);
ps_rk_hist = zeros(nx, ny, nt);
t_rk_hist  = zeros(nx, ny, nz, nt);

for i=1:nt
    time = ns + nf*(i-1);
    for j=1:ny
        for k=1:nx
            [vals, rank]  = sort([squeeze(ps_true(k,j,:,time)) squeeze(ps_prior(k,j,:,i))'],'ascend');
            ps_rk_hist(k,j,i) = find(rank==1)-1; % rank is between 0 and n_ens
        end
    end
    
    for j=1:nz
        for k=1:ny
            for p=1:nx
                [vals, rank]  = sort([squeeze(t_true(p,k,j,:,time)) squeeze(t_prior(p,k,j,:,i))'],'ascend');
                t_rk_hist(p,k,j,i) = find(rank==1)-1; % rank is between 0 and n_ens
            end
        end
    end
end

[nx, ny, nz, ne, nt] = size(u_prior);
u_rk_hist = zeros(nx, ny, nz, nt);
v_rk_hist = zeros(nx, ny, nz, nt);

for i=1:nt
    time = ns + nf*(i-1);
   
    for j=1:nz
        for k=1:ny
            for p=1:nx
                [vals, rank]  = sort([squeeze(u_true(p,k,j,:,time)) squeeze(u_prior(p,k,j,:,i))'],'ascend');
                u_rk_hist(p,k,j,i) = find(rank==1)-1; % rank is between 0 and n_ens
                
                [vals, rank]  = sort([squeeze(v_true(p,k,j,:,time)) squeeze(v_prior(p,k,j,:,i))'],'ascend');
                v_rk_hist(p,k,j,i) = find(rank==1)-1; % rank is between 0 and n_ens
            end
        end
    end
    
end
%}


%% Save files:

output_name = ['state_diag_',exp_name,'_',da_config,'.mat'];
save(output_name, ...
     'ps_rmse_prior', 'ps_rmse_post', 'ps_spread_prior', 'ps_spread_post', 'ps_rk_hist', ...
     't_rmse_prior',  't_rmse_post',  't_spread_prior',  't_spread_post',  't_rk_hist' , ...
     'u_rmse_prior',  'u_rmse_post',  'u_spread_prior',  'u_spread_post',  'u_rk_hist' , ...  
     'v_rmse_prior',  'v_rmse_post',  'v_spread_prior',  'v_spread_post',  'v_rk_hist'     )
%}
