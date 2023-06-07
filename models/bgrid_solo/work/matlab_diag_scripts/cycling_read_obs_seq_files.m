% cycling obs_seq_filesssss
% combine all the cycling files

%% parameters:
% for each file (manually check file to modify below)
n_header = 65;  % # of lines to skip
n_obs    = 300; % # of observations
n_block  = 67;  % # of lines for each obs block
n_ens    = 25;  % # of ensemble members

% for cycling info:
t_start  = 81;   % # of 21600 seconds (e.g., t_start = 1 -> 0 day 0 sec, t_start = 5 -> 1 day 0 sec)
t_end    = 1437;  % # of obs times
t_freq   = 4;    % obs frequency ( = t_freq* 21600 second)

n_t      = (t_end - t_start)/t_freq;
%% RUN section:

truth = zeros(n_obs, n_t);
obs   = zeros(n_obs, n_t);
prior = zeros(n_obs, n_ens, n_t);
post  = zeros(n_obs, n_ens, n_t);

t_ct=1;

for tt=t_start:t_freq:t_end

time = (tt-1)*21600;
dd   = floor(time/86400);
ss   = mod(time, 86400);
%disp(['Start the ',num2str(tt),' iteration. time= ',num2str(dd),' day ',num2str(ss),' sec.'])    
    
fn_prior = ['obs_seq.prior_',num2str(dd),'_',num2str(ss)]
fn_post  = ['obs_seq.post_',num2str(dd),'_',num2str(ss)];

fid_prior = fopen(fn_prior);
fid_post  = fopen(fn_post);

% skip the headers info
for i=1:n_header
tline = fgetl(fid_prior);
tline = fgetl(fid_post);
end

for i=1:n_obs
    tline = fgetl(fid_prior); % header
    tline = fgetl(fid_post); % header
%     disp(tline)
    
    obs(i,t_ct)   = str2double(fgetl(fid_prior));
    truth(i,t_ct) = str2double(fgetl(fid_prior));
    
    tline = fgetl(fid_post); % obs
    tline = fgetl(fid_post); % truth
   
    for j=1:4
       tline = fgetl(fid_prior); % lines to skip
       tline = fgetl(fid_post ); % lines to skip
    end
    
    for j=1:n_ens
       tline = fgetl(fid_prior); % lines to skip
       prior(i,j,t_ct) = str2double(fgetl(fid_prior));
       
       tline = fgetl(fid_post); % lines to skip
       post(i,j,t_ct)  = str2double(fgetl(fid_post));
    end
    
    for j=1:10
       tline = fgetl(fid_prior); % lines to skip       
       tline = fgetl(fid_post ); % lines to skip
    end
end

fclose(fid_prior);
fclose(fid_post);

t_ct = t_ct + 1;

end
%}

%% save output file into mat

 output_name = 'obs_seq_tmp.mat';

PFF_prior = prior;
PFF_post  = post;
%EAKF_prior = prior;
%EAKF_post  = post;
save(output_name,'truth','obs','PFF_prior','PFF_post')

%save(output_name,'truth','obs','PFF_prior','PFF_post','prior','post')

%}
