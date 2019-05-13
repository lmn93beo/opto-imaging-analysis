%% New beginnings

%clear
clc

correct_artifact = 0;
%% Load Fluoresence data and determine **exact** framerate;
% Get the necessary files and folders
[fileL_flo, dirL_flo] = uigetfile('*.csv', ...
        'Select fluorescence file, laser'); %laser files

cd(dirL_flo);
[fileC_flo, dirC_flo] = uigetfile('*.csv', ...
    'Select fluorescence file, control'); % control files

% .xml files
[fileL_xml, dirL_xml] = uigetfile('*.xml', ...
    'Select xml file, laser'); % control files
[fileC_xml, dirC_xml] = uigetfile('*.xml', ...
    'Select xml file, control'); % control files

% Stim files .mat
[fileL_mat, dirL_mat] = uigetfile('*.mat', ...
    'Select mat file, laser'); % control files
[fileC_mat, dirC_mat] = uigetfile('*.mat', ...
    'Select mat  file, control'); % control files


%%
rawF = importdata('Results_control.csv');
rawF = rawF.data(:,2:end)';

for c = 1:size(rawF,1)
        [freq xi] = ksdensity(rawF(c,:));
        [~,idx] = max(freq);
        baseline = xi(idx);
        %baseline = prctile(data.raw_F(c,:),10);
        dff(c,:) = (rawF(c,:)-baseline)/baseline*100;
%       data.DFF(c,:) = zscore(data.DFF(c,:));
end

%%
cfgfiles = [dir('*.xml')];

if isempty(cfgfiles)
    disp('Could not find config file, using default framerate')
    keyboard
else
    raw_xml = fileread(cfgfiles(2).name);
end

sep = strfind(raw_xml,'<Frame relativeTime');
frameNum = numel(sep);
for i = 1:frameNum
    curr_idx = strfind(raw_xml(sep(i):sep(i)+100),'"');
    toGet = sep(i)+curr_idx(1);
    if i == 1 %because first entry for "relativetime" is 0
        frametime(i) = str2num(raw_xml(1,toGet));
    else
        frametime(i) = str2num(raw_xml(1,toGet:toGet+8));
    end
end

exactPeriod = (frametime(end) - frametime(1))/frameNum;
avgFR = 1/exactPeriod;
%% Load file with opto stim info
%cd('C:\Users\Rafiq Huda\Dropbox (MIT)\2P4 Opto Laser\Data\');
[stim_file folder] = uigetfile('*.mat','Select the file with the optogenetic stimuli data');
load([folder stim_file]);

%% Split into trials

ts = data.trial_start';
ts_frame = fix((ts * avgFR)+1); %trial start time in frames

dum = repmat(frametime,numel(ts),1);
[~,ix] = min(abs(dum - ts),[],2);

dt = [-5 15];
dt_frame = round(dt*avgFR);

% ix_range = 0:98;
ix_range = repmat(dt_frame(1):dt_frame(2)-1,numel(ix),1); 
ix_range = ix_range + ix;
ix_range(ix_range < 1) = 1;

nCells = size(dff,1);
nTrials = size(ix_range,1);

for i = 1:nCells
    for ii = 1:nTrials
        trials{i}(ii,:) = dff(i,ix_range(ii,:));
        if correct_artifact
            artifact_frames = ceil(avgFR * data.stimType(data.stimOrder(ii),1));
            trials{i}(ii,1:artifact_frames) = NaN;
        end        
    end
    temp = trials{i}(:,1:abs(round(dt(1)*avgFR)));
    temp = reshape(temp, 1,size(temp,1)*size(temp,2));
    all_mean = repmat(nanmean(temp),size(trials{i},1), size(trials{i},2));
    all_std = repmat(nanstd(temp),size(trials{i},1), size(trials{i},2));
    trials_zscore{i} = (trials{i} - all_mean)./all_std;
end
%%

% 
% trialDur_frame = fix([diff(ts) mean(diff(ts))]*avgFR); %duration for last trial unknown -- use average duration
% numcells = size(dff,1);
% 
% for i = 1:numcells;
%     for ii = 1:numel(ts_frame);
%         pad = [];
%         curr_idx = ts_frame(ii):(ts_frame(ii)-10)+trialDur_frame(ii)-10;
%         if curr_idx(end) > frameNum
%             framesOver = curr_idx(end) - frameNum;
%             pad = nan(1,framesOver);
%             curr_idx = curr_idx(1:end-framesOver);
%         end
%         trials{i}(ii,:) = [dff(i,curr_idx) pad];
%     end
% end 
% 
stimOrder = data.stimOrder;
num_stims = unique(stimOrder);

for i = 1:nCells
    curr = trials_zscore{i};
    for ii = 1:numel(num_stims);
        curr_idx = stimOrder == ii;
        trial_types{i}{ii} = curr(curr_idx,:);
    end
end

%% Plot
close all; 

for i = 1:nCells
    curr = trial_types{i};
    for ii = 1:numel(num_stims)
        toPlot = nanmean(curr{ii});
        pre = nanmean(curr{ii}(:,20:24),2);
        post = nanmean(curr{ii}(:,29:34),2);
        inh_sig(i,ii) = signrank(pre,post) & mean(post) < mean(pre) ;
        exc_sig(i,ii) = signrank(pre,post) & mean(post) > mean(pre) ;
%         ylim([min(toPlot)-0.1 max(max(toPlot))]);
        hold on;
%         if inh_sig(i,ii) == 1
%             figure(i);
%             plot_shadedaverage(curr{ii},-4.8:0.2:15,[0 0 0],'-');
%             xlim([-2 5]);
%             ylim([-0.2 2]); shg
%         end
    end
end 

%%
% inhibited = find(inh_sig)';
% % inhibited = 1:nCells;
% excited = find(exc_sig)';
% 
% for i = 1:numel(inhibited)
%     inhibited_z(i,:) = nanmean(trial_types{inhibited(i)}{1},1);
% end
% 
% for i = 1:numel(excited)
%     excited_z(i,:) = nanmean(trial_types{excited(i)}{1},1);
% end
% 
% plot_shadedaverage(inhibited_z,-4.8:0.2:15,[0 0 0],'-');
% xlim([-1 5])
% ylim([-1 0.5]);
% 
% figure;
% plot_shadedaverage(excited_z,-4.8:0.2:15,[0 0 0],'-');
% xlim([-1 5])
% ylim([-0.5 1]);
%% Look at changes in correlation structure


