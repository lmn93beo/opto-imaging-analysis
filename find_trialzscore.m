function [trials, trials_zscore] = find_trialzscore(data, dff, frametime, avgFR, options)
dt = options.dt;
correct_artifact = options.correct_artifact;


ts = data.trial_start';
%ts_frame = fix((ts * avgFR)+1); %trial start time in frames

dum = repmat(frametime,numel(ts),1);
[~,ix] = min(abs(dum - ts),[],2);

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