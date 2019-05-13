function [pre, post, inh_sig, exc_sig] = analyze_trialtypes(trial_types, num_stims)

for i = 1:numel(trial_types)
    curr = trial_types{i};
    for ii = 1:numel(num_stims)
        pre_cell(ii,:) = nanmean(curr{ii}(:,20:24),2);
        post_cell(ii,:) = nanmean(curr{ii}(:,29:34),2);
        inh_sig(i,ii) = signrank(pre_cell(ii,:), post_cell(ii,:)) < 0.05;
        exc_sig(i,ii) = signrank(pre_cell(ii,:), post_cell(ii,:)) < 0.05;
    end
    pre{i} = pre_cell;
    post{i} = post_cell;
end 