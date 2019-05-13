function dff = get_dff_from_raw(raw)

% Find df/f for all cells, control
for c = 1:size(raw,1)
        [freq, xi] = ksdensity(raw(c,:));
        [~,idx] = max(freq);
        baseline = xi(idx);
        %baseline = prctile(data.raw_F(c,:),10);
        dff(c,:) = (raw(c,:)-baseline)/baseline*100;
%       data.DFF(c,:) = zscore(data.DFF(c,:));
end
