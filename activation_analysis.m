data = xlsread('Results.csv');

% Eliminate first time
data = data(2:end,:);


% Plot
plot(data(:, 1));

trials = reshape(data(:, 1), [], 20);