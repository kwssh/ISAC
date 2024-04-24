folder_name = '0422';
folder_velocity = '50';
folder_numuser = '1';
folder_pattern_1 = strcat(folder_name, '\');

folder_pattern_ISAC = strcat(folder_name, '_-*_', folder_velocity);
folder_pattern_commonly = strcat(folder_name, '_commonly_', folder_numuser, '_', folder_velocity);
folder_pattern_SF = strcat(folder_name, '_SF_-*_', folder_velocity);

current_path = strcat(pwd, '\figure\', folder_pattern_1);
folder_ISAC_path = strcat(current_path, folder_pattern_ISAC);
folder_commonly_path = strcat(current_path, folder_pattern_commonly);
folder_SF_path = strcat(current_path, folder_pattern_SF);

folders_ISAC = dir(folder_ISAC_path); 
folder_ISAC_names = {folders_ISAC.name};

folders_SF = dir(folder_SF_path); 
folder_SF_names = {folders_SF.name};

pattern = '-?\d+';

x_ISAC = zeros(1, length(folder_ISAC_names));
x_SF = zeros(1, length(folder_ISAC_names));

y_ISAC = zeros(1, length(folder_ISAC_names));
y_commonly = zeros(1, length(folder_ISAC_names));
y_SF = zeros(1, length(folder_ISAC_names));

for i = 1:length(folder_ISAC_names)
    numbers_ISAC = regexp(folder_ISAC_names{i}, pattern, 'match');
    x_ISAC(i) = str2double(numbers_ISAC{2});

    numbers_SF = regexp(folder_SF_names{i}, pattern, 'match');
    x_SF(i) = str2double(numbers_SF{2});
end

folder_ISAC_paths = cellfun(@(x) fullfile(current_path, x), folder_ISAC_names, 'UniformOutput', false);
folder_SF_paths = cellfun(@(x) fullfile(current_path, x), folder_SF_names, 'UniformOutput', false);

for i = 1:length(folder_ISAC_paths)

    cd(folder_ISAC_paths{i});

    line_ISAC = readlines('result.txt');

    numbers = regexp(line_ISAC(end-1), '\d+\.\d+', 'match');
    rate = str2double(numbers{1});
    y_ISAC(i) = rate;

    cd('..');
    cd('..');
    cd('..');
end

[x_ISAC_sorted, idx] = sort(x_ISAC);
y_ISAC_sorted = y_ISAC(idx);


for i = 1:length(folder_SF_paths)

    cd(folder_SF_paths{i});

    line_SF = readlines('result.txt');

    numbers = regexp(line_SF(end-1), '\d+\.\d+', 'match');
    rate = str2double(numbers{1});
    y_SF(i) = rate;


    cd('..');
    cd('..');
    cd('..');
end

[x_SF_sorted, idx] = sort(x_SF);
y_SF_sorted = y_SF(idx);


cd(folder_commonly_path);
line_commonly = readlines('result.txt');
numbers = regexp(line_commonly(end-1), '\d+\.\d+', 'match');
rate = str2double(numbers{1});
y_commonly = repmat(rate, size(y_commonly));
cd('..');
cd('..');
cd('..');


grid on
hold on
box on
plot(x_ISAC_sorted, y_commonly, 'o-', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'red');
plot(x_ISAC_sorted, y_ISAC_sorted, 'o-', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'green');
plot(x_SF_sorted, y_SF_sorted, 'o-', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'blue');
plot(x_ISAC_sorted, y_commonly, 'o-', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'red');
xlabel('Sensing Threshold (dBm)');
ylabel('Sum Rate (bps/Hz)');
xticks(min(x_ISAC_sorted):3:max(x_ISAC_sorted));
ylim([0 rate + 3])
legend('Communication only', 'Optimizing UAV trajectory', 'Straight forward');
