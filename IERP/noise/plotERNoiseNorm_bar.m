clear;
clc;

data_root = "D:\data\flow betweenness\IERP\noise";
output_dir = fileparts(mfilename('fullpath'));

N_vec = [20, 50, 100];
p = 0.4;
noise_vec = 0:0.1:1;

num_N = numel(N_vec);
num_noise = numel(noise_vec);

ierp_mean = nan(num_noise, num_N);
ierp_std = nan(num_noise, num_N);
benchmark_mean = nan(num_noise, num_N);
benchmark_std = nan(num_noise, num_N);

availability_N = zeros(num_N * num_noise, 1);
availability_eta = zeros(num_N * num_noise, 1);
availability_expected_files = zeros(num_N * num_noise, 1);
availability_present_files = zeros(num_N * num_noise, 1);
availability_expected_rows = zeros(num_N * num_noise, 1);
availability_loaded_rows = zeros(num_N * num_noise, 1);
availability_ierp_rows = zeros(num_N * num_noise, 1);
availability_benchmark_rows = zeros(num_N * num_noise, 1);
availability_status = strings(num_N * num_noise, 1);

issue_N = zeros(0, 1);
issue_eta = zeros(0, 1);
issue_inputpara = zeros(0, 1);
issue_type = strings(0, 1);
issue_expected_rows = zeros(0, 1);
issue_actual_rows = zeros(0, 1);
issue_file = strings(0, 1);

availability_index = 0;

for count_N = 1:num_N
    N = N_vec(count_N);

    if ismember(N, [20, 50])
        inputpara_vec = 1:2;
        expected_rows_per_file = 500;
    else
        inputpara_vec = 1:50;
        expected_rows_per_file = 20;
    end

    for count_noise = 1:num_noise
        noise_amplitude = noise_vec(count_noise);
        norm_values = zeros(0, 1);
        benchmark_values = zeros(0, 1);
        present_files = 0;
        loaded_rows = 0;
        has_issue = false;

        for inputpara = inputpara_vec
            filename = sprintf( ...
                "IERP_N%dERp%.4f_noise%.2f_simu%d.txt", ...
                N, p, noise_amplitude, inputpara);
            filepath = fullfile(data_root, filename);

            if ~isfile(filepath)
                has_issue = true;
                [issue_N, issue_eta, issue_inputpara, issue_type, ...
                    issue_expected_rows, issue_actual_rows, issue_file] = ...
                    append_issue(issue_N, issue_eta, issue_inputpara, ...
                    issue_type, issue_expected_rows, issue_actual_rows, ...
                    issue_file, N, noise_amplitude, inputpara, ...
                    "MissingFile", expected_rows_per_file, 0, filepath);
                continue;
            end

            try
                results = readmatrix(filepath);
            catch read_error
                has_issue = true;
                [issue_N, issue_eta, issue_inputpara, issue_type, ...
                    issue_expected_rows, issue_actual_rows, issue_file] = ...
                    append_issue(issue_N, issue_eta, issue_inputpara, ...
                    issue_type, issue_expected_rows, issue_actual_rows, ...
                    issue_file, N, noise_amplitude, inputpara, ...
                    "ReadError: " + string(read_error.identifier), ...
                    expected_rows_per_file, 0, filepath);
                continue;
            end

            present_files = present_files + 1;
            actual_rows = size(results, 1);
            loaded_rows = loaded_rows + actual_rows;

            if actual_rows ~= expected_rows_per_file
                has_issue = true;
                [issue_N, issue_eta, issue_inputpara, issue_type, ...
                    issue_expected_rows, issue_actual_rows, issue_file] = ...
                    append_issue(issue_N, issue_eta, issue_inputpara, ...
                    issue_type, issue_expected_rows, issue_actual_rows, ...
                    issue_file, N, noise_amplitude, inputpara, ...
                    "UnexpectedRowCount", expected_rows_per_file, ...
                    actual_rows, filepath);
            end

            if size(results, 2) < 8
                has_issue = true;
                [issue_N, issue_eta, issue_inputpara, issue_type, ...
                    issue_expected_rows, issue_actual_rows, issue_file] = ...
                    append_issue(issue_N, issue_eta, issue_inputpara, ...
                    issue_type, issue_expected_rows, issue_actual_rows, ...
                    issue_file, N, noise_amplitude, inputpara, ...
                    "FewerThanEightColumns", expected_rows_per_file, ...
                    actual_rows, filepath);
                continue;
            end

            current_ierp = results(:, 4);
            current_benchmark = results(:, 8);
            norm_values = [norm_values; current_ierp(isfinite(current_ierp))]; %#ok<AGROW>
            benchmark_values = [benchmark_values; ...
                current_benchmark(isfinite(current_benchmark))]; %#ok<AGROW>
        end

        if ~isempty(norm_values)
            ierp_mean(count_noise, count_N) = mean(norm_values);
            ierp_std(count_noise, count_N) = std(norm_values);
        end

        if ~isempty(benchmark_values)
            benchmark_mean(count_noise, count_N) = mean(benchmark_values);
            benchmark_std(count_noise, count_N) = std(benchmark_values);
        end

        availability_index = availability_index + 1;
        availability_N(availability_index) = N;
        availability_eta(availability_index) = noise_amplitude;
        availability_expected_files(availability_index) = numel(inputpara_vec);
        availability_present_files(availability_index) = present_files;
        availability_expected_rows(availability_index) = ...
            numel(inputpara_vec) * expected_rows_per_file;
        availability_loaded_rows(availability_index) = loaded_rows;
        availability_ierp_rows(availability_index) = numel(norm_values);
        availability_benchmark_rows(availability_index) = ...
            numel(benchmark_values);

        if isempty(norm_values) && isempty(benchmark_values)
            availability_status(availability_index) = "MissingDataPoint";
        elseif has_issue
            availability_status(availability_index) = "PartialDataPoint";
        else
            availability_status(availability_index) = "CompleteDataPoint";
        end
    end
end

availability_table = table( ...
    availability_N, availability_eta, availability_expected_files, ...
    availability_present_files, availability_expected_rows, ...
    availability_loaded_rows, availability_ierp_rows, ...
    availability_benchmark_rows, availability_status, ...
    'VariableNames', {'N', 'Eta', 'ExpectedFiles', 'PresentFiles', ...
    'ExpectedRows', 'LoadedRows', 'ValidIERPRows', ...
    'ValidBenchmarkRows', 'Status'});

issue_table = table( ...
    issue_N, issue_eta, issue_inputpara, issue_type, ...
    issue_expected_rows, issue_actual_rows, issue_file, ...
    'VariableNames', {'N', 'Eta', 'Inputpara', 'Issue', ...
    'ExpectedRows', 'ActualRows', 'File'});

availability_file = fullfile(output_dir, ...
    "IERP_noise_data_availability.csv");
issue_file_csv = fullfile(output_dir, "IERP_noise_missing_data.csv");
issue_file_txt = fullfile(output_dir, "IERP_noise_missing_data.txt");
writetable(availability_table, availability_file);
writetable(issue_table, issue_file_csv);
write_issue_report(issue_file_txt, availability_table, issue_table);

fig = figure;
fig.Position = [100, 100, 650, 360];
hold on;

colors = ["#D08082", "#6FB494", "#D9B382"];
y_floor = 1e-3;

bar_values = nan(num_noise, 2 * num_N);
bar_errors = nan(num_noise, 2 * num_N);
bar_labels = strings(1, 2 * num_N);
for count_N = 1:num_N
    ierp_column = 2 * count_N - 1;
    benchmark_column = 2 * count_N;

    bar_values(:, ierp_column) = ierp_mean(:, count_N);
    bar_values(:, benchmark_column) = benchmark_mean(:, count_N);
    bar_errors(:, ierp_column) = ierp_std(:, count_N);
    bar_errors(:, benchmark_column) = benchmark_std(:, count_N);

    bar_labels(ierp_column) = sprintf('RGP, $N=%d$', N_vec(count_N));
    bar_labels(benchmark_column) = sprintf( ...
        'LSM, $N=%d$', N_vec(count_N));
end

bar_plot_values = bar_values;
valid_values = isfinite(bar_plot_values);
bar_plot_values(valid_values) = max(bar_plot_values(valid_values), y_floor);

bar_handles = bar(noise_vec, bar_plot_values, 'grouped');
for count_N = 1:num_N
    ierp_column = 2 * count_N - 1;
    benchmark_column = 2 * count_N;

    bar_handles(ierp_column).FaceColor = colors(count_N);
    bar_handles(ierp_column).FaceAlpha = 0.85;
    bar_handles(ierp_column).EdgeColor = colors(count_N);
    bar_handles(ierp_column).BaseValue = y_floor;

    bar_handles(benchmark_column).FaceColor = colors(count_N);
    bar_handles(benchmark_column).FaceAlpha = 0.35;
    bar_handles(benchmark_column).EdgeColor = colors(count_N);
    bar_handles(benchmark_column).LineWidth = 1.2;
    bar_handles(benchmark_column).BaseValue = y_floor;
end

for bar_index = 1:numel(bar_handles)
    x_positions = bar_handles(bar_index).XEndPoints;
    y_positions = bar_plot_values(:, bar_index);
    upper_errors = bar_errors(:, bar_index);
    valid_error = isfinite(x_positions(:)) & isfinite(y_positions) & ...
        isfinite(upper_errors);

    errorbar(x_positions(valid_error), y_positions(valid_error), ...
        zeros(sum(valid_error), 1), upper_errors(valid_error), ...
        'k', ...
        'LineStyle', 'none', ...
        'LineWidth', 1, ...
        'CapSize', 5, ...
        'HandleVisibility', 'off');
end

bar_upper_values = bar_plot_values + max(bar_errors, 0);
finite_bar_values = bar_upper_values(isfinite(bar_upper_values));
if isempty(finite_bar_values)
    y_max = 1;
else
    y_max = max(finite_bar_values);
end

xlabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$|D-\Omega|$', ...
    'Interpreter', 'latex', 'FontSize', 16);
xlim([min(noise_vec) - 0.05, max(noise_vec) + 0.05]);
xticks(noise_vec);
xticklabels(compose('%.1f', noise_vec));
xtickangle(0);
box on;

lgd = legend(bar_handles, cellstr(bar_labels), ...
    'Interpreter', 'latex', ...
    'FontSize', 14, ...
    'Location', 'northwest', ...
    'NumColumns', 2, ...
    'Box', 'on');
lgd.Position(1) = lgd.Position(1) - 0.04; % 向左
lgd.Position(2) = lgd.Position(2) + 0.015; % 向上

lgd.ItemTokenSize = [12, 10];
lgd.Color = 'none';

ax = gca;
ax.FontSize = 12;
ax.YScale = 'log';
ylim([y_floor, 10 ^ ceil(log10(max(y_max, 1)))]);

pdf_file = fullfile(output_dir, "IERP_ER_noise_norm.pdf");
png_file = fullfile(output_dir, "IERP_ER_noise_norm.png");
fig_file = fullfile(output_dir, "IERP_ER_noise_norm.fig");
exportgraphics(fig, pdf_file, ...
    'BackgroundColor', 'none', 'ContentType', 'vector');
exportgraphics(fig, png_file, ...
    'BackgroundColor', 'white', 'Resolution', 600);
savefig(fig, fig_file);

fprintf("Created plot and data reports in:\n%s\n", output_dir);
disp(availability_table);

disp("mfilename fullpath =")
disp(mfilename('fullpath'))

disp("output_dir =")
disp(output_dir)

disp("pdf_file =")
disp(pdf_file)
disp("png_file =")
disp(png_file)
disp("fig_file =")
disp(fig_file)


function [issue_N, issue_eta, issue_inputpara, issue_type, ...
    issue_expected_rows, issue_actual_rows, issue_file] = ...
    append_issue(issue_N, issue_eta, issue_inputpara, issue_type, ...
    issue_expected_rows, issue_actual_rows, issue_file, N, eta, ...
    inputpara, type, expected_rows, actual_rows, filepath)

issue_N(end + 1, 1) = N;
issue_eta(end + 1, 1) = eta;
issue_inputpara(end + 1, 1) = inputpara;
issue_type(end + 1, 1) = type;
issue_expected_rows(end + 1, 1) = expected_rows;
issue_actual_rows(end + 1, 1) = actual_rows;
issue_file(end + 1, 1) = filepath;
end

function write_issue_report(filepath, availability_table, issue_table)
fid = fopen(filepath, 'w');
if fid == -1
    error("Unable to create issue report: %s", filepath);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, "IERP noise data availability report\n");
fprintf(fid, "Generated: %s\n\n", string(datetime('now')));

fprintf(fid, "Missing data points (no valid rows loaded):\n");
missing_points = availability_table( ...
    availability_table.Status == "MissingDataPoint", :);
if isempty(missing_points)
    fprintf(fid, "None\n");
else
    for row = 1:height(missing_points)
        fprintf(fid, "N=%d, eta=%.1f\n", ...
            missing_points.N(row), missing_points.Eta(row));
    end
end

fprintf(fid, "\nPartial data points:\n");
partial_points = availability_table( ...
    availability_table.Status == "PartialDataPoint", :);
if isempty(partial_points)
    fprintf(fid, "None\n");
else
    for row = 1:height(partial_points)
        fprintf(fid, ...
            ['N=%d, eta=%.1f: files %d/%d, rows %d/%d, ' ...
            'IERP valid rows=%d, benchmark valid rows=%d\n'], ...
            partial_points.N(row), partial_points.Eta(row), ...
            partial_points.PresentFiles(row), ...
            partial_points.ExpectedFiles(row), ...
            partial_points.LoadedRows(row), ...
            partial_points.ExpectedRows(row), ...
            partial_points.ValidIERPRows(row), ...
            partial_points.ValidBenchmarkRows(row));
    end
end

fprintf(fid, "\nFile-level issues:\n");
if isempty(issue_table)
    fprintf(fid, "None\n");
else
    for row = 1:height(issue_table)
        fprintf(fid, ...
            ['N=%d, eta=%.1f, inputpara=%d, issue=%s, ' ...
            'rows=%d/%d, file=%s\n'], ...
            issue_table.N(row), issue_table.Eta(row), ...
            issue_table.Inputpara(row), issue_table.Issue(row), ...
            issue_table.ActualRows(row), issue_table.ExpectedRows(row), ...
            issue_table.File(row));
    end
end
end
