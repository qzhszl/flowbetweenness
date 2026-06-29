clear;
clc;

data_root = "D:\data\flow betweenness\IERP\noise";
output_dir = fileparts(mfilename('fullpath'));

N_vec = [20, 50, 100];
p = 0.4;
noise_vec = 0:0.1:1;

num_N = numel(N_vec);
num_noise = numel(noise_vec);

rgp_mean = nan(num_noise, num_N);
rgp_std = nan(num_noise, num_N);
lsm_mean = nan(num_noise, num_N);
lsm_std = nan(num_noise, num_N);

availability_N = zeros(num_N * num_noise, 1);
availability_eta = zeros(num_N * num_noise, 1);
availability_expected_files = zeros(num_N * num_noise, 1);
availability_present_files = zeros(num_N * num_noise, 1);
availability_expected_rows = zeros(num_N * num_noise, 1);
availability_loaded_rows = zeros(num_N * num_noise, 1);
availability_rgp_rows = zeros(num_N * num_noise, 1);
availability_lsm_rows = zeros(num_N * num_noise, 1);
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
    normalization_factor = nchoosek(N, 2);

    if ismember(N, [20, 50])
        inputpara_vec = 1:2;
        expected_rows_per_file = 500;
    else
        inputpara_vec = 1:20;
        expected_rows_per_file = 50;
    end

    for count_noise = 1:num_noise
        noise_amplitude = noise_vec(count_noise);
        rgp_values = zeros(0, 1);
        lsm_values = zeros(0, 1);
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

            current_rgp = results(:, 2) / normalization_factor;
            current_lsm = results(:, 6) / normalization_factor;
            rgp_values = [rgp_values; current_rgp(isfinite(current_rgp))]; %#ok<AGROW>
            lsm_values = [lsm_values; current_lsm(isfinite(current_lsm))]; %#ok<AGROW>
        end

        if ~isempty(rgp_values)
            rgp_mean(count_noise, count_N) = mean(rgp_values);
            rgp_std(count_noise, count_N) = std(rgp_values);
        end

        if ~isempty(lsm_values)
            lsm_mean(count_noise, count_N) = mean(lsm_values);
            lsm_std(count_noise, count_N) = std(lsm_values);
        end

        availability_index = availability_index + 1;
        availability_N(availability_index) = N;
        availability_eta(availability_index) = noise_amplitude;
        availability_expected_files(availability_index) = numel(inputpara_vec);
        availability_present_files(availability_index) = present_files;
        availability_expected_rows(availability_index) = ...
            numel(inputpara_vec) * expected_rows_per_file;
        availability_loaded_rows(availability_index) = loaded_rows;
        availability_rgp_rows(availability_index) = numel(rgp_values);
        availability_lsm_rows(availability_index) = numel(lsm_values);

        if isempty(rgp_values) && isempty(lsm_values)
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
    availability_loaded_rows, availability_rgp_rows, ...
    availability_lsm_rows, availability_status, ...
    'VariableNames', {'N', 'Eta', 'ExpectedFiles', 'PresentFiles', ...
    'ExpectedRows', 'LoadedRows', 'ValidRGPRows', ...
    'ValidLSMRows', 'Status'});

issue_table = table( ...
    issue_N, issue_eta, issue_inputpara, issue_type, ...
    issue_expected_rows, issue_actual_rows, issue_file, ...
    'VariableNames', {'N', 'Eta', 'Inputpara', 'Issue', ...
    'ExpectedRows', 'ActualRows', 'File'});

availability_file = fullfile(output_dir, ...
    "IERP_noise_L_data_availability.csv");
issue_file_csv = fullfile(output_dir, "IERP_noise_L_missing_data.csv");
issue_file_txt = fullfile(output_dir, "IERP_noise_L_missing_data.txt");
writetable(availability_table, availability_file);
writetable(issue_table, issue_file_csv);
write_issue_report(issue_file_txt, availability_table, issue_table);

fig = figure;
fig.Position = [100, 100, 650, 360];
hold on;

colors = ["#D08082", "#6FB494", "#D9B382"];
rgp_handles = gobjects(num_N, 1);
lsm_handles = gobjects(num_N, 1);
legend_handles = gobjects(2 * num_N, 1);
legend_labels = strings(2 * num_N, 1);

for count_N = 1:num_N
    rgp_handles(count_N) = errorbar( ...
        noise_vec, rgp_mean(:, count_N), rgp_std(:, count_N), ...
        'Color', colors(count_N), ...
        'LineStyle', '-', ...
        'LineWidth', 1.5, ...
        'Marker', 'o', ...
        'MarkerSize', 5, ...
        'MarkerFaceColor', 'white', ...
        'CapSize', 5);

    lsm_handles(count_N) = errorbar( ...
        noise_vec, lsm_mean(:, count_N), lsm_std(:, count_N), ...
        'Color', colors(count_N), ...
        'LineStyle', '--', ...
        'LineWidth', 1.5, ...
        'Marker', 's', ...
        'MarkerSize', 5, ...
        'MarkerFaceColor', colors(count_N), ...
        'CapSize', 5);

    legend_idx = 2 * count_N - 1;
    legend_handles(legend_idx) = rgp_handles(count_N);
    legend_handles(legend_idx + 1) = lsm_handles(count_N);
    legend_labels(legend_idx) = ...
        sprintf('RGP, $N=%d$', N_vec(count_N));
    legend_labels(legend_idx + 1) = ...
        sprintf('LSM, $N=%d$', N_vec(count_N));
end

xlabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$L/{N \choose 2}$', 'Interpreter', 'latex', 'FontSize', 16);
xlim([min(noise_vec) - 0.05, max(noise_vec) + 0.05]);
xticks(noise_vec);
xticklabels(compose('%.1f', noise_vec));
xtickangle(0);
box on;

lgd = legend(legend_handles, cellstr(legend_labels), ...
    'Interpreter', 'latex', ...
    'FontSize', 12, ...
    'Location', 'northeast', ...
    'NumColumns', 1, ...
    'Box', 'on');
lgd.ItemTokenSize = [12, 10];
lgd.Color = 'none';

ax = gca;
ax.FontSize = 12;
ax.YScale = 'linear';

all_upper = [rgp_mean + rgp_std; lsm_mean + lsm_std];
finite_upper = all_upper(isfinite(all_upper) & all_upper > 0);
all_lower = [rgp_mean - rgp_std; lsm_mean - lsm_std];
finite_lower = all_lower(isfinite(all_lower) & all_lower > 0);
if ~isempty(finite_upper) && ~isempty(finite_lower)
    y_min = 10 ^ floor(log10(min(finite_lower)));
    y_max = 10 ^ ceil(log10(max(finite_upper)));
    if y_min == y_max
        y_min = y_min / 10;
        y_max = y_max * 10;
    end
    ylim([y_min, y_max]);
end

pdf_file = fullfile(output_dir, "IERP_ER_noise_L.pdf");
png_file = fullfile(output_dir, "IERP_ER_noise_L.png");
fig_file = fullfile(output_dir, "IERP_ER_noise_L.fig");
exportgraphics(fig, pdf_file, ...
    'BackgroundColor', 'none', 'ContentType', 'vector');
exportgraphics(fig, png_file, ...
    'BackgroundColor', 'white', 'Resolution', 600);
savefig(fig, fig_file);

fprintf("Created plot and data reports in:\n%s\n", output_dir);
disp(availability_table);

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
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, "IERP noise L data availability report\n");
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
            'RGP valid rows=%d, LSM valid rows=%d\n'], ...
            partial_points.N(row), partial_points.Eta(row), ...
            partial_points.PresentFiles(row), ...
            partial_points.ExpectedFiles(row), ...
            partial_points.LoadedRows(row), ...
            partial_points.ExpectedRows(row), ...
            partial_points.ValidRGPRows(row), ...
            partial_points.ValidLSMRows(row));
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
