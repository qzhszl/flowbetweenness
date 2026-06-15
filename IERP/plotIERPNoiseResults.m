clear;
clc;

N_vec = [20, 50, 100, 200];
N_vec = [20];

p = 0.4;
noise_vec = [0, 0.02, 0.04, 0.06, 0.08, 0.10, ...
             0.12, 0.14, 0.16, 0.18, 0.20];
expected_samples = 1000;

data_mean = zeros(length(noise_vec), length(N_vec));
data_std = zeros(length(noise_vec), length(N_vec));

for countN = 1:length(N_vec)
    N = N_vec(countN);

    if N == 20 || N == 50
        inputpara_vec = 1;
    else
        inputpara_vec = 1:20;
    end

    for count_noise = 1:length(noise_vec)
        noise_amplitude = noise_vec(count_noise);
        norm_values = zeros(0, 1);

        for inputpara = inputpara_vec
            filename = sprintf( ...
                "D:\\data\\flow betweenness\\IERP\\noise\\IERP_N%dERp%.4f_noise%.2f_simu%d.txt", ...
                N, p, noise_amplitude, inputpara);

            results = readmatrix(filename);
            if size(results, 2) < 4
                error("File %s contains fewer than four columns.", filename);
            end

            % The fourth column is Norm_output.
            norm_values = [norm_values; results(:, 4)]; %#ok<AGROW>
        end

        if length(norm_values) ~= expected_samples
            error(["Expected %d samples for N = %d and noise = %.2f, " ...
                   "but loaded %d samples."], ...
                  expected_samples, N, noise_amplitude, length(norm_values));
        end

        data_mean(count_noise, countN) = mean(norm_values);
        data_std(count_noise, countN) = std(norm_values);
    end
end

fig = figure;
fig.Position = [100, 100, 900, 600];
hold on;

colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", ...
          "#6FB494", "#D9B382"];

for countplot = 1:length(N_vec)
    errorbar(noise_vec, data_mean(:, countplot), data_std(:, countplot), ...
        'o-', ...
        'Color', colors(countplot), ...
        'LineWidth', 4, ...
        'MarkerSize', 10, ...
        'CapSize', 8);
end

ax = gca;
ax.FontSize = 30;
xticks(noise_vec);
xlim([noise_vec(1), noise_vec(end)]);
xlabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 40);
ylabel('$\|D-\Omega\|$', 'Interpreter', 'latex', 'FontSize', 40);

lgd = legend({'$N = 20$', '$N = 50$', '$N = 100$', '$N = 200$'}, ...
    'Interpreter', 'latex', ...
    'Location', 'best', ...
    'FontSize', 24);
lgd.NumColumns = 2;

box on;
hold off;

script_dir = fileparts(mfilename('fullpath'));
picname = fullfile(script_dir, "IERP_ER_noise_norm.pdf");
exportgraphics(fig, picname, ...
    'BackgroundColor', 'none', ...
    'Resolution', 600);
