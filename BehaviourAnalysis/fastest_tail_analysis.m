%% Looming response time analysis: fastest-tail (q10), bootstrap CI,
%  critical threshold T* (first significant)

clear; clc;

%% Data (ms)
rt_leg = [21.6 17.6 54 20.8 41.2 43.6 18.8 59.2 88 38 16.8 57.2 77.6 ...
      22.4 127.2 44.4 59.6 114.8 76 25.6 44.4 25.6 39.2 35.2 36.4 ...
      106.8 132.4 34 17 32 127 16 18 42 48]';

rt_antena = [19 23 24 23 26 19 24 18 24 25 25 25 22 24 26 26 27 28 31 25 ...
      28 27 26 27 27 24 23 32 28 23 26 28 26 31 30 29 25 26 28 29 ...
      31 33 25 27 39 24 26 26 26 24 24 28 22 24 23 22 22 21 28 30 ...
      16 30 27 27 30 27 24 22 31 24 26 23 30 28 26 28 14 26 27]';

rt = rt_leg;

n = numel(rt);

% Parameters
q = 10;              % 10th percentile = fastest 10%
B = 50000;           % bootstrap iterations
alpha = 0.05;        % significance level
rng(0);              % reproducibility

%%  Point estimates
q10    = prctile(rt, q);
rt_min = min(rt);
rt_med = median(rt);

%% Bootstrap distribution of q10
boot_q10 = zeros(B,1);
for b = 1:B
    s = rt(randsample(n, n, true));  % resample trials with replacement
    boot_q10(b) = prctile(s, q);
end

% Two-sided 95% bootstrap CI for q10
ci_q10 = prctile(boot_q10, [100*alpha/2, 100*(1-alpha/2)]);

% "First significant" threshold T*
Tcrit    = prctile(boot_q10, 100*(1-alpha));  
T_report = ceil(Tcrit);                         % report conservatively

% One-sided bootstrap p-value for H0: q10 >= T_report vs H1: q10 < T_report
p_one_sided = mean(boot_q10 >= T_report);

%%  Print results
fprintf('--- Looming response time summary ---\n');
fprintf('n = %d trials\n', n);
fprintf('Minimum RT = %.1f ms\n', rt_min);
fprintf('Median  RT = %.1f ms\n', rt_med);
fprintf('q10 (10th percentile) = %.2f ms\n', q10);
fprintf('Bootstrap 95%% CI for q10: [%.2f, %.2f] ms\n', ci_q10(1), ci_q10(2));
fprintf('Critical threshold T* (alpha=%.2f): %.2f ms\n', alpha, Tcrit);
fprintf('(one-sided bootstrap test, p = %.3f)\n', p_one_sided);
fprintf('Reported threshold (rounded up): T = %d ms\n', T_report);
fprintf('One-sided bootstrap p-value for q10 < %d ms: p = %.5f\n', T_report, p_one_sided);

%%  Figure: Histogram (LEFT, blue) + ECDF (RIGHT, black with black axis)
figure('Color','w'); hold on;

ax = gca;

% LEFT axis: Histogram--bluee
yyaxis left
nbins = max(10, round(sqrt(n)) + 5);
edges = linspace(min(rt), max(rt), nbins);

h = histogram(rt, edges);     % counts
h.FaceColor = 'b';
h.EdgeColor = 'b';
h.FaceAlpha = 0.18;
h.EdgeAlpha = 0.25;
ylabel('Count');
ax.YAxis(1).Color = 'b';      % left axis colour matches histogram

% RIGHT Axis: ECDF --black
yyaxis right
ylim([0 1]);
yl = ylim;

% CI band (black tint), put behind ECDF by drawing before the line
patch([ci_q10(1) ci_q10(2) ci_q10(2) ci_q10(1)], [yl(1) yl(1) yl(2) yl(2)], ...
      [0 0 0], 'FaceAlpha', 0.08, 'EdgeColor', 'none');

[f,x] = ecdf(rt);
pECDF = plot(x, f, 'LineWidth', 2);
pECDF.Color = 'k';
ylabel('Empirical CDF');

% Reference lines (black)
xline(q10, '--r', 'LineWidth', 1.8, ...
    'Label', sprintf('q10 = %.2f ms', q10), ...
    'LabelVerticalAlignment','bottom');

xline(T_report, ':r', 'LineWidth', 2, ...
    'Label', sprintf('T* = %d ms', T_report), ...
    'LabelVerticalAlignment','bottom');

% Make the entire right axis styling black (ticks, axis line)
ax.YAxis(2).Color = 'k';

% Cosmetics
xlabel('Response time (ms)');
grid on; box on;

