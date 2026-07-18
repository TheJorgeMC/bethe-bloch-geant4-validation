%GENERATECURVESWITHUNCERTAINTY Re-generates the dE/dx dataset over the
%same energy grid as generateCurves.m, adding Delta_S_mass_MeV_cm2_g,
%S_mass_relative_uncertainty_pct, Delta_S_linear_MeV_cm, and
%S_linear_relative_uncertainty_pct columns.
%
%   Output files (in matlab/analytic/output/):
%       bethe_analytic_with_uncertainty_<material>.csv
%       bethe_analytic_with_uncertainty_all.csv
%       dEdx_with_uncertainty.png  -- two-panel figure: (left) S(T)
%       with +/-1 sigma shaded band, (right) relative uncertainty (%)
%       vs. energy.
%
%   Run this script directly (F5) from matlab/analytic/.

outputDir = fullfile(fileparts(mfilename('fullpath')), 'output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    graphics_toolkit('gnuplot');
end

energies = betheBlochEnergyGrid();
materials = betheBlochMaterials();
matKeys = {'water', 'aluminium', 'copper', 'lead'};

allTables = {};

for i = 1:numel(matKeys)
    key = matKeys{i};
    mat = materials.(key);

    S_mass = betheBlochMassStoppingPower(energies, mat);
    S_linear = betheBlochLinearStoppingPower(energies, mat);
    dS_mass = betheBlochPropagateMassUncertainty(energies, key);
    dS_linear = betheBlochPropagateLinearUncertainty(energies, key);
    vc = betheBlochValidityCheck(energies);

    T = table(repmat({mat.name}, numel(energies), 1), energies(:), S_mass(:), ...
        dS_mass(:), 100*dS_mass(:)./S_mass(:), S_linear(:), dS_linear(:), ...
        100*dS_linear(:)./S_linear(:), vc.beta(:), vc.gamma(:), vc.beta_gamma(:), ...
        vc.valid_uncorrected_bethe(:), ...
        'VariableNames', {'material', 'T_MeV', 'S_mass_MeV_cm2_g', ...
        'Delta_S_mass_MeV_cm2_g', 'S_mass_relative_uncertainty_pct', ...
        'S_linear_MeV_cm', 'Delta_S_linear_MeV_cm', 'S_linear_relative_uncertainty_pct', ...
        'beta', 'gamma', 'beta_gamma', 'valid_uncorrected_bethe'});

    outPath = fullfile(outputDir, sprintf('bethe_analytic_with_uncertainty_%s.csv', key));
    writetable(T, outPath);
    fprintf('Written %s (%d rows)\n', outPath, height(T));

    allTables{end+1} = T; %#ok<AGROW>
end

combined = vertcat(allTables{:});
combinedPath = fullfile(outputDir, 'bethe_analytic_with_uncertainty_all.csv');
writetable(combined, combinedPath);
fprintf('Written %s (%d rows)\n', combinedPath, height(combined));

% --- Two-panel figure ---
fig = figure('Position', [100, 100, 1400, 600]);
ax1 = subplot(1, 2, 1); hold(ax1, 'on');
ax2 = subplot(1, 2, 2); hold(ax2, 'on');

colors = lines(numel(matKeys));

for i = 1:numel(matKeys)
    key = matKeys{i};
    mat = materials.(key);
    sub = combined(strcmp(combined.material, mat.name), :);

    Tg = sub.T_MeV;
    S = sub.S_mass_MeV_cm2_g;
    dS = sub.Delta_S_mass_MeV_cm2_g;
    rel = sub.S_mass_relative_uncertainty_pct;

    plot(ax1, Tg, S, 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', mat.name);
    fill(ax1, [Tg; flipud(Tg)], [S-dS; flipud(S+dS)], colors(i,:), ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    plot(ax2, Tg, rel, 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', mat.name);
end

for ax = [ax1, ax2]
    xline(ax, 4.68, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    set(ax, 'XScale', 'log');
    xlim(ax, [3, 1000]);
    xlabel(ax, 'Proton Kinetic Energy, T [MeV]');
    grid(ax, 'on');
end

set(ax1, 'YScale', 'log');
ylabel(ax1, 'Mass stopping power, -(1/\rho)(dE/dx) [MeV cm^2 g^{-1}]');
title(ax1, 'Bethe (no corrections) with \pm1\sigma uncertainty band');
legend(ax1, 'show', 'Location', 'northeast');

ylabel(ax2, 'Relative uncertainty on S, \DeltaS/S [%]');
title(ax2, 'Uncertainty budget dominated by \DeltaI');
legend(ax2, 'show', 'Location', 'northeast');

figPath = fullfile(outputDir, 'dEdx_with_uncertainty.png');
saveas(fig, figPath);
fprintf('Saved %s\n', figPath);

% --- Uncertainty budget summary, mirrors Python __main__ printout ---
fprintf('\nUncertainty budget summary at T=100 MeV (variance fraction, %%):\n');
fprintf('%-12s%8s%8s%8s%8s%14s\n', 'Material', 'K', 'ZA', 'me', 'I', 'Total rel.%');
for i = 1:numel(matKeys)
    key = matKeys{i};
    mat = materials.(key);
    S = betheBlochMassStoppingPower(100.0, mat);
    [dS, b] = betheBlochPropagateMassUncertainty(100.0, key);
    totalVar = b.K^2 + b.ZA^2 + b.me^2 + b.I^2;
    fprintf('%-12s%8.3f%8.3f%8.3f%8.3f%14.3f\n', mat.name, ...
        100*b.K^2/totalVar, 100*b.ZA^2/totalVar, 100*b.me^2/totalVar, ...
        100*b.I^2/totalVar, 100*dS/S);
end
