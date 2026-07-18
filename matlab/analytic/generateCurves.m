%GENERATECURVES Generates dE/dx vs. Energy curves (relativistic Bethe,
%without corrections) for the four project materials, covering the
%clinical range of 3-300 MeV, with an extended margin up to 1000 MeV
%for visualization purposes only (outside the clinical scope).
%
%   Output files (in matlab/analytic/output/):
%       bethe_analytic_<material>.csv
%       bethe_analytic_all.csv
%       dEdx_vs_energy.png
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

figure('Position', [100, 100, 800, 600]);
hold on;

for i = 1:numel(matKeys)
    key = matKeys{i};
    mat = materials.(key);

    S_mass = betheBlochMassStoppingPower(energies, mat);
    S_linear = betheBlochLinearStoppingPower(energies, mat);
    vc = betheBlochValidityCheck(energies);

    T = table(repmat({mat.name}, numel(energies), 1), energies(:), S_mass(:), ...
        S_linear(:), vc.beta(:), vc.gamma(:), vc.beta_gamma(:), vc.valid_uncorrected_bethe(:), ...
        'VariableNames', {'material', 'T_MeV', 'S_mass_MeV_cm2_g', ...
        'S_linear_MeV_cm', 'beta', 'gamma', 'beta_gamma', 'valid_uncorrected_bethe'});

    outPath = fullfile(outputDir, sprintf('bethe_analytic_%s.csv', key));
    writetable(T, outPath);
    fprintf('Written %s (%d rows)\n', outPath, height(T));

    allTables{end+1} = T; %#ok<AGROW>

    loglog(energies, S_mass, 'LineWidth', 2, 'DisplayName', mat.name);
end

combined = vertcat(allTables{:});
combinedPath = fullfile(outputDir, 'bethe_analytic_all.csv');
writetable(combined, combinedPath);
fprintf('Written %s (%d rows)\n', combinedPath, height(combined));

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Proton Kinetic Energy, T [MeV]');
ylabel('Mass stopping power, -(1/\rho)(dE/dx) [MeV cm^2 g^{-1}]');
title('Relativistic Bethe (no corrections) - 4 Materials');
xlim([3, 1000]);
legend('show', 'Location', 'northeast');
grid on;
hold off;

figPath = fullfile(outputDir, 'dEdx_vs_energy.png');
saveas(gcf, figPath);
fprintf('Saved %s\n', figPath);
