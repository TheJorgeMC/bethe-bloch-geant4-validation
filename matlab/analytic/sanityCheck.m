%SANITYCHECK Proton de 100 MeV en agua -- comparacion contra NIST PSTAR.
%
%   NOTE: the expected value changed from earlier drafts of this project.
%   With WATER.I_eV = 78 eV (ICRU 90, current), S is ~7.2545 MeV cm^2/g,
%   not 7.2909 (which was the value obtained with the older ICRU 37/49
%   I=75 eV). The Python companion notebook was executed and gives
%   exactly S = 7.2545 MeV cm^2/g at T=100 MeV in water -- this script
%   should reproduce that number to floating-point precision.

materials = betheBlochMaterials();
water = materials.water;

T_test = 100.0;
[beta, gamma] = betheBlochBetaGamma(T_test, []);
S = betheBlochMassStoppingPower(T_test, water);

fprintf('Proton %.1f MeV in water:\n', T_test);
fprintf('    beta = %.5f, gamma = %.5f,     beta*gamma = %.4f\n', beta, gamma, beta*gamma);
fprintf('    S (no corrections) = %.4f [MeV cm^2/g]\n', S);
fprintf('NIST PSTAR Reference Value (with corrections) ~ 7.289,\n');
fprintf('see ''pstar_water.csv'' file for comparison with NIST PSTAR values.\n');
diff_pct = abs(S - 7.289) / 7.289 * 100;
fprintf('Difference from NIST PSTAR value: %.2f%%\n', diff_pct);
if diff_pct < 2
    fprintf('Within 2%% of NIST PSTAR value.\n');
else
    fprintf('Over 2%% of NIST PSTAR value.\n');
end

expected_python = 7.2545;
rel_diff = abs(S - expected_python) / expected_python;
fprintf('\nDifference vs. Python result (%.4f): %.6f%%\n', expected_python, rel_diff*100);
if rel_diff < 1e-4
    fprintf('-> OK: MATLAB matches Python within floating-point error.\n');
else
    fprintf('-> WARNING: check Python/MATLAB cross-validation (Fase 3).\n');
end
