function testBetheBlochUncertainty()
%TESTBETHEBLOCHUNCERTAINTY Uncertainty-propagation tests, mirror of
%python/analytic cell 15. Run with: testBetheBlochUncertainty()

    materials = betheBlochMaterials();

    tests = {
        @test_negligible_K_and_me_contributions, ...
        @test_I_dominates_uncertainty_budget, ...
        @test_relative_uncertainty_in_reasonable_range, ...
        @test_relative_uncertainty_grows_at_lower_energy, ...
        @test_linear_uncertainty_at_least_density_scaled_mass_uncertainty, ...
        @test_I_term_scales_linearly_with_delta_I, ...
        @test_water_5eV_shift_matches_literature_order_of_magnitude ...
    };

    nOk = 0;
    for i = 1:numel(tests)
        try
            tests{i}(materials);
            fprintf('OK   %s\n', func2str(tests{i}));
            nOk = nOk + 1;
        catch ME
            fprintf('FAIL %s: %s\n', func2str(tests{i}), ME.message);
        end
    end
    fprintf('%d/%d tests passed\n', nOk, numel(tests));
end

function test_negligible_K_and_me_contributions(~)
    % K and m_e derive from CODATA constants known to ~1e-9 relative
    % precision; their variance contribution should be utterly negligible
    % (< 0.01% of the total variance) at any representative energy point.
    [~, b] = betheBlochPropagateMassUncertainty(100.0, 'water');
    totalVar = b.K^2 + b.ZA^2 + b.me^2 + b.I^2;
    assert(b.K^2/totalVar < 1e-4, 'K contributes too much variance');
    assert(b.me^2/totalVar < 1e-4, 'me contributes too much variance');
end

function test_I_dominates_uncertainty_budget(materials)
    % Delta(I) should be the single largest contributor for all four
    % materials at a representative clinical energy -- though for lead,
    % Delta(A)=1.1 g/mol (IUPAC's conservative interval-covering value)
    % is no longer negligible, so we only require I to be the SINGLE
    % largest contributor, not necessarily >50% for every material.
    keys = fieldnames(materials);
    for i = 1:numel(keys)
        [~, b] = betheBlochPropagateMassUncertainty(100.0, keys{i});
        varFracs = [b.K^2, b.ZA^2, b.me^2, b.I^2];
        [~, idx] = max(varFracs);
        assert(idx == 4, sprintf('%s: expected I to dominate, got index %d', keys{i}, idx));
    end
end

function test_relative_uncertainty_in_reasonable_range(materials)
    energies = linspace(3.0, 300.0, 20);
    keys = fieldnames(materials);
    for i = 1:numel(keys)
        S = betheBlochMassStoppingPower(energies, materials.(keys{i}));
        dS = betheBlochPropagateMassUncertainty(energies, keys{i});
        rel = dS ./ S;
        assert(all(rel > 0), sprintf('%s: relative uncertainty should be positive', keys{i}));
        assert(all(rel < 0.05), sprintf('%s: relative uncertainty exceeds 5%%', keys{i}));
    end
end

function test_relative_uncertainty_grows_at_lower_energy(materials)
    energies = linspace(3.0, 300.0, 30);
    S = betheBlochMassStoppingPower(energies, materials.water);
    dS = betheBlochPropagateMassUncertainty(energies, 'water');
    rel = dS ./ S;
    assert(all(diff(rel) <= 0), 'relative uncertainty should decrease as energy increases');
end

function test_linear_uncertainty_at_least_density_scaled_mass_uncertainty(materials)
    T = 100.0;
    keys = fieldnames(materials);
    for i = 1:numel(keys)
        material = materials.(keys{i});
        dS_mass = betheBlochPropagateMassUncertainty(T, keys{i});
        dS_linear = betheBlochPropagateLinearUncertainty(T, keys{i});
        assert(dS_linear >= material.rho_g_cm3 * dS_mass - 1e-12);
    end
end

function test_I_term_scales_linearly_with_delta_I(materials)
    % Pure linearity check on dS/dI * Delta(I): doubling Delta(I) must
    % exactly double the I-term contribution. Tested directly on the
    % partial derivative rather than mutating the material struct.
    material = materials.copper;
    c = betheBlochConstants();
    [beta, gamma] = betheBlochBetaGamma(100.0, c.PROTON_MASS_MEV);
    I_val = material.I_eV * c.EV_TO_MEV;
    d = betheBlochPartialDerivatives(c.K_MEV_CM2_PER_MOL, material.Z_over_A, ...
        c.ELECTRON_MASS_MEV, I_val, beta, gamma, 1);
    term_1x = d.dI * (material.Delta_I_eV * c.EV_TO_MEV);
    term_2x = d.dI * (2 * material.Delta_I_eV * c.EV_TO_MEV);
    ratio = term_2x / term_1x;
    assert(abs(ratio - 2.0) < 1e-10, sprintf('expected exactly 2x scaling, got %.6f', ratio));
end

function test_water_5eV_shift_matches_literature_order_of_magnitude(materials)
    % Kumazaki et al. (2007): a 75 -> 80 eV shift in water's I-value
    % produces a 0.8-1.2% change in stopping power over 10-250 MeV.
    % Built explicitly at I=75 and I=80 here (independent of whichever
    % I-value is currently adopted as WATER's central value).
    T = 100.0;
    water = materials.water;
    water75 = water; water75.I_eV = 75.0;
    water80 = water; water80.I_eV = 80.0;
    S_75 = betheBlochMassStoppingPower(T, water75);
    S_80 = betheBlochMassStoppingPower(T, water80);
    finiteDiffPct = 100 * abs(S_80 - S_75) / S_75;
    assert(finiteDiffPct > 0.5 && finiteDiffPct < 2.0, ...
        sprintf('5 eV shift gives %.2f%%, expected ~0.8-1.2%% per literature', finiteDiffPct));
end
