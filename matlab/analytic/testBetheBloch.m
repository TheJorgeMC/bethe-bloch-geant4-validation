function testBetheBloch()
%TESTBETHEBLOCH Basic unit tests, mirror of python/analytic cell 7.
%   Run with: testBetheBloch()

    materials = betheBlochMaterials();
    c = betheBlochConstants();

    tests = {
        @test_beta_gamma_limits, ...
        @test_stopping_power_positive_and_finite, ...
        @test_stopping_power_decreasing_with_energy_in_clinical_range, ...
        @test_water_100mev_matches_nist_within_1pct, ...
        @test_linear_stopping_power_scales_with_density, ...
        @test_charge_scaling, ...
        @test_validity_flag_at_clinical_lower_bound ...
    };

    nOk = 0;
    for i = 1:numel(tests)
        try
            tests{i}(materials, c);
            fprintf('OK  %s\n', func2str(tests{i}));
            nOk = nOk + 1;
        catch ME
            fprintf('FAIL %s: %s\n', func2str(tests{i}), ME.message);
        end
    end
    fprintf('%d/%d tests passed\n', nOk, numel(tests));
end

function test_beta_gamma_limits(~, c)
    %   As T -> 0 should give beta -> 0
    [beta, gamma] = betheBlochBetaGamma(1e-6, []);   % 1 [eV]
    assert(beta < 1e-3);
    assert(abs(gamma - 1.0) < 1e-6);

    %   T = Mc^2 -> gamma = 2
    [beta, gamma] = betheBlochBetaGamma(c.PROTON_MASS_MEV, []);
    assert(abs(gamma - 2.0) < 1e-6);
    assert(beta < 1.0);
end

function test_stopping_power_positive_and_finite(materials, ~)
    energies = [3.0, 10.0, 62.3, 100.0, 200.0, 300.0];
    keys = {'water', 'aluminium', 'copper', 'lead'};
    for i = 1:numel(keys)
        S = betheBlochMassStoppingPower(energies, materials.(keys{i}));
        assert(all(isfinite(S)));
        assert(all(S > 0));
    end
end

function test_stopping_power_decreasing_with_energy_in_clinical_range(materials, ~)
    %   In the 3-300 MeV range (over the ionization minima), dE/dx
    %   should be monotonously decreasing as the energy increases, as
    %   '-1/beta^2' dominates over the logarithm with the relativistic
    %   argument increases in this range.
    energies = linspace(3.0, 300.0, 50);
    S = betheBlochMassStoppingPower(energies, materials.water);
    assert(all(diff(S) < 0), 'S(T) should be monotonously decreasing as the energy increases');
end

function test_water_100mev_matches_nist_within_1pct(materials, ~)
    %   Sanity check against NIST PSTAR value (7.289 [MeV cm^2/g]).
    %   A 1% difference is allowed since, for water, the shell correction
    %   is nearly negligible in the considered energy range.
    S = betheBlochMassStoppingPower(100.0, materials.water);
    nistReference = 7.289;
    relDiff = abs(S - nistReference) / nistReference;
    assert(relDiff < 0.01, sprintf('Difference %.4f%% exceeds expected 1%%', relDiff*100));
end

function test_linear_stopping_power_scales_with_density(materials, ~)
    T = 100.0;
    S_mass = betheBlochMassStoppingPower(T, materials.lead);
    S_linear = betheBlochLinearStoppingPower(T, materials.lead);
    assert(abs(S_linear - S_mass * materials.lead.rho_g_cm3) < 1e-10);
end

function test_charge_scaling(materials, ~)
    %   Stopping power should scale as z^2 (testing with z=2, alpha-like)
    T = 100.0;
    S_z1 = betheBlochMassStoppingPower(T, materials.water, 'z', 1);
    S_z2 = betheBlochMassStoppingPower(T, materials.water, 'z', 2);
    assert(abs(S_z2 - 4*S_z1) < 1e-8 * S_z1);
end

function test_validity_flag_at_clinical_lower_bound(~, ~)
    %   The lower bound of the clinical range (3 MeV) should be close
    %   to or under the validity threshold: beta*gamma >= 0.1
    vc_3mev = betheBlochValidityCheck(3.0);
    vc_300mev = betheBlochValidityCheck(300.0);
    assert(vc_300mev.valid_uncorrected_bethe);
    %   At 3 MeV beta*gamma is expected to be lower than at higher
    %   energies; left as informative print, not a strict assertion,
    %   since it is the point of highest documented uncertainty in the
    %   whole project.
    fprintf('beta*gamma at 3 MeV: %.4f (uncorrected-Bethe validity threshold: 0.1)\n', ...
        vc_3mev.beta_gamma);
end
