function energies = betheBlochEnergy()
%BETHEBLOCHENERGYGRID Energy grid: 44 log-spaced points in the clinical
%range of 3-300 MeV, with an extended margin up to 1000 MeV (10
%additional points) for visualization purposes only. Identical to
%python/analytic build_energy_grid(), for direct cross-validation.

    clinical = logspace(log10(3.0), log10(300.0), 44);
    extension = logspace(log10(300.0), log10(1000.0), 10);
    energies = [clinical, extension];

end
