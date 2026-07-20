function energies = betheBlochEnergyGrid()
%BETHEBLOCHENERGYGRID Loaded directly from the NIST PSTAR standard
%grid (nist_data/energy_grid.csv) -- see
%python/analytic/generate_curves.py for the Python twin of this
%function. Both MUST read the same file.

    thisDir = fileparts(mfilename('fullpath'));
    gridPath = fullfile(thisDir, '..', '..', 'nist_data', 'energy_grid.csv');
    T = readtable(gridPath);
    energies = T.energy_MeV';

end