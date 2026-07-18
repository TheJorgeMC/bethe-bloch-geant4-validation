function w = betheBlochWaterComposition()
%BETHEBLOCHWATERCOMPOSITION Liquid water Bragg-rule constituents.
%
%   Composition by mass fraction (NIST STAR, matno = 276):
%   H w=0.111894, O w=0.888106.
%   Atomic weights (IUPAC Table 2, "List of Elements in Atomic Number
%   Order", https://iupac.qmul.ac.uk/AtWt/): H = 1.0080(2), O = 15.999(1).
%
%   Used both to compute WATER.Z_over_A (via betheBlochMaterials) and to
%   propagate Delta(Z/A) for water (via betheBlochDeltaZA), so the same
%   source values feed both the central value and its uncertainty.

    w.A_H = 1.0080;       w.Delta_A_H = 0.0002;   % IUPAC: H = 1.0080(2)
    w.A_O = 15.999;       w.Delta_A_O = 0.001;    % IUPAC: O = 15.999(1)
    w.w_H = 0.111894;     w.w_O = 0.888106;       % NIST STAR mass fractions
    w.Z_H = 1;            w.Z_O = 8;

end
