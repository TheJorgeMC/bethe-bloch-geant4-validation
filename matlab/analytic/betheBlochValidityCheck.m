function vc = betheBlochValidityCheck(T_MeV, particleMassMeV)
%BETHEBLOCHVALIDITYCHECK Returns a struct with beta, gamma, beta*gamma,
%and a logical flag indicating whether the point is within the
%validity range of the Bethe formula WITHOUT corrections
%(beta*gamma >= 0.1, i.e. T >= 4.68 MeV for protons). See Sec.
%"Project Scope and Future Work" of the derivation document.

    if nargin < 2
        particleMassMeV = [];
    end
    [beta, gamma] = betheBlochBetaGamma(T_MeV, particleMassMeV);
    bg = beta .* gamma;

    vc.beta = beta;
    vc.gamma = gamma;
    vc.beta_gamma = bg;
    vc.valid_uncorrected_bethe = bg >= 0.1;

end
