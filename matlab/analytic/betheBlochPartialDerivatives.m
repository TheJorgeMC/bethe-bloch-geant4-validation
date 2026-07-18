function d = betheBlochPartialDerivatives(K, ZA, me, I, beta, gamma, z)
%BETHEBLOCHPARTIALDERIVATIVES Partial derivatives of the mass stopping
%power S(K, ZA, me, I; beta, gamma, z), used for uncertainty
%propagation.
%
%   S = K * z^2 * ZA * (1/beta^2) *
%       [ ln(2*me*beta^2*gamma^2/I) - beta^2 ]
%
%   These four expressions were derived symbolically in the Python
%   companion notebook (SymPy, see cell 12) and are HARDCODED here
%   rather than computed via a symbolic toolbox, since this project
%   deliberately uses only base MATLAB (no toolboxes required, per the
%   Fase 0 checklist). The printed SymPy output was:
%
%       dS/dK  = ZA*z^2*(-beta^2 + log(2*beta^2*gamma^2*me/I))/beta^2
%       dS/dZA = K*z^2*(-beta^2 + log(2*beta^2*gamma^2*me/I))/beta^2
%       dS/dme = K*ZA*z^2/(beta^2*me)
%       dS/dI  = -K*ZA*z^2/(I*beta^2)
%
%   Returns a struct d with fields dK, dZA, dme, dI (each the same
%   shape as beta/gamma, i.e. broadcastable over an energy grid).

    bracket = log(2.0 .* beta.^2 .* gamma.^2 .* me ./ I) - beta.^2;

    d.dK  = ZA .* z.^2 .* bracket ./ beta.^2;
    d.dZA = K .* z.^2 .* bracket ./ beta.^2;
    d.dme = K .* ZA .* z.^2 ./ (beta.^2 .* me);
    d.dI  = -K .* ZA .* z.^2 ./ (I .* beta.^2);

end
