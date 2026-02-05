function f = f_OMA_model(Re,eD)
% THIS FILE IS TO EXTRACT THE COEFFICIENTS OF THE OMA REGRESSION, IT
% IS NOT A PART OF THE MAINSCRIPT.
%F_OMA_MODEL OMA correlation: friction factor f(Re,e/D).
%
% Re, eD can be scalars or arrays; they will be broadcast to a common size.

% Order must match coeffnames: Avisc, A, B, C, D, eta, k, n, beta, Retr, r, p, Aturb
Avisc = 0.9801;      
A     = 0.5045;      
B     = 0.0;      
C     = 0.0318;    
D     = 0.7407;      
eta   = -0.2240;    
k     = 5.4044;     
n     = 0.5267;     
beta  = 0.0035;   
Retr  = 52.8455;     
r     = 0.8501;     
p     = 1.2709;     
Aturb = 0.7684;
% --------------------------------------------------

% --- broadcast Re, eD to the same shape ---
if isscalar(Re) && ~isscalar(eD)
    Re = Re .* ones(size(eD));
elseif ~isscalar(Re) && isscalar(eD)
    eD = eD .* ones(size(Re));
elseif ~isequal(size(Re),size(eD))
    error('f_OMA_model:Re_eD_size','Re and eD must be same size or scalar.');
end

Re = Re(:);
eD = eD(:);

% --- EXACT formula from fittype ---
z   = -beta .* (eD.^r) .* (Re - Retr ./ (eD.^p));
den = 1 + exp(z);

term1 = Avisc .* (Re.^eta) .* (A + B.*eD.^k).^2 .* exp(z) ./ den;
term2 = Aturb .* (A + B.*eD.^k) .* (C + D.*eD.^n) ./ den;

f_vec = term1 + term2;

% sanitize
f_vec(~isfinite(f_vec) | f_vec <= 0) = NaN;

f = reshape(f_vec,size(Re));

end
