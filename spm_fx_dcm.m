function [y] = spm_fx_dcm(x,u,P,M)
% state equation for a dynamic [bilinear/nonlinear/Balloon] model of fMRI
% responses
% FORMAT [y] = spm_fx_dcm(x,u,P,M)
% x      - state vector
%   x(:,1) - neuronal acivity                         u
%   x(:,2) - vascular signal                          s
%   x(:,3) - rCBF                                  ln(f)
%   x(:,4) - venous volume                         ln(v)
%   x(:,5) - deoyxHb                               ln(q)
% y      - dx/dt
%___________________________________________________________________________
%
% References for hemodynamic & neuronal state equations:
% 1. Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
%    changes during brain activation: The Balloon model. MRM 39:855-864,
%    1998.
% 2. Friston KJ, Mechelli A, Turner R, Price CJ. Nonlinear responses in
%    fMRI: the Balloon model, Volterra kernels, and other hemodynamics.
%    Neuroimage 12:466-477, 2000.
% 3. Stephan KE, Kasper L, Harrison LM, Daunizeau J, den Ouden HE,
%    Breakspear M, Friston KJ. Nonlinear dynamic causal models for fMRI.
%    Neuroimage 42:649-662, 2008.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston & Klaas Enno Stephan
% $Id: spm_fx_dcm.m 3705 2010-02-01 20:51:28Z karl $


% hemodynamic parameters
%--------------------------------------------------------------------------
%   H(1) - signal decay                                   d(ds/dt)/ds)
%   H(2) - autoregulation                                 d(ds/dt)/df)
%   H(3) - transit time                                   (t0)
%   H(4) - exponent for Fout(v)                           (alpha)
%   H(5) - resting oxygen extraction                      (E0)
%   H(6) - ratio of intra- to extra-vascular components   (epsilon)
%          of the gradient echo signal

% Neuronal motion
%==========================================================================

% effective intrinsic connectivity with bilinear (input) terms
%--------------------------------------------------------------------------
for i = 1:size(P.B,3)
    P.A = P.A + u(i)*P.B(:,:,i);
end

% and nonlinear (state) terms
%--------------------------------------------------------------------------
for i = 1:size(P.D,3)
    P.A = P.A + x(i,1)*P.D(:,:,i);
end

% implement differential state equation y = dx/dt (neuronal)
%--------------------------------------------------------------------------
y         = x;
y(:,1)    = exp(P.b)*(P.A*x(:,1) + P.C*u);


% Hemodynamic motion
%==========================================================================

% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------
x(:,3:5)  = exp(x(:,3:5));

% Fout = f(v) - outflow
%--------------------------------------------------------------------------
fv        = x(:,4).^(1./P.H(:,4));

% e = f(f) - oxygen extraction
%--------------------------------------------------------------------------
ff        = (1 - (1 - P.H(:,5)).^(1./x(:,3)))./P.H(:,5);

% implement differential state equation y = dx/dt (hemodynamic)
%--------------------------------------------------------------------------
y(:,2)    = x(:,1)  - P.H(:,1).*x(:,2) - P.H(:,2).*(x(:,3) - 1);
y(:,3)    = x(:,2)./x(:,3);
y(:,4)    = (x(:,3) - fv)./(P.H(:,3).*x(:,4));
y(:,5)    = (ff.*x(:,3) - fv.*x(:,5)./x(:,4))./(P.H(:,3).*x(:,5));
y         = y(:);
