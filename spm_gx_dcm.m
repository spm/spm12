function [y] = spm_gx_dcm(x,u,P,M)
% Simulated BOLD response to input
% FORMAT [y] = spm_gx_dcm(x,u,P,M)
% y          - BOLD response (%)
% x          - state vector     (see spm_fx_dcm)
% P          - Parameter vector (see spm_fx_dcm)
% M          - model specification structure (see spm_nlsi)
%__________________________________________________________________________
%
% This function implements the BOLD signal model described in: 
%
% Stephan KE, Weiskopf N, Drysdale PM, Robinson PA, Friston KJ (2007)
% Comparing hemodynamic models with DCM. NeuroImage 38: 387-401.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston & Klaas Enno Stephan
% $Id: spm_gx_dcm.m 3705 2010-02-01 20:51:28Z karl $
 
 
% Biophysical constants for 1.5T
%==========================================================================
 
% time to echo (TE) (default 0.04 sec)
%--------------------------------------------------------------------------
try, TE = M.TE; catch, TE = 0.04; end
 
% resting venous volume
%--------------------------------------------------------------------------
V0  = 100*0.04;
 
% slope r0 of intravascular relaxation rate R_iv as a function of oxygen 
% saturation Y:  R_iv = r0*[(1-Y)-(1-Y0)]
%--------------------------------------------------------------------------
r0  = 25;   % [Hz]
 
% frequency offset at the outer surface of magnetized vessels
%--------------------------------------------------------------------------
nu0 = 40.3; % [Hz]
 
% estimated region-specific resting oxygen extraction fractions
%--------------------------------------------------------------------------
E0  = P.H(:,5);
 
% estimated region-specific ratios of intra- to extra-vascular components of
% the gradient echo signal (prior mean = 1, log-normally distributed 
% scaling factor) 
%--------------------------------------------------------------------------
ep  = exp(P.H(:,6));
 
%-Coefficients in BOLD signal model
%==========================================================================
k1  = 4.3.*nu0.*E0.*TE;
k2  = ep.*r0.*E0.*TE;
k3  = 1 - ep;
 
%-Output equation of BOLD signal model
%==========================================================================
v   = exp(x(:,4));
q   = exp(x(:,5));
y   = V0*(k1.*(1 - q) + k2.*(1 - q./v) + k3.*(1 - v));
