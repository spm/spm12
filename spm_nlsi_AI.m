function [Ep] = spm_nlsi_AI(M,Y,U)
% Bayesian inversion of a linear-nonlinear model of the form F(p)*G(g)'
% FORMAT [Ep,Eg,Cp,Cg,S,F,L]= spm_nlsi_N(M,U,Y)
%
% Generative model
%__________________________________________________________________________
% 
% M.IS - IS(p,M,U) A prediction generating function name; usually an 
%        integration scheme for state-space models of the form
%
%        M.f  - f(x,u,p,M) - state equation:  dxdt = f(x,u)
%
%        that returns hidden states - x; however, it can be any nonlinear
%        function of the inputs u. I.e., x = IS(p,M,U)
%
% M.G  - G(g,M) - linear observer: y = (x - M.x')*G(g,M)'
%
% M.FS - function name f(y,M) - feature selection
%        This [optional] function performs feature selection assuming the
%        generalized model y = FS(y,M) = FS(x*G',M) + X0*P0 + e
%
% M.x  - The expansion point for the states (i.e., the fixed point)
%
% M.P  - starting estimates for model parameters [ states - optional]
% M.Q  - starting estimates for model parameters [ observer - optional]
%
% M.pE - prior expectation  - of model parameters - f(x,u,p,M)
% M.pC - prior covariance   - of model parameters - f(x,u,p,M)
%
% M.gE - prior expectation  - of model parameters - G(g,M)
% M.gC - prior covariance   - of model parameters - G(g,M)
%
% M.hE - prior expectation  - E{h}   of log-precision parameters
% M.hC - prior covariance   - Cov{h} of log-precision parameters
%
% U.u  - inputs
% U.dt - sampling interval
%
% Y.y  - {[ns,nx],...} - [ns] samples x [nx] channels x {trials}
% Y.X0 - Confounds or null space
% Y.dt - sampling interval for outputs
% Y.Q  - error precision components
%
%
% Parameter estimates
%--------------------------------------------------------------------------
% Ep  - (p x 1)         conditional expectation  E{p|y}
% Cp  - (p x p)         conditional covariance   Cov{p|y}
%
% Eg  - (p x 1)         conditional expectation  E{g|y}
% Cg  - (p x p)         conditional covariance   Cov{g|y}
%
% S   - (v x v)         [Re]ML estimate of error Cov{e(h)}
%
% log evidence
%--------------------------------------------------------------------------
% F   - [-ve] free energy F = log evidence = p(y|m)
% 
%     L(1) = - ey'*iS*ey/2;             accuracy of states
%     L(2) = - ep'*ipC*ep/2;            accuracy of parameters (f)
%     L(3) = - eg'*igC*eg/2;            accuracy of parameters (g)
%     L(4) = - eu'*iuC*eu/2;            accuracy of parameters (u)
%     L(5) = - eh'*ihC*eh/2;            accuracy of precisions (u)
%     L(6) = - ns*nr*log(8*atan(1))/2;  constant
%     L(7) = - nq*spm_logdet(S)/2;      precision
%     L(8) = spm_logdet(ibC*Cb)/2;      parameter complexity
%     L(9) = spm_logdet(ihC*Ch)/2;      precision complexity
%
%__________________________________________________________________________
% Returns the moments of the posterior p.d.f. of the parameters of a
% nonlinear model specified by IS(P,M,U) under Gaussian assumptions. Usually,
% IS would be an integrator of a dynamic MIMO input-state-output model 
%
%              dx/dt = f(x,u,p)
%              y     = G(g)*x  + X0*B + e
%
% The E-Step uses a Fisher-Scoring scheme and a Laplace
% approximation to estimate the conditional expectation and covariance of P
% If the free-energy starts to increase, a Levenberg-Marquardt scheme is
% invoked.  The M-Step estimates the precision components of e, in terms
% of [Re]ML point estimators of the log-precisions.
% An optional feature selection can be specified with parameters M.FS
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_nlsi_AI.m 7679 2019-10-24 15:54:07Z spm $
 

% setup and initialise
%==========================================================================
GRAPH = 1;
n     = 4;


% sufficient statistics of model (P) and functional parameters (B)
%--------------------------------------------------------------------------
Ep.p  = M.pE;
Ep.h  = M.hE;
Cp.p  = M.pC;
Cp.h  = M.hC;

% priors and posteriors of vectorised model parameters
%--------------------------------------------------------------------------
EP    = spm_vec(Ep);
CP    = spm_cat(spm_diag({Cp.p,Cp.h}));

% initial functional parameters
%--------------------------------------------------------------------------
np    = spm_length(EP);
d     = {1};
for i = 1:n   
   bE{i} = zeros(d{:});
   bC{i} = eye(spm_length(bE{i}));
   d     = [{np} d];
end

% functional priors
%--------------------------------------------------------------------------
PP    = spm_inv(CP);
bE{1} = zeros(1,1);
bE{2} = zeros(np,1);
bE{3} = speye(np,np);
bC{1} = exp(16);
bC{2} = speye(np,np)*2;
bC{3} = diag(spm_vec(speye(np,np)) + 1)*2;

% priors and posteriors of vectorised functional parameters
%--------------------------------------------------------------------------
BE    = spm_vec(bE);
BC    = spm_cat(spm_diag(bC));
EB    = BE;
CB    = BC;
eB    = bE;

% amplitude of higher order terms and expected free energy functional
%--------------------------------------------------------------------------
spm_Q = @(p) p'*p + exp(-4);
spm_E = @(G,CB,Q) log(Q)/2 - log(G'*CB*G + Q)/2;
spm_F = @(p,b,W) spm_G(p,numel(b))'*spm_vec(b); % + spm_logdet(W'*spm_H(p,b)*W)/2;


% sample
%--------------------------------------------------------------------------
nc    = 1024;                                 % number of samples
h     = 1;
TL    = - M.L(spm_unvec(EP,Ep),Y,M);
for i = 1:128
    
    % get candidate sample from model posterior
    %======================================================================
    
    % sample points
    %----------------------------------------------------------------------
    W     = spm_sqrtm(PP);
    U     = spm_inv(W);
    for j = 1:nc
        q(:,j) = rand(np,1) - 1/2;
    end
    q   = kron([1/32,1/3,1/2],spm_en(spm_sigma(np)));
    nc  = size(q,2);
    
    
    % expected free energy E (with distance from current sample Q)
    %----------------------------------------------------------------------
    for j = 1:nc
%         F(j) = spm_F(q(:,j),eB,W);
        E(j) = spm_E(spm_G(q(:,j),n),CB,h*spm_Q(q(:,j)));
    end
    
    % select and evaluate sample
    %----------------------------------------------------------------------
    [E,j]  = min(E);
    p(:,i) = EP + U*q(:,j);    
    L(i,1) = - M.L(spm_unvec(p(:,i),Ep),Y,M);
    
    % invert functional model
    %======================================================================
    
    % variance component (higher order terms)
    %----------------------------------------------------------------------
    for j = 1:i
        qj     = W*(p(:,j) - EP);
        Q(j,j) = spm_Q(qj);
        X(:,j) = spm_G(qj,n);
    end   
    
    % parametric empirical Bayes
    %----------------------------------------------------------------------
    P{1}.X = X';
    P{2}.X = BE;
    P{1}.C = {Q};
    P{2}.C = BC;
    
    C     = spm_PEB(L,P);
    EB    = C{2}.E;
    CB    = C{2}.C;
    h     = C{1}.h
    eB    = spm_unvec(EB,bE);
    
    % model parameters
    %----------------------------------------------------------------------
    H     = spm_H(EP,eB);
    PP    = W'*spm_sqrtm(H'*H)*W;
    CP    = spm_inv(PP);
    TP    = EP - U*eB{2};
%     [E,j] = min(F);
%     EP    = EP + U*(q(:,j));
    tL    = - M.L(spm_unvec(TP,Ep),Y,M);
    if tL < TL
        EP = TP;
        TL = tL;
    end

    % ensure predicted location is better than the current sampling
    %----------------------------------------------------------------------
    if L(i,1) < TL
        EP = p(:,i);
        TL = L(i,1);
    end
    
    
    % evaluate (positive) free energy and convergence
    %----------------------------------------------------------------------
    FE(i) = - spm_logdet(CP)/2 - M.L(spm_unvec(EP,Ep),Y,M);
    
    % graphics if requested
    %======================================================================
    if GRAPH
        
        % get range of selected parameter space
        %------------------------------------------------------------------
        s     = 3;
        r     = 64;
        j     = 1;
        x     = linspace((EP(j) - s*sqrt(CP(j,j))),(EP(j) + s*sqrt(CP(j,j))),r);
        jE(i) = EP(j);
        jC(i) = CP(j,j);
        
        % true energy
        %------------------------------------------------------------------
        for k = 1:r
            kp    = EP;
            kp(j) = x(k);
            pk    = W*(kp - EP);
            G     = spm_G(pk,n);
            Q     = spm_Q(pk)*h;
            kL(k) = -M.L(spm_unvec(kp,Ep),Y,M);
            kF(k) = spm_E(G,CB,Q);
            kE(k) = G'*EB;
            kC(k) = G'*CB*G + Q;
        end
        kE        = kE - min(kL);
        kL        = kL - min(kL);
        
        % plot
        %------------------------------------------------------------------
        spm_figure('GetWin','AI (convergence)'); clf
        subplot(2,2,1)
        spm_plot_ci(kE,kC,x), hold on
        plot(x,kL,'r',p(j,:),zeros(1,i),'.g','MarkerSize',8), hold off
        set(gca,'YLim',[-8 64]),axis square
        
        subplot(2,2,2)
        plot(x,spm_softmax(-kL(:)),'r'), hold on
        plot(x,spm_softmax(-kE(:)),'b'), hold on
        plot(x,spm_softmax(-kF(:)),'g'), hold off
        axis square
        
        subplot(2,2,3)
        spm_plot_ci(jE(:),jC(:),1:i)
        axis square
        
        subplot(2,2,4)
        plot(L,X'*EB,'o',L,L,':')
        axis square
        
        disp(full(eB{3}))
        disp(full(eB{2}))

        
    end
        
end

if GRAPH
    
    % get range of selected parameter space
    %------------------------------------------------------------------
    i     = 1;
    j     = np;
    xi    = linspace((EP(i) - s*sqrt(CP(i,i))),(EP(i) + s*sqrt(CP(i,i))),r);
    xj    = linspace((EP(j) - s*sqrt(CP(j,j))),(EP(j) + s*sqrt(CP(j,j))),r);

    
    % true energy
    %------------------------------------------------------------------
    for u = 1:r
        for v = 1:r
            kp      = EP;
            kp(i)   = xi(u);
            kp(j)   = xj(v);
            kL(u,v) = -M.L(spm_unvec(kp,Ep),Y,M);
            kE(u,v) = spm_G(W*(kp - EP),n)'*EB;
        end
    end
    sL    = reshape(spm_softmax(-kL(:)),r,r);
    sE    = reshape(spm_softmax(-kE(:)),r,r);
    
    % plot
    %------------------------------------------------------------------
    spm_figure('GetWin','AI (bivariate)'); clf
    subplot(2,2,1), imagesc(xj,xi,kL), title('Energy'),    axis square
    subplot(2,2,2), imagesc(xj,xi,kE), title('Estimate'),  axis square
    subplot(2,2,3), imagesc(xj,xi,sL), title('Posterior'), axis square
    subplot(2,2,4), imagesc(xj,xi,sE), title('Estimate'),  axis square
    
    
    
end

return


% auxiliary functions
%==========================================================================
function G  = spm_G(p,n)
% Vectorisation of tensor products
% FORMAT G  = spm_G(p,n)
% 
% p   - parameter vector
% n   - expansion order
%__________________________________________________________________________
% Copyright (C) 2013-2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_nlsi_AI.m 7679 2019-10-24 15:54:07Z spm $
%--------------------------------------------------------------------------
G{1}  = 1;
for i = 2:n
    G{i} = spm_cross(G{i - 1},p)/factorial(i - 1);
end
G     = spm_vec(G);

return

% auxiliary functions
%==========================================================================
function H  = spm_H(p,B)
% Hessian of expansion
% FORMAT H  = spm_H(p,b)
% 
% p   - parameter vector
% b   - functional parameters (cell format)
%__________________________________________________________________________
% Copyright (C) 2013-2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_nlsi_AI.m 7679 2019-10-24 15:54:07Z spm $
%--------------------------------------------------------------------------
H     = B{3};
x     = {p};
for i = 4:numel(B)
    H = H + spm_dot(B{i},x)/factorial(i - 3);
    x = [x{:},{p}];
end

return




% auxiliary functions
%==========================================================================
function K  = spm_sigma(n)
% Sigma points for sampling
% FORMAT K  = spm_sigma(n)
% 
% n   - number of parameters
%__________________________________________________________________________
% Copyright (C) 2013-2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_nlsi_AI.m 7679 2019-10-24 15:54:07Z spm $
%--------------------------------------------------------------------------
k     = [0 1 -1];
K     = k;
for i = 2:n
    K = [kron(K,ones(1,n));kron(ones(1,n^(i - 1)),k)];
end
return
















% get log likelihood function assuming NLSI form for generative model
%==========================================================================

function spm_nlsi_GN_L(P,U,M)
% options
%--------------------------------------------------------------------------
try, M.nograph; catch, M.nograph = 0;  end
try, M.Nmax;    catch, M.Nmax    = 64; end
try, M.Gmax;    catch, M.Gmax    = 8;  end
try, M.Hmax;    catch, M.Hmax    = 4;  end

% figure (unless disabled)
%--------------------------------------------------------------------------
if ~M.nograph
    Fsi = spm_figure('GetWin','AI');
end
 
% check integrator
%--------------------------------------------------------------------------
try
    IS = M.IS;
catch
    IS = 'spm_int_U';
end
 
% check observer has not been accidentally specified
%--------------------------------------------------------------------------
try
    M = rmfield(M,'g');
end
 
% composition of feature selection and prediction (usually an integrator)
%--------------------------------------------------------------------------
if isfield(M,'FS')
 
    % FS(y,M)
    %----------------------------------------------------------------------
    try
        y  = feval(M.FS,Y.y,M);
        try
            FS = inline([M.FS '(y,M)'],'y','M');
        catch
            FS = M.FS;
        end
 
    % FS(y)
    %----------------------------------------------------------------------
    catch
        y  = feval(M.FS,Y.y);
        FS = inline([M.FS '(y)'],'y','M');
 
    end
else
 
    % y
    %----------------------------------------------------------------------
    y  = Y.y;
    FS = inline('y','y','M');
end
 
% size of data (usually samples x channels)
%--------------------------------------------------------------------------
if iscell(y)
    
    % concatenate samples over cell, ensuring the same for predictions
    %----------------------------------------------------------------------
    ns = size(y{1},1);
    y  = spm_cat(y(:));             
    IS = inline(['spm_cat(' IS '(P,M,U))'],'P','M','U');
    
else
    ns = size(y,1);
end
ny   = length(spm_vec(y));
nr   = ny/ns;                            % number of samples and responses
M.ns = ns;                               % store in M.ns for integrator
 
% initial states
%--------------------------------------------------------------------------
try
    M.x;
catch
    try
        M.n;
    catch
        M.n = 0;
    end
    M.x = sparse(M.n,1);
end
 
% input
%--------------------------------------------------------------------------
try
    U;
catch
    U = [];
end
 
% initial parameters
%--------------------------------------------------------------------------
try
    spm_vec(M.P) - spm_vec(M.pE);
    if ~M.nograph
        fprintf('\n(state) parameter initialisation successful\n')
    end
catch
    M.P = M.pE;
end
try
    spm_vec(M.Q) - spm_vec(M.gE);
    if ~M.nograph
        fprintf('\n(observer) parameter initialisation successful\n')
    end
catch
    M.Q = M.gE;
end
 
 
% time-step
%--------------------------------------------------------------------------
try
    Y.dt;
catch
    Y.dt = 1;
end
 
 
% precision components Q
%--------------------------------------------------------------------------
try
    Q = Y.Q;
catch
    Q = spm_Ce(ns*ones(1,nr));
end
nh    = length(Q);          % number of precision components
nt    = length(Q{1});       % number of time bins
nq    = nr*ns/nt;           % for compact Kronecker form of M-step

 
% confounds (if specified)
%--------------------------------------------------------------------------
try
    if isempty(Y.X0)
        Y.X0 = sparse(ns,0);
    end
    dgdu = kron(speye(nr,nr),Y.X0);
catch
    dgdu = sparse(ns*nr,0);
end
 
% hyperpriors - expectation (and initialize hyperparameters)
%--------------------------------------------------------------------------
try
    hE  = M.hE;
    if length(hE) ~= nh
        hE = hE + sparse(nh,1);
    end
catch
    hE  = sparse(nh,1) - log(var(spm_vec(y))) + 4;
end
h       = hE;

% hyperpriors - covariance
%--------------------------------------------------------------------------
try
    ihC = spm_inv(M.hC);
    if length(ihC) ~= nh
        ihC = ihC*speye(nh,nh);
    end
catch
    ihC = speye(nh,nh)*exp(4);
end


% unpack prior covariances
%--------------------------------------------------------------------------
if isstruct(M.pC); M.pC = spm_diag(spm_vec(M.pC)); end
if isstruct(M.gC); M.gC = spm_diag(spm_vec(M.gC)); end
if isvector(M.pC); M.pC = spm_diag(M.pC); end
if isvector(M.gC); M.gC = spm_diag(M.gC); end

% dimension reduction of parameter space
%--------------------------------------------------------------------------
Vp    = spm_svd(M.pC,0);
Vg    = spm_svd(M.gC,0);
np    = size(Vp,2);                   % number of parameters (f)
ng    = size(Vg,2);                   % number of parameters (g)
nu    = size(dgdu,2);                 % number of parameters (u)

 
% prior moments
%--------------------------------------------------------------------------
pE    = M.pE;
gE    = M.gE;
uE    = sparse(nu,1);
 
% second-order moments (in reduced space)
%--------------------------------------------------------------------------
sw    = warning('off','all');
pC    = Vp'*M.pC*Vp;
gC    = Vg'*M.gC*Vg;
uC    = speye(nu,nu)*exp(16);
ipC   = spm_inv(pC);                           % p - state parameters
igC   = spm_inv(gC);                           % g - observer parameters
iuC   = spm_inv(uC);                           % u - fixed parameters
ibC   = spm_cat(spm_diag({ipC,igC,iuC}));      % all parameters
bC    = speye(size(ibC))*exp(-16);

 
% initialize conditional density
%--------------------------------------------------------------------------
Ep    = M.P;
Eg    = M.Q;
Eu    = spm_pinv(dgdu)*spm_vec(y);

% expansion point
%--------------------------------------------------------------------------
if ~isempty(M.x)
    x0 = ones(size(y,1),1)*spm_vec(M.x)';
else
    x0 = 0;
end


% EM
%==========================================================================
warning(sw); sw = warning('off','all');
criterion       = [0 0 0 0];

C.F   = -Inf;                                   % free energy
v     = -4;                                     % log ascent rate
dgdp  = zeros(ny,np);
dgdg  = zeros(ny,ng);
dFdh  = zeros(nh,1);
dFdhh = zeros(nh,nh);

 
% Optimize p: parameters of f(x,u,p)
%==========================================================================
EP     = [];
for ip = 1:M.Nmax
 
    % time
    %----------------------------------------------------------------------  
    Ti = tic;
    
    % predicted hidden states (x) and dxdp
    %----------------------------------------------------------------------
    [dxdp,x] = spm_diff(IS,Ep,M,U,1,{Vp});  
    
    % check for initial iterations and dissipative dynamics
    %----------------------------------------------------------------------
    if all(isfinite(spm_vec(x)))
        Gmax = M.Gmax;
        if ip < 8
            vg = -4;
        else
            vg = 2;
        end
    else
        Gmax = 0;
    end
    
       
    % Optimize g: parameters of G(g)
    %======================================================================
    for ig = 1:Gmax
        
        % prediction yp = G(g)*x
        %------------------------------------------------------------------
        [dGdg,G] = spm_diff(M.G,Eg,M,1,{Vg});
        yp       = FS((x - x0)*G',M);
        
        % prediction errors - states
        %==================================================================
        ey    = spm_vec(y)  - spm_vec(yp) - dgdu*Eu;
 
        % prediction errors - parameters
        %------------------------------------------------------------------
        ep    = Vp'*(spm_vec(Ep) - spm_vec(pE));
        eg    = Vg'*(spm_vec(Eg) - spm_vec(gE));
        eu    =      spm_vec(Eu) - spm_vec(uE);
 
        % gradients
        %------------------------------------------------------------------
        for i = 1:np
            dgdp(:,i) = spm_vec(FS(dxdp{i}*G',M));
        end
        try
            for i = 1:ng
                dgdg(:,i) = spm_vec(FS((x - x0)*dGdg{i}',M));
            end
        catch
            dgdg = FS((x - x0)*dGdg,M);
        end
  
        % Optimize F(h): parameters of iS(h)
        %==================================================================        
        dgdb   = [dgdp dgdg dgdu];           
        for ih = 1:M.Hmax
 
            % precision
            %--------------------------------------------------------------
            iS    = speye(nt,nt)*exp(-32);
            for i = 1:nh
                iS = iS + Q{i}*exp(h(i));
            end
            S     = spm_inv(iS);
            iS    = kron(speye(nq),iS);
            dFdbb = dgdb'*iS*dgdb + ibC;
            Cb    = spm_inv(dFdbb) + bC;
            
            % precision operators for M-Step
            %--------------------------------------------------------------
            for i = 1:nh
                P{i}  = Q{i}*exp(h(i));
                PS{i} = P{i}*S;
                P{i}  = kron(speye(nq),P{i});
            end
 
            % derivatives: dLdh = dL/dh,...
            %--------------------------------------------------------------
            for i = 1:nh
                dFdh(i,1)      =   trace(PS{i})*nq/2 ...
                                 - real(ey'*P{i}*ey)/2 ...
                                 - spm_trace(Cb,dgdb'*P{i}*dgdb)/2;
                for j = i:nh
                    dFdhh(i,j) = - spm_trace(PS{i},PS{j})*nq/2;
                    dFdhh(j,i) =   dFdhh(i,j);
                end
            end
 
            % add hyperpriors
            %--------------------------------------------------------------
            eh    = h     - hE;
            dFdh  = dFdh  - ihC*eh;
            dFdhh = dFdhh - ihC;
            Ch    = spm_inv(-dFdhh);
            
            % M-Step: update ReML estimate of h
            %--------------------------------------------------------------
            dh    = spm_dx(dFdhh,dFdh,{4});
            h     = h + min(max(dh,-2),2);
 
            % convergence
            %--------------------------------------------------------------
            if dFdh'*dh < exp(-2), break, end
 
        end
 
        % E-step: optimise F(g,u)
        %==================================================================
        
        % update gradients and curvature - counfounds
        %------------------------------------------------------------------
        dFdu  =  dgdu'*iS*ey   - iuC*eu;
        dFduu = -dgdu'*iS*dgdu - iuC;

        % Conditional updates of confounds (u)
        %------------------------------------------------------------------
        du    = spm_dx(dFduu,dFdu,{4});
        Eu    = Eu + du;
        
        % update gradients and curvature - parameters
        %------------------------------------------------------------------
        dFdg  =  dgdg'*iS*ey   - igC*eg;
        dFdgg = -dgdg'*iS*dgdg - igC;
 
        % Conditional updates of parameters (g)
        %------------------------------------------------------------------
        dg    = spm_dx(dFdgg,dFdg,{vg});
        Eg    = spm_unvec(spm_vec(Eg) + Vg*dg,Eg);
         
        % convergence
        %------------------------------------------------------------------
        dG    = dFdg'*dg;
        if ig > 1 && dG < exp(-2), break, end
        
    end
    
    % optimise objective function: F(p) = log-evidence - divergence
    %======================================================================
    L(1) = - ey'*iS*ey/2;            % accuracy
    L(2) = - ep'*ipC*ep/2;           % complexity
    L(3) = - eg'*igC*eg/2;           % complexity
    L(4) = - eu'*iuC*eu/2;           % complexity
    L(5) = - eh'*ihC*eh/2;           % complexity
    L(6) = - ns*nr*log(8*atan(1))/2; % accuracy
    L(7) = - nq*spm_logdet(S)/2;     % accuracy
    L(8) = spm_logdet(ibC*Cb)/2;     % complexity
    L(9) = spm_logdet(ihC*Ch)/2;     % complexity
    F    = sum(L);
    
    % record increases and reference log-evidence for reporting
    %----------------------------------------------------------------------
    try
        F0;
        fprintf(' actual: %.3e (%.2f sec)\n',full(F - C.F),toc(Ti))
    catch
        F0 = F;
    end
     
    % if F has increased, update gradients and curvatures for E-Step
    %----------------------------------------------------------------------
    if F > C.F || ip < 4
        
        % update gradients and curvature
        %------------------------------------------------------------------
        dFdp  =  dgdp'*iS*ey   - ipC*ep;
        dFdpp = -dgdp'*iS*dgdp - ipC;
 
        % decrease regularization
        %------------------------------------------------------------------
        v     = min(v + 1/2,4);
        str   = 'EM(+)';
 
        % accept current estimates
        %------------------------------------------------------------------
        C.Cb  = Cb;                               % conditional covariance
        C.Ep  = Ep;                               % and expectations
        C.Eg  = Eg;
        C.Eu  = Eu;
        C.h   = h;
        C.F   = F;
        C.L   = L;   
        
    else
 
        % reset expansion point
        %------------------------------------------------------------------
        Cb    = C.Cb;                             % conditional covariance
        Ep    = C.Ep;                             % and expectations
        Eg    = C.Eg;
        Eu    = C.Eu;
        h     = C.h;
 
        % and increase regularization
        %------------------------------------------------------------------
        v     = min(v - 2,-4);
        str   = 'EM(-)';
 
    end
 
    % Optimize p: parameters of f(x,u,p)
    %======================================================================
    dp    = spm_dx(dFdpp,dFdp,{v});
    Ep    = spm_unvec(spm_vec(Ep) + Vp*dp,Ep);
    
    % diagnostic
    %----------------------------------------------------------------------
    % EP(:,end + 1) = spm_vec(Ep);
 
    
    % subplot times
    %----------------------------------------------------------------------
    try
        if length(Y.pst) == size(yp,1)
            yt = Y.pst;
        else
            yt = (1:size(yp,1))*Y.dt*1000;
        end
    catch
        yt = (1:size(yp,1))*Y.dt*1000;
    end
    
    
    % graphics
    %----------------------------------------------------------------------
    if exist('Fsi', 'var')
        spm_figure('Select', Fsi)
        
        % subplot prediction
        %------------------------------------------------------------------
        subplot(3,1,1)
        try
            plot(yt,x)
            xlabel('time (ms)')
            set(gca,'XLim',[yt(1) yt(end)])
        catch
            plot(x), spm_axis tight
        end
        title(sprintf('%s: %i','E-Step: hidden states',ip))
        grid on
        
        subplot(3,1,2)
        plot(yt,yp),                        hold on
        plot(yt,yp + spm_unvec(ey,yp),':'), hold off
        xlabel('time (ms)')
        set(gca,'XLim',[yt(1) yt(end)])
        title('E-Step: response and prediction')
        grid on
 
        % subplot parameters - f(P)
        %------------------------------------------------------------------
        subplot(3,2,5)
        bar(full(Vp*ep))
        xlabel('parameter f(x)')
        title('conditional [minus prior] expectation')
        grid on
        
        % subplot parameters - g(G)
        %------------------------------------------------------------------
        subplot(3,2,6)
        bar(full(Vg*eg))
        xlabel('parameter (g(x))')
        title('conditional [minus prior] expectation')
        grid on
        drawnow
        
    end
 
    % convergence
    %----------------------------------------------------------------------
    dF  = dFdp'*dp;
    ig  = max([0 ig]);
    fprintf('%-6s: %-2i (%i,%i) %4s %-6.3e %6s %6.3e ',str,ip,ig,ih,'F:',full(C.F - F0),'dF predicted:',full(dF))
    criterion = [(dF < 1e-1) criterion(1:end - 1)];
    if all(criterion), fprintf(' convergence\n'), break, end
    
end
if exist('Fsi', 'var')
    spm_figure('Focus', Fsi)
end

% outputs
%--------------------------------------------------------------------------
Ep     = C.Ep;
Eg     = C.Eg;
Cp     = Vp*C.Cb((1:np),     (1:np)     )*Vp';
Cg     = Vg*C.Cb((1:ng) + np,(1:ng) + np)*Vg';
F      = C.F;
L      = C.L;
warning(sw);

% diagnostic
%--------------------------------------------------------------------------
% save('spm_nlsi_N_Ep','EP')

return
