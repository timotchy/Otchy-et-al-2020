function [modelstrct,z,u,eta,W0,sigma0,psi0] = timing_var_EM(X,Q,varargin)
% [modelstrct,z,u,eta,W0,sigma0,psi0] = timing_var_EM(X,Q,varargin)
% Expectation-maximization algorithm to estimate parameters and latent
% variables for the timing variability model described in Glaze & Troyer (2012).
%
%
% NOTE: In manuscript Sigma and Psi refer to component covariance
% matrices, while for effiency in the code, sigma and psi simply represent the diagonals
% of those matrices
%
%
% FUNCTION ARGUMENTS
% X: NxK matrix of data, N = # of samples, K is # of dimensions (i.e. unique time
% intervals per sequence)
%
% Q: Dimensionality of global variable.
%
% W0, sigma0, psi0 (optional): Optional arguments giving initial conditions for
% parameter estimates. If none are given they are generated here.
%
% conv_crit, loglkflg, maxiter (optional): conv_crit is the criterion for model convergence. 
% If loglkflg=1, conv_crit is interpreted as a percentage decrease in log-likelihood. If loglkflg=0,
% conv_crit is interpeted as a percentage change in model parameters
% (although see psiopt argument). maxiter is the maximum number of EM
% interactions allowed. Defaults are conv_crit=0.01%, loglkflg=0, maxiter=10000. 
%
% rotopt (optional): If rotopt=1, W is rotated to maximize the sum of
% values in the first column while the remaining columns are orthogonalized
% and sorted by eigenvalue (described in main manuscript). Default is 1.
%
% psiopt (optional): If psiopt=1, independent parameters are used (along with others) to
% evaluate convergence (provided that loglkflg=0). Default is 1.
%
% Optional arguments can be specified in any combination or order but
% require strings preceding each (as in other Matlab functions). For example:
%
% >> [W,sigma,psi,loglk,z,u,eta] = timing_var_EM(X,3,'sigma0',sigma0tmp,'loglkflg',1);
%
% will return estimates for data X, global dimension=3, with initial
% conditions on sigma specified by vector sigma0tmp; other initial conditions determined in
% this function; log-likelihood used to evaluate convergence; the
% default 0.01% change in log-likelihood used as the convergence
% criterion; and all other arguments set by default.
%
%
% FUNCTION RETURNS
% 
% modelstrct: Structure containing model parameter estimates and
% information. Fields represent the following:
%
% modelstrct.W: KxQ matrix of factor weights, K is # of dimensions in data set, Q is #
% of global dimensions.
%
% modelstrct.sigma: Kx1 vector of jitter variances.
%
% modelstrct.psi: Kx1 vector of independent variances.
%
% modelstrct.loglk: Log likelihood of the model.
%
% modelstrct.iter: Number of iterations needed for convergence.
%
% modelstrct.degflg: Flag indicating whether routine was terminated due to a
% degenerate matrix.
%
%
% z: NxQ matrix with global latent variables.
%
% u: NxK-1 matrix of jitter latent variables.
%
% eta: NxK matrix of independent latent variables.
%
% W0,sigma0,psi0: Initial parameter estimates.
%




% First get dimensions of X, subtract off the mean, generate covariance 
% and differencing matrices so any required initial conditions can be generated.

[N,K] = size(X);
X = X - repmat(mean(X),N,1);
Sdata = X'*X/(N);
vr = diag(Sdata);
D = -diff(eye(K))';

% Get user-secified arguments and assign values to the remaining
allargs = {'conv_crit','loglkflg','maxiter','rotopt','psiopt'};
args2timing_vars(allargs,varargin);

% Initialize parameter estimates if not given by user
initvars = {'W0','sigma0','psi0'};
for varind = 1:length(initvars)
    argind = find(strcmp(varargin,initvars{varind}));
    if isempty(argind)
        switch initvars{varind}
            case 'W0'
                W0 = repmat(sqrt(vr),1,Q).*rand(K,Q);
            case 'sigma0'
                sigma0 = .5*mean(vr)*rand(K-1,1);
            case 'psi0'
                psi0 = vr.*rand(K,1);           
        end
    else
        eval([varargin{argind} ' = varargin{argind+1};'])
    end
  
end

% assign initial conditions to parameter estimates
W = W0;
sigma = sigma0;
psi = psi0;

% positive deteriminants of component covariance matrices will be checked
% for in each iteration (matrices are diagonal)
detsigma = prod(sigma);
detpsi = prod(psi);

% deltaflg = 1 indicates lack of convergence
deltaflg = 1;

% degflg indicates degeneracy of at least 1 non-global matrix (expected global
% variance = 1 always)
degflg = detsigma > 0 & detpsi > 0;

% concatenate variables for efficiency later on (described in manuscript)
A = [W D];
Qones = ones(Q,1);
phi = [Qones;sigma];

% throughout S is the estimated covariance matrix
S = A*diag(phi)*A' + diag(psi);

if loglkflg
    % computed once for effiency
    Klog2pi = K*log(2*pi);
    
    loglk = -N/2*(Klog2pi + log(det(S)) + trace(inv(S)*Sdata));
end


% if parameter convergence is required, and user chooses to ignore psi or
% sets 0 dimensions to the global variable, then set respective default changes 
% across EM iteractions to be 0 so those parameters are effectively ignored 
% when convergence is checked for
if ~loglkflg
    if psiopt==0
        psidiff = 0;
    end
    
    if Q==0
        Wdiff = 0;
    end
end

iter = 1;

while iter < maxiter && degflg && deltaflg
    
    % EM algorithm for estimating parameters as described in Glaze & Troyer
    % (2012)
    
    % E-step
    invpsi = diag(1./psi);
    G = inv(diag(1./phi) + A'*invpsi*A);
    
    % this step involves a lot of multiplication by 0 and can probably be made much
    % more efficient
    
    v = X*invpsi*A*G;

    % M-step
    
    % These products are used to estimate parameters, computed once for efficiency
    Xv = X'*v;
    v2sm = v'*v + N*G;
    
    % If parameter convergence is required, save previous values
    if ~loglkflg
        sigma_prev = sigma;
        W_prev = W;
        psi_prev = psi;
    end
    
    sigma = diag(v2sm(Q+1:end,Q+1:end)) / N;
    W = (Xv(:,1:Q)-D*v2sm(Q+1:end,1:Q)) * inv(v2sm(1:Q,1:Q));
    
    A = [W D];
    phi = [Qones;sigma];
    
    psi = diag(Sdata + A*v2sm*A'/N - 2*Xv*A'/N);
    
    if loglkflg
        loglk_prev = loglk;
        S = A*diag(phi)*A' + diag(psi);
        loglk = -N/2*(Klog2pi + log(det(S)) + trace(inv(S)*Sdata));
    end
    
    
    % check for convergence
    
    if loglkflg
        deltaflg = ~(loglk - loglk_prev > 0 && 100*abs((loglk - loglk_prev)/loglk_prev) < conv_crit);
    else
        
        sigmadiff = 100 * max(abs(1 - sqrt(sigma) ./ sqrt(sigma_prev)));
        
        % only include psi if that option is chosen by user or default
        if psiopt
            psidiff = 100 * max(abs(1 - sqrt(psi) ./ sqrt(psi_prev)));
        end
        
        % only include W if Q>0
        if Q > 0
            Wdiff = 100 * max(abs(1 - W(:) ./ W_prev(:)));
        end
        
        deltaflg = ~(sigmadiff < conv_crit && Wdiff < conv_crit && psidiff < conv_crit);
    end
    
    % re-check for degeneracy in model
    detsigma = prod(sigma);
    detpsi = prod(psi);
    
    degflg = detsigma > 0 & detpsi > 0;
    iter = iter + 1;
    
end

% If parameters used for convergence compute log-likelihood now
if ~loglkflg
    S = A*diag(phi)*A' + diag(psi);
    loglk = -N/2*(K*log(2*pi) + log(det(S)) + trace(inv(S)*Sdata));
end

% Rotate W if option is chosen and Q>0
if rotopt==1 && Q>0
    W = rotate_max(W);
end

modelstrct.W = W;
modelstrct.sigma = sigma;
modelstrct.psi = psi;
modelstrct.loglk = loglk;
modelstrct.iter = iter;
modelstrct.degflg = degflg;

% If latent variables required for output, compute them (v not used in case
% W has been rotated, requiring similar transformation for z)
if nargout > 1
    [z,u,eta] = timing_latvar(X,W,sigma,psi);
end