function [z,u,eta] = timing_latvar(X,W,sigma,psi)
% [z,u,eta] = timing_latvar(X,W,sigma,psi)
%
% Computes latent variables for Glaze & Troyer (2012) timing variability
% model. 

[N,K] = size(X);
[Ktmp,Q] = size(W);

% Make sure X is 0-mean.
X = X - repmat(mean(X),N,1);

D = -diff(eye(K))';
phi = [ones(Q,1);sigma]; 
invpsi = diag(1./psi);

A = [W D];
G = inv(diag(1./phi) + A'*invpsi*A);
v = X*invpsi*A*G;

z = v(:,1:Q);
u = v(:,Q+1:end);
eta = X - v*A';