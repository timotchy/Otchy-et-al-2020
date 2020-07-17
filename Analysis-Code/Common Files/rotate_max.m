function V = rotate_max(W)

[K,Q] = size(W);
One = ones(K,1);
V1 = W*W'*One / norm(W'*One);

[Vr,Lambda] = eig(W*W' - V1*V1');
[lambda,sortinds] = sort(diag(Lambda),'descend');
matinds = sortinds(1:Q-1);
Vr = Vr(:,matinds) * sqrt(Lambda(matinds,matinds));

V = [V1 Vr];
