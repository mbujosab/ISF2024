function zd = propaga(Phi, Gam, H, D, x1, u)

n = size(u,1);
r = size(u,2);
if isempty(Gam), Gam = zeros(size(Phi,1),r); end
if isempty(D), D = zeros(size(H,1),r); end
zd = zeros(n,size(H,1));

for i = 1:size(u,1)
    zd(i,:) = (H*x1 + D*u(i,:)')';
    x1 = Phi*x1 + Gam*u(i,:)';
end    