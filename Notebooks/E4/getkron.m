function ki = getkron(Phi,H)

n = size(Phi,1);
m = size(H,1);
ki = zeros(m,1);
nk = 0;
k = (1:m)';
m2 = m;

if n == 0, return; end

C = obsv(Phi,H);

for i=1:n
    ik = ones(m2,1) == 1;
    for j=1:m2
        r = rank(C(1:(i-1)*m+k(j),:));
        if r > nk
           ki(j) = ki(j) + 1;
           nk = nk + 1;
        else
           ik(j) = 0;
        end        
    end
    if nk == n, break; end
    m2 = sum(ik);
    k = k(ik);
end          

