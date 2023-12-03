function  svc = nid1(Nm, i, nvec, S1);


N = Nm(1); m = Nm(2);

Ct = exp(-2)*N^(-.87)*i^(1.63);  % NIDC

svc = zeros(size(nvec,2),2);
   for k = 1:size(nvec,2)
       dn = 2*nvec(k)*m;     % A lo Hannan-Deistler, sin exogenas
       svc(nvec(k)+1,1) = S1(nvec(k)+1)^2+log(N)*dn/N;
       svc(nvec(k)+1,2) = S1(nvec(k)+1)^2 + dn*Ct;  
   end
   
svc = svc(nvec+1,:);