function [P, Q, R] = siderv(Phi, H, E, G, L0)
%   [P, Q, R] = siderv(Phi, H, E, G, L0)
%
   P = lyapns(Phi, Phi-E*H, G*E');
   [U S V] = svd(P);
   P = U*S*U';
   iE = pinv(E);
   Q = iE*(P-Phi*P*Phi')*iE';
   R = L0 - H*P*H' - Q;
%
end        