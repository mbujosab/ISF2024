function [E,Q,U,iU,P] = sstoinn(Phi, E, H, C, Q, S, R)
%
% SSTOINN  - Transforms a general SS noise structure to an innovations structure.
%        [E,Q,U,iU,P] = sstoinn(Phi, E, H, C, Q, S, R)
%
%   U = cholp(Q)
%   iU = inv(U')
%   P  < Solution of ARE (algebraic Riccati equation)

   global E4OPTION

   scaleb  = E4OPTION(2);
   zeps    = E4OPTION(15);

   l = size(Phi,1);
   m = size(H,1);
   N = C*R*C';
   iN = pinv(C*R*C');
   ESCt = E*S*C';
   M = E*Q*E' - ESCt*iN*ESCt';
   A = Phi - ESCt*iN*H;
   if trace(M) < zeps && ~sum(abs(eig(A)) > 1)
      P=zeros(l);
   else
      [U,S,V] = svd([N;-H']);
      T = [U(m+1:l+m,m+1:l+m)'*A' zeros(l); -M   eye(l) ];
      Y = [U(m+1:l+m,m+1:l+m)' U(1:m,m+1:l+m)'*H; zeros(l) A ];


      %% begin MBB
      if any(exist ('OCTAVE_VERSION'))
	[AAAA, BBBB, QQQQ, ZZZZ, W, WWWW, d] = qz (T, Y);
      else
	[W,d] = eig(T,Y);
	d = diag(d);
      end			
      %% End MBB
      [e,index] = sort(abs(d));
      WW = W(:,index(1:l));
      P = real(WW(l+1:2*l,:)*pinv(WW(1:l,:)));
   end
   
   Q = H*P*H'+N;
   if scaleb, U = cholp(Q, abs(Q));
   else       U = cholp(Q); end
   iU = eye(m)/U';
   E = (Phi*P*H' + ESCt)*iU'*iU;
