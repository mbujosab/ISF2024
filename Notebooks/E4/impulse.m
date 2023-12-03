function [irf] = impulse(Phi, Gam, E, H, D, C, Q, S, R, hor)
% PLOTIRF - Displays standarized graphs of an impulse response function
%	irf = impulse(Phi, Gam, E, H, D, C, Q, S, R, hor)
% hor	> horizon of the response function.
% tit	> (mx?) matrix which contains the names of the endogenous variables;
% irf	< matrix of standarized impulse response function. Each column is the
%	  response of each endogenous variable vs. each noise. If there are
%	  any exogenous variable, the responses of the endogenous will be
%	  the noises.
% 
% 20/8/98
% Copyright (c) Jaime Terceiro, 1998



% IMPULSE RESPONSE FUNCTION FOR NOISES

% Covariance matrix

V = [Q S; S' R];

if Q == S & Q == R
  V = Q;
end


% Cholesky Decomposition

K = chol(V);
P = K';
clear K;


% Eigenvalues and eigenvectors

[W,L]=eig(V);


% Z matrix

Z = inv(W)*P;


% Impulse Response Function

[m_states,l_exog] = size(Gam);
[m_states,k_noises] = size(E);
[n_endog,l_exog] = size(D);
[n_endog,m_states]  = size(H);

Resp(1:n_endog,1:k_noises) = H*E*P;
Response = Resp;

for time_counter = 2 : hor,
  Resp = H*Phi^(time_counter-1)*E*P;
  Response = [ Response
	       Resp    ];
end

clear Resp;

Result = Response;

if n_endog >=2

  Result = [];
  Result(1,1:k_noises) = Response(1,:);

  for n = 2 : n_endog,
    Resp(1,:) = Response(n,:);
    Result = [Result Resp];
  end

  clear Resp;

  for t = 2 : hor,
    Respaux = [];
      for n = 1 : n_endog,
        Resp = Response(t*n-(t-1)*(n-n_endog),:);
        Respaux = [Respaux Resp];
      end
    Result = [ Result
	       Respaux ];
  end
end

clear Resp;
clear Respaux;
Impulse = [];



% COMPOSITE MODELS

if Gam ~= [];

  % Impulse Response Function

  clear Resp;
  clear Response;
  clear Respaux;
  for l = 1 : l_exog,
    e = zeros(l_exog,1);
    e(l,1) = 1;
    Resp = [];
    Resp(1,:) = D*e;
    Response = Resp;

    for time_counter = 2 : hor,
      Resp = H*Phi^(time_counter-1)*Gam*e;
      Response = [ Response
		   Resp    ];
    end
    Imp_mat = [Imp_mat Response];
  end
  
  Impulse = Imp_mat;
  clear Resp;

  if n_endog >=2
    Impulse = [];
    Impulse(1,:) = Imp_mat(1,:);

    for n = 2 : n_endog,
      Resp(1,:) = Imp_mat(n,:);
      Impulse = [Impulse Resp];
    end
    clear Resp;

    for t = 2 : hor,
      Respaux = [];

        for n = 1 : n_endog,
          Resp = Imp_mat(t*n-(t-1)*(n-n_endog),:);
          Respaux = [Respaux Resp];
        end
      Impulse = [ Impulse
		  Respaux ];
    end
  end
    clear Resp;
    clear Respaux;
end


% PLOTTING THE FUNCTIONS

irf = [Result Impulse];
for n = 1 : n_endog,
  for k = 1 : k_noises,
    lab = ['Impulse Response Function of ' 'endogenous' num2str(n) ' vs. noise' num2str(k)];
    label = [ label
	      lab ];
  end
end
for i = 1 : n_endog*k_noises,
  figure; whitebg('w');
  plot(irf(:,i));
  title(label(i,:));
  grid;
end

if size(irf,2) >= n_endog*k_noises+1
  for n = 1 : n_endog,
    for l = 1 : l_exog,
      lab = ['Impulse Response Function of ' 'endogenous' num2str(n) ' vs. exogenous' num2str(l)];
      label = [ label
		lab ];
    end
  end
  for i = n_endog*k_noises+1 : size(irf,2),
    for l = 1 : n_endog*l_exog,
      figure; whitebg('w');
      plot(irf(:,i));
      title(label(l,:));
      grid;
    end
  end
end

irf