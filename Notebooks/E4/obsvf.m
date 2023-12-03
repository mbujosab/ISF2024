function [abar,bbar,cbar,t,k] = obsvf(a,b,c,tol)
%OBSVF	Observability staircase form.
%	[ABAR,BBAR,CBAR,T,K] = OBSVF(A,B,C) returns a decomposition
%	into the observable/unobservable subspaces.
%	[ABAR,BBAR,CBAR,T,K] = OBSVF(A,B,C,TOL) uses tolerance TOL.
%
%	If Ob=OBSV(A,C) has rank r <= n, then there is a similarity
%	transformation T such that
%
%      Abar = T * A * T' ,  Bbar = T * B  ,  Cbar = C * T' .
%
%	and the transformed system has the form
%
%	      | Ano   A12|           |Bno|
%	Abar =  ----------  ,  Bbar =  ---  ,  Cbar = [ 0 | Co].
%	      |  0    Ao |           |Bo |
%
%	                                           -1           -1
%	where (Ao,Bo) is controllable, and Co(sI-Ao) Bo = C(sI-A) B.
%
%	See also: OBSV.

%	Author : R.Y. Chiang  3-21-86
%	Revised 5-27-86 JNL
%	Copyright (c) 1986-93 by the MathWorks, Inc.

% Use CTRBF and duality:

if nargin == 4
	[aa,bb,cc,t,k] = ctrbf(a',c',b',tol);
else
	[aa,bb,cc,t,k] = ctrbf(a',c',b');
end
abar = aa'; bbar = cc'; cbar = bb';
