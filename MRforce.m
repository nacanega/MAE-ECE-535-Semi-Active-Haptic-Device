function Force = MRforce(h2,r4,r5,varargin)
%MRforce takes in the r4, r5, h2, and any optional B applied and returns
% the force on the spool due to the MR fluid
% INPUTS:
%      h2 - h2 [m]
%      r4 - r4 [m]
%      r5 - r5 [m]
% varagin - {1} B [T] if present
% OUTPUT:
%   Force - the force in Newtons on the spool
load("handoutData.mat","gamDotEta","Btau")

if nargin == 4
    Bmr = varargin{1};
    T1 = 1000*interp1(Btau(:,1),Btau(:,2),Bmr,"linear");
else
    v = 1;
    gammaDot = v./(0.5*(r4+r5).*log(r5./r4));
    eta = interp1(gamDotEta(:,1),gamDotEta(:,2),gammaDot,"linear","extrap");
    T1 = eta.*v./(r5-r4);
end

A = pi*(r4+r5)*h2;

Force = T1*A;

end