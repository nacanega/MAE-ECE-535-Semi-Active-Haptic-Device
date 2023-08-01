function Force = MRforce(h0,h2,r2,r4,r5,varargin)
%MRforce takes in the r4, r5, h2, and any optional B applied and returns
% the force on the spool due to the MR fluid
% INPUTS:
%      h0 - h0 [m]
%      h2 - h2 [m]
%      r2 - r2 [m]
%      r4 - r4 [m]
%      r5 - r5 [m]
% varagin - {1} B [T] if present
% OUTPUT:
%   Force - the force in Newtons on the spool
load("handoutData.mat","gamDotEta","Btau")

v = 1;

gammaDot1 = v./(r4.*log(r5./r4));
eta1 = interp1(gamDotEta(:,1),gamDotEta(:,2),gammaDot1,"linear","extrap");
T1 = eta1.*gammaDot1;

gammaDot2 = v./(r2.*log(r5./r2));
eta2 = interp1(gamDotEta(:,1),gamDotEta(:,2),gammaDot2,"linear","extrap");
T2 = eta2.*gammaDot2;

if nargin == 6
    Bmr = varargin{1};
    T1 = T1 + 1000*interp1(Btau(:,1),Btau(:,2),Bmr,"linear");
end

A1 = 4*pi*r4*h2;
A2 = 2*pi*r2*h0;

Force = T1*A1 + T2*A2;

end