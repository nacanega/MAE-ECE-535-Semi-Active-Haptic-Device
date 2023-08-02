function [F,varargout] = objfunFEMM(x,varargin)
%objfunFEMM Is the multi-objective function defining the project which
% takes in a normalized design string and returns the objectives at that
% location as well as constraint values if needed
% INPUTS:
%         x = 10x1 normalized design string
%  varargin = {1} showFEMM (true/1 or false/0)
% OUTPUT:
%         F = 1x3 vector of objective function values
% varargout = additional outputs
%           - {1} 2x constraint value matrix
%                 [B1-B6 Fmin Fmax tauLR Vc]
%                 row 1 shows value of constraining terms
%                 row 2 shows whether lower constraint is active
%                 row 3 shows whether upper constraint is active
%           - {2} designVector [h,r,q,Ic]
%           - {3} Additonal terms [Rfix,Mass1,Mass2,Pc,Vwaste];
K = 20;
Vnom = 12; % Nominal Voltage

if nargin == 2 && varargin{1}
    hideFEMM = false;
else
    hideFEMM = true;
end

% Denormalize in mm
[h,r,q,Ic] = denormalizeDesign(x,false);

% Evaluate using FEMM
[B,Vc,Rc,Lc,lc,mu,Phi,Nturns] = evalFEMM(h,r,q,Ic,hideFEMM);

% Convert lengths to meters
h = h*1e-3; r = r*1e-3; q = q*1e-3;

% Extract dVs
h0 = h(:,1); h1 = h(:,2); h2 = h(:,3);
r0 = r(:,1); r1 = r(:,2); r2 = r(:,3); 
r3 = r(:,4); r4 = r(:,5); r5 = r(:,6); r6 = r(:,7);

% Find Forces
Fmin = MRforce(h0,h2,r2,r4,r5);
Fmax = MRforce(h0,h2,r2,r4,r5,B(:,4));

% Adds a series resistor to get 12V total potential difference if needed
if Vc < Vnom
    Rfix = (Vnom - Vc)./Ic;
else
    Rfix = 0;
end

% Calculate circuit time constant
tauLR = Lc./(Rc+Rfix);

% Determine lagrange multipliers for constraints
lambda = constrFunc(B,Fmin,Fmax,tauLR,Vc);

% Penalty term
penaltyTerm = sum((K*lambda).^2);

% Objectives

% Mass
[massSpool, massTube] = calcMass(lc,h0,h1,h2,r0,r1,r3,r4,r5,r6);
% Power
Pc = Vc.*Ic;
% Wasted Volume
Vwaste = 2*pi*(r4-q/2) .* (r4-r2).*h0;

F(1) = massSpool + penaltyTerm;
F(2) = Pc + penaltyTerm;
F(3) = -Fmax + penaltyTerm;

if nargout > 1
    calcParams = [B,Fmin,Fmax,tauLR,Vc];
    cMat = [calcParams;lambda(1:2:19);lambda(2:2:20)];
    varargout{1} = cMat;
    designVector = [h,r,q,Ic];
    varargout{2} = designVector;
    varargout{3} = [Rfix,massSpool,massTube,Pc,Vwaste,mu,Phi,Nturns];
end

end