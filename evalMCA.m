function [B,Vc,Rc,Lc,lc,varargout] = evalMCA(h,r,q,Icoil,varargin)
%evalMCA takes in the parameters and returns the values of interest for
% constraint evaluation and optimization
% INPUTS:
%          h = 3x1 or 1x3 vector containing h0, h1, and h2 [m]
%          r = 7x1 or 1x7 vector containing r0 - r6 [m]
%          q = [m]
%      Icoil = coil current [A]
%   varargin = {1} headless (true or false) or (1 or 0)
% OUTPUT:
%          B = 1x5 array of Magnetic Flux Densities (Br,Bz) [T]
%         Vc = Voltage Drop, [V]
%         Rc = Resistance, [Ohms]
%         Lc = Inductance [H]
%         lc = length of Wire [m]
%  varargout = Variable outputs
%            - {1} mu = 1x5 array of Magnetic Permeabilities [H/m]
%            - {2} Phi = Flux [Wb]
%            - {3} Nturns = Number of Turns
Rperm = 0.1037/0.3048;
mu0 = 4*pi*1e-7;
murFe = 5000;
murMR = 5;
muFe = mu0*murFe;
muMR = mu0*murMR;

% Extract parameters
h0 = h(:,1); h1 = h(:,2); h2 = h(:,3);
r0 = r(:,1); r1 = r(:,2); r2 = r(:,3); 
r3 = r(:,4); r4 = r(:,5); r5 = r(:,6); r6 = r(:,7);

% Calculate Coil Geometry (input is in millimeters)
[~,Nturns,lc,RZ] = calcCoil(h0*1e3,r2*1e3,q*1e3);

% Calculate Coil Resistance and Voltage drop
Rc = Rperm*lc;
Vc = Rc.*Icoil;

% MCA %%%%%%%%

% Calculate lengths
l1 = h0 + h1;
l2 = r3-(r1+r0)/2;
l3 = r4-r3;
l4 = r5-r4;
l5 = (r6 - r5)/2;
l6 = h0 + h2;

% Calculate areas
A1 = pi*(r1.^2-r0.^2);
A2 = 2*pi*h1.*(2*r3+r1+r0)/4;
A3 = pi*h2.*(r4+r3);
A4 = pi*h2.*(r5+r4);
A5 = 2*pi*h2.*(r6+3*r5)/4;
A6 = pi*(r6.^2-r5.^2);

% Solve for Reluctance
Rmag = l1./(muFe.*A1) + (2*l2)./(muFe.*A2) + (2*l3)./(muFe.*A3) + (2*l4)./(muMR.*A4) + (2*l5)./(muFe.*A5) + l6./(muFe.*A6);

% Calculate magnetic mmf
Fmmf = Nturns.*Icoil;

% Solve for  B fields using Biot-Savart Field intensities (to avoid single layer ampere assumption)
%coilR = RZ(1);
%coilZ = RZ(2);

%Jcoil = (Nturns*Icoil)/(coilR*coilZ);

% [Hcoil,HcoilR,HcoilZ] = biotRingH(Jcoil,r1,coilZ,coilR);
% 
% Bcoil = Hcoil*muFe;
% BcoilR = HcoilR*muFe;
% BcoilZ = HcoilZ*muFe;

% Solve for Magnetic Flux
Phi = Fmmf./Rmag;

B1 = Phi./A1; 
B2 = Phi./A2; 
B3 = Phi./A3; 
B4 = Phi./A4;
B5 = Phi./A5; 
B6 = Phi./A6; 

B = [B1, B2, B3, B4, B5, B6];

% Solve for Inductance
% NI = Rmag*Phi -> Phi = NI/Rmag
% L = N*Phi/I -> L = N^2/Rmag
Lc = (Nturns.^2)./Rmag;

if nargout > 5
    mu = ones(size(B)).*[muFe,muFe,muFe,muMR,muFe,muFe];
    varargout{1} = mu;
    varargout{2} = Phi;
    varargout{3} = Nturns;
end

end