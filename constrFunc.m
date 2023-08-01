function [Y] = constrFunc(B,Foff,Fon,tauLR,Vc)
%constrFunc Takes in the nonlinear constraints needed for optimization and
% then checks them against specified limits, returning an array with the
% associated "costs"
%   INPUTS:
%        B = 6x1 Magnetic flux densities [T]
%     Fmin = Force when magnet off [N]
%     Fmax = Force when magnet on [N]
%    tauLR = Circuit time constant [s]
%       Vc = Voltage drop over coil [V]
%   OUTPUT:
%        c = array of costs
BminFe = 0; BmaxFe = 2.2; % Target 1.6 or less 
BminMR = 0; BmaxMR = 1; % Target 0.6 or less
Foffmin = 0; Foffmax = 5;
Fonmin = 40; Fonmax = 120;
taumin = 0; taumax = 0.1;
Vmin = 0; Vmax = 13;
Vnom = 12;

% Constraint Evaluation
g1 = BminFe - B(1);
g2 = B(1)/BmaxFe - 1;
g3 = BminFe - B(2);
g4 = B(2)/BmaxFe - 1;
g5 = BminFe - B(3);
g6 = B(3)/BmaxFe - 1;
g7 = BminMR - B(4);
g8 = B(4)/BmaxMR - 1;
g9 = BminFe - B(5);
g10 = B(5)/BmaxFe - 1;
g11 = BminFe - B(6);
g12 = B(6)/BmaxFe - 1;
g13 = Foffmin - Foff;
g14 = Foff/Foffmax - 1;
g15 = 1 - Fon/Fonmin;
g16 = Fon/Fonmax - 1;
g17 = taumin - tauLR;
g18 = tauLR/taumax - 1;
g19 = Vmin - Vc;
if Vc > Vnom
    g20 = Vc/Vmax - 1;
else
    g20 = 0;
end

Y = zeros(1,20);

% UPPER LAGRANGE MULTIPLIERS
% Magnetic Fields
Y(1) = max(0,g1);
Y(3) = max(0,g3);
Y(5) = max(0,g5);
Y(7) = max(0,g7);
Y(9) = max(0,g9);
Y(11) = max(0,g11);

% Forces
Y(13) = max(0,g13);
Y(15) = max(0,g15);

% Time constant
Y(17) = max(0,g17);

% Voltage
Y(19) = max(0,g19);

% LOWER LAGRANGE MULTIPLIERS
% Magnetic Fields
Y(2) = max(0,g2);
Y(4) = max(0,g4);
Y(6) = max(0,g6);
Y(8) = max(0,g8);
Y(10) = max(0,g10);
Y(12) = max(0,g12);

% Forces
Y(14) = max(0,g14);
Y(16) = max(0,g16);

% Time constant
Y(18) = max(0,g18);

% Voltage
Y(20) = max(0,g20);

end