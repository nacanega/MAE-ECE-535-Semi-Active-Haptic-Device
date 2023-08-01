function [coilPts, Nturns, lenWire, varargout] = calcCoil(h0,r1,q,varargin)
%calcCoil Calculates the four corners of the coil rectangle for FEMM as
% well as the number of turns and length of wire
% INPUTS:
%       h0 = h0 [mm]
%       r1 = r1 [mm]
%       q = q [mm]
% varargin = {2} diameter of magnet wire [mm]
% OUTPUT:
%  coilPts = 4x2 (x,y) points of coil [mm]
%   Nturns = Number of turns of wire [Turns]
%  lenWire = Length of wire [m]
% varargout = {1} Magnetic field intensity and average intensity from Biot-Savart

if nargin == 4
    d = varargin{1};
else
    d = 0.254;
end

s = sqrt(3)/2;

w = q;

if w > d
    m = floor((w-d)/(d*s)) + 1;
else
    m = 1;
end

n = floor(h0/d);
boxW = d*(1+s*(m-1));

Nturns = n*m;

if abs(mod(h0,d) - d/2) > eps(d) 
    flag = true;
    boxH = d*n;
    Nturns = Nturns - floor(m/2);
else
    boxH = d*(n + 0.5);
    flag = false;
end

coilPts = [
           r1, -boxH/2;
    r1 + boxW, -boxH/2;
    r1 + boxW,  boxH/2;
           r1,  boxH/2];

if nargout == 4
    varargout{1} = 1e-3*[boxW,boxH];
end

r = r1 + d/2;
lenWire = 2*pi*r*n;

for i = 2:m
    r = r + s*d;
    if flag && mod(i,2) == 0
        lenWire = lenWire + 2*pi*r*(n-1);
    else
        lenWire = lenWire + 2*pi*r*n;
    end
end

lenWire = lenWire/1000; % Convert to meters