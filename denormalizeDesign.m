function [h,r,q,I] = denormalizeDesign(x,mScale)
%denormalizeDesign Takes in the normalized design string and outputs a
%standard converted string concatenated horizontally or split
% INPUTS:
%         x - (1x10) [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10]
%                    [h0 h1 z2 w1 w2 w3 w4 w5  q  Ic]
%    mScale - (1 or true) or (0 or false) convert to meters
% OUTPUT:
%      hrqI - (1x12) [h0 h1 h2 r0 r1 r2 r3 r4 r5 r6 q I]
%              OR (1x3) h = [h0 h1 h2]
% varargout - {1} (1x7) r = [r0 r1 r2 r3 r4 r5 r6]
%             {2} (1x1) q
%             {3} (1x1) I 

if mScale
    s = 1e-3;
else
    s = 1;
end

d = 0.254;
r0 = 1*s;

% We need to do this in two steps since q depends on w2 and w3
% xLB
xLB = [ 5, 1, 0, 1, 0, d, 1, 1, d, 0.05];
% xUB
xUB = [35,10,15,10,10,10, 3,10,20, 1];

x(1:8) = xLB(1:8) + (xUB(1:8)-xLB(1:8)).*x(1:8);
xUB(9) = x(5) + x(6);
x(9:10) = xLB(9:10) + (xUB(9:10)-xLB(9:10)).*x(9:10);
x(1:9) = x(1:9)*s;

% Preallocate
[m,~] = size(x);
r = zeros(m,7);

% Copy Values
h = [x(1:2), x(2) + x(3)];

r(1) = r0;           % r0
r(2) = r(1) + x(4);  % r1
r(3) = r(2);         % r2
r(4) = r(3) + x(5);  % r3
r(5) = r(4) + x(6);  % r4
r(6) = r(5) + x(7);  % r5
r(7) = r(6) + x(8);  % r6

q = x(9); 
I = x(10);

end