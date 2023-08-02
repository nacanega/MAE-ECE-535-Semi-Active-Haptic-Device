function [B,Vc,Rc,Lc,lc,varargout] = evalFEMM(h,r,q,Icoil,varargin)
%evalFEMM takes in the parameters and returns the values of interest for
% constraint evaluation and optimization
% INPUTS:
%          h = 3x1 or 1x3 vector containing h0, h1, and h2 [mm]
%          r = 7x1 or 1x7 vector containing r0 - r6 [mm]
%          q = [mm]
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

if nargin == 5 && varargin{1}
    hideWindow = true;
else
    hideWindow = false;
end

% Extract parameters
h0 = h(1); h1 = h(2); h2 = h(3);
r0 = r(1); r1 = r(2); %r2 = r(3); 
r3 = r(4); r4 = r(5); r5 = r(6); r6 = r(7);

% Define points
spool = [
    r0, 0;
    r0, h0/2;
    r0, h0/2 + h1;
    r3, h0/2 + h1;
    r3, h0/2 + h2;
    r4, h0/2 + h2;
    r4, h0/2;
    r1, h0/2;
    r1, 0];

if h2 == h1
    spool(4,:) = [];
end

tube = [
    r5, -2*(h0+2*h2);
    r6, -2*(h0+2*h2);
    r6,  2*(h0+2*h2);
    r5,  2*(h0+2*h2)];

mrfluid = [
    r5, -2*(h0+2*h2);
     0, -2*(h0+2*h2);
     0,  2*(h0+2*h2);
    r5,  2*(h0+2*h2)];

spool = [spool; [1 -1].*flipud(spool(1:end-1,:))];

% Calculate Coil Geometry
[coil,Nturns,~] = calcCoil(h0,r1,q);

% Add FEMM to path
%addpath("c:\\femm42\\mfiles");

% Open FEMM
openfemm(hideWindow)

if ~hideWindow
    main_maximize();
end

% Create New Document
newdocument(0);

fileNum = 0;
s = sprintf('tempAnalysis%d.FEM',fileNum);
while isfile(s)
    fileNum = randi(1000);
    s = sprintf('tempAnalysis%d.FEM',fileNum);
end

% Temporarily save
mi_saveas(s)

% Problem Setup
mi_probdef(0,'millimeters','axi',1e-8,0,30,0);
mi_smartmesh(1);

% Load Materials
mi_getmaterial('Air');

mi_getmaterial('Pure Iron');
mi_modifymaterial('Pure Iron',0,'Iron Spool');
% mi_clearbhpoints('Iron Spool');
% mi_modifymaterial('Iron Spool',1,5000);
% mi_modifymaterial('Iron Spool',2,5000);

mi_getmaterial('Pure Iron');
mi_modifymaterial('Pure Iron',0,'Iron Tube');
% mi_clearbhpoints('Iron Tube');
% mi_modifymaterial('Iron Tube',1,5000);
% mi_modifymaterial('Iron Tube',2,5000);

mi_getmaterial('MRX-336AG');
% mi_clearbhpoints('MRX-336AG');
% mi_modifymaterial('MRX-336AG',1,5);
% mi_modifymaterial('MRX-336AG',2,5);

mi_getmaterial('30 AWG');
mi_modifymaterial('30 AWG',0,'30 AWG Coil');

% Define Geometry

% Define Spool
mi_addnode(spool(1,1), spool(1,2))
for i = 2:length(spool)
    mi_addnode(spool(i,1), spool(i,2));
    mi_addsegment(spool(i-1,1), spool(i-1,2), spool(i,1), spool(i,2));
end
mi_addsegment(spool(end,1), spool(end,2), spool(1,1), spool(1,2));
mi_addblocklabel((r0+r1)/2,(h0+h1)/2);
mi_selectlabel((r0+r1)/2,(h0+h1)/2);
mi_setblockprop('Iron Spool',1,0,'<None>',0,0,1);
mi_clearselected();

% Define Tube
mi_addnode(tube(1,1), tube(1,2))
for i = 2:length(tube)
    mi_addnode(tube(i,1), tube(i,2));
    mi_addsegment(tube(i-1,1), tube(i-1,2), tube(i,1), tube(i,2));
end
mi_addsegment(tube(end,1), tube(end,2), tube(1,1), tube(1,2));
mi_addblocklabel((r5+r6)/2,-(h0+h1)/2);
mi_selectlabel((r5+r6)/2,-(h0+h1)/2);
mi_setblockprop('Iron Tube',1,0,'<None>',0,0,1);
mi_clearselected();

% Define MR Fluid
for i = 2:length(mrfluid)-1
    mi_addnode(mrfluid(i,1), mrfluid(i,2));
    mi_addsegment(mrfluid(i,1),mrfluid(i,2),mrfluid(i-1,1),mrfluid(i-1,2));
end
mi_addsegment(mrfluid(i,1),mrfluid(i,2),mrfluid(end,1),mrfluid(end,2));
mi_addblocklabel(r0/2,(h0+2*h2));
mi_selectlabel(r0/2,(h0+2*h2));
mi_setblockprop('MRX-336AG',1,0,'<None>',0,0,1);
mi_clearselected();

% Define Coil
mi_addcircprop('Coil',Icoil,1)

mi_addnode(coil(1,1),coil(1,2));
for i = 2:length(coil)
    mi_addnode(coil(i,1),coil(i,2));
    mi_addsegment(coil(i-1,1),coil(i-1,2),coil(i,1),coil(i,2));
end
mi_addsegment(coil(end,1),coil(end,2),coil(1,1),coil(1,2));
mi_addblocklabel(mean(coil(:,1)),0);
mi_selectlabel(mean(coil(:,1)),0);
mi_setblockprop('30 AWG Coil',1,0,'Coil',0,0,Nturns)
mi_clearselected();

% Add air
mi_addblocklabel(r6+r1,h0+2*h2);
mi_selectlabel(r6+r1,h0+2*h2);
mi_setblockprop('Air',1,0,'<None>',0,0,1);
mi_clearselected();

% Add boundary
mi_makeABC(7,max(3*(h0+2*h2),1.5*r6),0,0,0);

% TempSave
mi_zoom(0,-h0-h2,2*r6,h0+h2)

% Create Mesh
%fprintf('Started Creating Mesh for %s\n',s)
mi_createmesh();
%fprintf('Finished Creating Mesh for %s\n',s)

% Perform Analysis
mi_analyze();

% Load Solution
mi_loadsolution();
if ~hideWindow
    mo_zoom(0,-h0-h2,2*r6,h0+h2);
    mo_showdensityplot(1,0,0,2.2,'mag');
end

% Get B Values of Centers Sections to check saturation
Bpoints = [
    (r1+r0)/2,           0;
    (r1+r3)/2,   (h0+h1)/2;
    (r4+r3)/2,   (h0+h2)/2;
    (r5+r4)/2,   (h0+h2)/2;
    (r6+3*r5)/4, (h0+h2)/2;
    (r6+r5)/2,           0];
Bvals = zeros(size(Bpoints));
muvals = Bvals;

for i = 1:length(Bpoints)
    Bvals(i,:) = mo_getb(Bpoints(i,1),Bpoints(i,2));
end

B = vecnorm(Bvals,2,2)';

circProps = mo_getcircuitproperties('Coil');
Vc = circProps(2);
Rc = Vc/Icoil;
lc = Rc/Rperm;
Phi = circProps(3);
Lc = Phi/Icoil;

if nargout > 5
    for i = 1:length(Bpoints)
        muvals(i,:) = mo_getmu(Bpoints(i,1),Bpoints(i,2));
    end
    mu = vecnorm(muvals,2,2)';
    varargout{1} = mu;
    varargout{2} = Phi;
    varargout{3} = Nturns;
end
    
if hideWindow
    mo_close();
    mi_close();
    delete(s);
    delete(sprintf('tempAnalysis%d.ans',fileNum));
    closefemm();
end

end