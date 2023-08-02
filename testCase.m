Tests various cases to validate model
clear; clc; close all;
%%
h0 = 0.5;
h1 = 0.5;
z2 = 0.5;
w1 = 0.6;
w2 = 0.25;
w3 = 0.25;
w4 = 0;
w5 = 0.15;
q = 0.5;
I = 0.5;
xTestN = [h0,h1,z2,w1,w2,w3,w4,w5,q,I]
[h,r,q,I] = denormalizeDesign(xTestN,false)

[F1,c1,xD1,R1] = objfunMCA(xTestN)

[F2,c2,xD2,R2] = objfunFEMM(xTestN,true)

%%
clear; clc;
mu0 = 4e-7*pi;
murFe = 5000;
I = 8;%0.5;
N = 2000;%1600;

h0 = 100e-3;%20e-3;
d = 1e-4;%2.5e-4;

n = h0/d; % Number of coils per row
Bm = N/n;  % Number of rows

q = Bm*d;

r0 = 0e-3;
r1 = 39.95e-3;

core = [ 
    r0,0;
    r1,0;
    r1,h0;
    r0,h0];
coil = [
      r1,0;
    r1+q,0;
    r1+q,h0;
      r1,h0];

W = q;
H = h0;
A = W*H;
J = N*I/A;
RC = r1;
LC = h0;
TC = W;

[Bp,Br,Bz,Bm,BavgR,BavgZ,B0,B1] = biotRingH(J,RC,LC,TC,1e-3,murFe);


figure(1)
t = tiledlayout(1,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
patch(core(:,2),core(:,1),[0.5 0.5 0.5])
hold on
patch(coil(:,2),coil(:,1),[0.7 0.4 0.2])
patch(core(:,2),-core(:,1),[0.5 0.5 0.5])
patch(coil(:,2),-coil(:,1),[0.7 0.4 0.2])
xlabel("z, m")
ylabel("r, m")
legend(["Core","Coil"],"Location","best")
daspect([1 1 1])
axis tight

nexttile
contX = [        B1(:,:,1);  B0(:,:,1); B1(:,:,1)];
contY = [-flipud(B1(:,:,2)); B0(:,:,2); B1(:,:,2)];
contC = [ flipud(B1(:,:,3)); B0(:,:,3); B1(:,:,3)];
coreContour = imagesc(contC*(4e-7*pi));
set(gca,'YTick',[],'XTick',[])
cb = colorbar;
ylabel(cb,"B, T")
daspect([1,1,1])
colormap parula

Bbiot = mean(B0(1,50:51,3))
Bamp = mu0*N*I/h0
%% 
% Open FEMM
openfemm();

main_maximize();

% Create New Document
newdocument(0);

% Temporarily save
mi_saveas('testCase.FEM')

% Problem Setup
mi_probdef(0,'meters','axi',1e-8,0,30,0);
mi_smartmesh(1);

% Load Materials
mi_getmaterial('Air');

mi_getmaterial('Pure Iron');
mi_modifymaterial('Pure Iron',0,'Iron Core');

mi_getmaterial('0.25mm');
mi_modifymaterial('0.25mm',0,'0.25mm Coil');

coref = core;
coref(:,2) = core(:,2) - h0/2;

coilf = coil;
coilf(:,2) = coil(:,2) - h0/2;

% Define core
mi_addnode(coref(1,1), coref(1,2))
for i = 2:length(coref)
    mi_addnode(coref(i,1), coref(i,2));
    mi_addsegment(coref(i-1,1), coref(i-1,2), coref(i,1), coref(i,2));
end
mi_addsegment(coref(end,1), coref(end,2), coref(1,1), coref(1,2));
mi_addblocklabel(mean(coref(:,1)),mean(coref(:,2)));
mi_selectlabel(mean(coref(:,1)),mean(coref(:,2)));
mi_setblockprop('Iron Core',1,0,'<None>',0,0,1);%'Air',1,0,'<None>',0,0,1);
mi_clearselected();

% Define Coil
mi_addcircprop('Coil',I,1)

mi_addnode(coilf(1,1),coilf(1,2));
for i = 2:length(coilf)
    mi_addnode(coilf(i,1),coilf(i,2));
    mi_addsegment(coilf(i-1,1),coilf(i-1,2),coilf(i,1),coilf(i,2));
end
mi_addsegment(coilf(end,1),coilf(end,2),coilf(1,1),coilf(1,2));
mi_addblocklabel(mean(coilf(:,1)),mean(coilf(:,2)));
mi_selectlabel(mean(coilf(:,1)),mean(coilf(:,2)));
mi_setblockprop('0.25mm Coil',1,0,'Coil',0,0,N)
mi_clearselected();

% Add air
mi_addblocklabel(1.5*h0,0);
mi_selectlabel(1.5*h0,0);
mi_setblockprop('Air',1,0,'<None>',0,0,1);
mi_clearselected();

% Add boundary
mi_makeABC(7,3*h0,0,0,0);

% TempSave
mi_zoom(0,-h0,h0,h0)

% Create Mesh
mi_createmesh();

% Perform Analysis
mi_analyze();

% Load Solution
mi_loadsolution();

% Get B
Bfemm = mo_getb(0,0);
Bfemm = Bfemm(2)