% Tests various cases to validate model
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