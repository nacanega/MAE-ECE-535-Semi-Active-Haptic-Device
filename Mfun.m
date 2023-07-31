function [mass] = Massfun(lc,h0,h1,h2,r1,r3,r4)
rho_cu = 8960;          % kg/m^3
rho_fe = 7874;          % kg/m^3
mass = pi*rho_cu*0.25*d^2*lc + pi*rho_fe*(h0.*(r1.^2 - r0.^2) + 2*h1.*(r3.^2 - r0.^2) + 2*h2.*(r4.^2 - r3.^2));
end