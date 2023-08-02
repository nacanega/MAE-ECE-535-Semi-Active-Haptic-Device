function [B,BR,BZ,varargout] = biotRingH(J,RC,LC,TC,varargin)
%biotRingH takes in the current density and coil dimensions and returns the
% estimated axial field intensity at the middle of the ring
% INPUTS:
%        J = Coil Current Density [A/m^2]
%       RC = Coil Inner Radius (core radius) [m]
%       LC = Coil Length (along z) [m]
%       TC = Coil Thickness (along radius) [m]
% varargin = {1} deltaRZ the step size 
%            {2} relative permeability of core (otherwise assumes air)
% OUTPUT:
%         H = Field intensity at center of coil
%        HR = Average field intensity along radius
%        HZ = Average field intensity along z
% varargout = Other outputs
%             {1} HM = Average field intensity over entire core
%             {2} HavgR = All evaluated field intensities along R at center
%             {3} HavgZ = All evaluated field intensities along Z at R = 0
%             {4} Coordinates and Havg along axis in 3D array
%             {5} Coordinates and Havg along upper half in 3D array
mu0 = 4e-7*pi;

if nargin == 5
    deltaRZ = varargin{1}; % for best performance, make at least twice the wire diameter
else
    deltaRZ = 1e-3; %1 mm increment size
end
M = ceil(TC/deltaRZ);
N = ceil(LC/deltaRZ);
O = ceil(RC/deltaRZ);

RcoilSteps = linspace(RC,RC+TC,M+1);
RcoilSteps = 0.5*(RcoilSteps(2:end) + RcoilSteps(1:end-1));
ZSteps = linspace(0,LC,N+1);
ZSteps = 0.5*(ZSteps(2:end) + ZSteps(1:end-1));
RcoreSteps = linspace(0,RC,O+1);
RcoreSteps = 0.5*(RcoreSteps(2:end) + RcoreSteps(1:end-1));

[coilZ,coilR] = meshgrid(ZSteps,RcoilSteps);

dRcoil = (TC)/(M);
dZ = (LC)/(N);
%dRcore = (RC)/(O);

dAcoil = dRcoil*dZ;
%dAcore = dRcore*dZ;

if nargin == 6
    mur = varargin{2};
    Bfun = @(theta,Rc,Rp,Zc,Zp) ...
    ((J.*Rc*mu0*mur) ./ (4*pi)) .* (Rc - Rp.*cos(theta)) ...
    ./ (((Zp-Zc).^2 + Rc.^2 + Rp.^2 - 2.*Rc.*Rp.*cos(theta)).^(1.5));
else
    Bfun = @(theta,Rc,Rp,Zc,Zp) ...
    ((J.*Rc*mu0) ./ (4*pi)) .* (Rc - Rp.*cos(theta)) ...
    ./ (((Zp-Zc).^2 + Rc.^2 + Rp.^2 - 2.*Rc.*Rp.*cos(theta)).^(1.5));
end

% Point
Pr = 0; % Center of coil
Bfun1 = @(theta) Bfun(theta,coilR,Pr,coilZ,LC/2);
B = sum(integral(Bfun1,0,2*pi,"ArrayValued",true),"all")*dAcoil;

% for Pz = ZSteps - LC/2
%     for Rc = RcoilSteps
%         Hfun1 = @(theta) Hfun(theta,Rc,Pr,Pz);
%         dHPoint = integral(Hfun1,0,2*pi);
%         H = H + dHPoint*dAcoil;
%     end
% end

if nargout > 1
    % Along Z
    BavgZ = zeros(1,N);
    for i = 1:N
        Hfun2 = @(theta) Bfun(theta,coilR,Pr,coilZ,ZSteps(i));
        BavgZ(i) = sum(integral(Hfun2,0,2*pi,"ArrayValued",true),"all")*dAcoil;
        %     for Pz = ZSteps - ZSteps(i)
        %         for Rc = RcoilSteps
        %             Hfun1 = @(theta) Hfun(theta,Rc,Pr,Pz);
        %             dHPoint = integral(Hfun1,0,2*pi);
        %             HavgZ(i) = HavgZ(i) + dHPoint*dAcoil;
        %         end
        %     end
    end
    BZ = mean(BavgZ);
end

if nargout > 2
    % Along R
    Pr = RcoreSteps;
    BavgR = zeros(1,O);
    for i = 1:O
        Bfun3 = @(theta) Bfun(theta,coilR,Pr(i),coilZ,LC/2);
        BavgR(i) = sum(integral(Bfun3,0,2*pi,"ArrayValued",true),"all")*dAcoil;
        %     for Pz = ZSteps-LC/2
        %         for Rc = RcoilSteps
        %             Hfun1 = @(theta) Hfun(theta,Rc,Pr(i),Pz);
        %             dHPoint = integral(Hfun1,0,2*pi);
        %             HavgR(i) = HavgR(i) + dHPoint*dAcoil;
        %         end
        %     end
    end
    BR = mean(BavgR);
end

if nargout > 3
    % All
    [coreZ,coreR] = meshgrid(ZSteps,RcoreSteps);
    Bavg = zeros(O,N);
    for i = 1:O
        for j = 1:N
            Bfun4 = @(theta) Bfun(theta,coilR,Pr(i),coilZ,ZSteps(j));
            Bavg(i,j) = sum(integral(Bfun4,0,2*pi,"ArrayValued",true),"all")*dAcoil;
            %         for Pz = ZSteps - ZSteps(j)
            %             for Rc = RcoilSteps
            %                 Hfun1 = @(theta) Hfun(theta,Rc,Pr(i),Pz);
            %                 dHPoint = integral(Hfun1,0,2*pi);
            %                 Havg(i,j) = Havg(i,j) + dHPoint*dAcoil;
            %             end
            %         end
        end
        varargout{1} = mean(Bavg,"all");
    end
end

if nargout > 4
    varargout{2} = BavgR;
    if nargout > 5
        varargout{3} = BavgZ;
        if nargout > 6
            varargout{4} = cat(3,ZSteps,zeros(size(ZSteps)),BavgZ);
            if nargout > 7
                varargout{5} = cat(3,coreZ,coreR,Bavg);
            end
        end
    end
end

end