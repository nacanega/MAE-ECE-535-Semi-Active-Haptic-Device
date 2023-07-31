function [H,HR,HZ,varargout] = biotRingH(J,RC,LC,TC)
%biotRingH takes in the current density and coil dimensions and returns the
% estimated axial field intensity at the middle of the ring
% INPUTS:
%     J = Coil Current Density [A/m]
%    RC = Coil Inner Radius (core radius) [m]
%    LC = Coil Length (along z) [m]
%    TC = Coil Thickness (along radius) [m]
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
M = ceil(TC/2.54e-4);
N = ceil(LC/2.54e-4);
O = ceil(RC/2.54e-4);

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

Hfun = @(theta,Rc,Rp,Zp) ...
    ((J.*Rc) ./ (4*pi)) .* (Rc - Rp.*cos(theta)) ...
    ./ ((Zp.^2 + Rc.^2 + Rp.^2 - 2.*Rc.*Rp.*cos(theta)).^(1.5));

% Point
Pr = 0; % Center of coil
Hfun1 = @(theta) Hfun(theta,coilR,Pr,coilZ-LC/2);
H = sum(integral(Hfun1,0,2*pi,"ArrayValued",true),"all")*dAcoil;

% for Pz = ZSteps - LC/2
%     for Rc = RcoilSteps
%         Hfun1 = @(theta) Hfun(theta,Rc,Pr,Pz);
%         dHPoint = integral(Hfun1,0,2*pi);
%         H = H + dHPoint*dAcoil;
%     end
% end

if nargout > 1
    % Along Z
    HavgZ = zeros(1,N);
    for i = 1:N
        Hfun2 = @(theta) Hfun(theta,coilR,Pr,coilZ-ZSteps(i));
        HavgZ(i) = sum(integral(Hfun2,0,2*pi,"ArrayValued",true),"all")*dAcoil;
        %     for Pz = ZSteps - ZSteps(i)
        %         for Rc = RcoilSteps
        %             Hfun1 = @(theta) Hfun(theta,Rc,Pr,Pz);
        %             dHPoint = integral(Hfun1,0,2*pi);
        %             HavgZ(i) = HavgZ(i) + dHPoint*dAcoil;
        %         end
        %     end
    end
    HZ = mean(HavgZ);
end

if nargout > 2
    % Along R
    Pr = RcoreSteps;
    HavgR = zeros(1,O);
    for i = 1:O
        Hfun3 = @(theta) Hfun(theta,coilR,Pr(i),coilZ-LC/2);
        HavgR(i) = sum(integral(Hfun3,0,2*pi,"ArrayValued",true),"all")*dAcoil;
        %     for Pz = ZSteps-LC/2
        %         for Rc = RcoilSteps
        %             Hfun1 = @(theta) Hfun(theta,Rc,Pr(i),Pz);
        %             dHPoint = integral(Hfun1,0,2*pi);
        %             HavgR(i) = HavgR(i) + dHPoint*dAcoil;
        %         end
        %     end
    end
    HR = mean(HavgR);
end

if nargout > 3
    % All
    [coreZ,coreR] = meshgrid(ZSteps,RcoreSteps);
    Havg = zeros(O,N);
    for i = 1:O
        for j = 1:N
            Hfun4 = @(theta) Hfun(theta,coilR,Pr(i),coilZ-ZSteps(j));
            Havg(i,j) = sum(integral(Hfun4,0,2*pi,"ArrayValued",true),"all")*dAcoil;
            %         for Pz = ZSteps - ZSteps(j)
            %             for Rc = RcoilSteps
            %                 Hfun1 = @(theta) Hfun(theta,Rc,Pr(i),Pz);
            %                 dHPoint = integral(Hfun1,0,2*pi);
            %                 Havg(i,j) = Havg(i,j) + dHPoint*dAcoil;
            %             end
            %         end
        end
        varargout{1} = mean(Havg,"all");
    end
end

if nargout > 4
    varargout{2} = HavgR;
    if nargout > 5
        varargout{3} = cat(3,ZSteps,zeros(size(ZSteps)),HavgZ);
        if nargout > 6
            varargout{4} = cat(3,coreZ,coreR,Havg);
        end
    end
end

end