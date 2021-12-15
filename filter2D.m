function [Uf_c] = filter2D(U_DNS,filterType,coarseGrainingType, Delta, N_LES)
%filter2D Filters and coarse grains 2D square grids

% filterType: gaussian
% 'gaussian': Gaussian Filter (Implemented in spectral domain)
% 'boxSpectral', 'boxPhysical', 'box': Box or Top hat filter 
%  (Implemented both in spectral and physical domain
% 'spectral' : Sharp Spectral filter (Implemented in spectral domain)
% 'gaussian+boxSpectral' or 'boxSpectral+gaussian': Gaussian and Box Filters
% 'gaussian+boxPhysical' or 'boxPhysical+gaussian': Gaussian and Box Filters
% 'spectral+boxPhysical' : Spectral followed by boxPhysical

% Coarse graining
% 'spectral' : Removing high wave-numbers in spectral domain
% 'physical' : subsampling in physical domain

addpath('../eqsdiscovery/');

N_DNS = size(U_DNS);

%% DNS grid 2D Turbulence
Lx = 2*pi;
dx = Lx/N_DNS(1);
% x = linspace(0,Lx-dx,N_DNS(1));
kx = (2*pi/Lx)*[0:(N_DNS(1)/2) (-N_DNS(1)/2+1):-1];

[Ky,Kx] = meshgrid(kx,kx);
Ksq = Kx.^2 + Ky.^2;
invKsq = 1./Ksq;
invKsq(1,1) = 0;
Kabs = sqrt(Ksq);

% Coarse grid
kkx = (2*pi/Lx)*[0:(N_LES(1)/2) (-N_LES(1)/2+1):-1];
[KKy,KKx] = meshgrid(kkx,kkx);
KKsq = KKx.^2 + KKy.^2;

if N_DNS == N_LES
    coarseGrainingType = 'none';
end


%% Filtering

if strcmp(filterType,'gaussian')
    Gk = exp(-Ksq*Delta^2/24);
    Uf = ifft2(Gk.*fft2(U_DNS));

elseif strcmp(filterType,'box') || strcmp(filterType,'boxSpectral')
%     Gk = sin(0.5 * sqrt(Ksq) * Delta)./(0.5*sqrt(Ksq)*Delta);
%     Gk(1,1) = 1;
    Gkx = sin(0.5 .* Kx * Delta)./(0.5.*Kx*Delta);
    Gky = sin(0.5 .* Ky * Delta)./(0.5.*Ky*Delta);
    Gkx(1,:) = 1;
    Gky(:,1) = 1;
    Gk = Gkx .* Gky;
% 
    Uf = ifft2(Gk.*fft2(U_DNS));

elseif strcmp(filterType,'boxPhysicalYifei')
    % Physical

    NX = N_DNS(1);
    NNX = N_LES(1);

    u = U_DNS;

    ratio = 2*NX/NNX;
    
    % Periodic extension
    uext = zeros(NX+ratio,NX+ratio);
    % middle points
    uext(0.5*ratio+1:0.5*ratio+NX,0.5*ratio+1:0.5*ratio+NX) = u;
    % up-down extensions
    for idx = 1:ratio/2
        uext(idx,0.5*ratio+1:0.5*ratio+NX) = u(NX-(ratio/2-idx),:);
        uext(NX+idx+0.5*ratio,0.5*ratio+1:0.5*ratio+NX) = u(idx,:);
    end
    % left-right extensions
    for idx = 1:ratio/2
        uext(:,idx) = uext(:,NX+idx);
        uext(:,NX+idx+0.5*ratio) = uext(:,idx+ratio*0.5);
    end
    
    % Box filtering
    uFp = zeros(NX,NX);
    for idx = ratio/2+1:ratio/2+NX
        for idy = ratio/2+1:ratio/2+NX
            box = uext(idx-ratio/2:idx+ratio/2,idy-ratio/2:idy+ratio/2);
            uFp(idx-ratio/2,idy-ratio/2) = mean(box(:));
        end
    end

    Uf = uFp;

elseif strcmp(filterType,'gaussian+box') || strcmp(filterType,'box+gaussian') || ...
        strcmp(filterType,'gaussian+boxSpectral') || strcmp(filterType,'boxSpectral+gaussian')
    Gkx = sin(0.5 .* Kx * Delta)./(0.5.*Kx*Delta);
    Gky = sin(0.5 .* Ky * Delta)./(0.5.*Ky*Delta);
    Gkx(1,:) = 1;
    Gky(:,1) = 1;
    Gkbox = Gkx .* Gky;

    Gkgaussian = exp(-Ksq*Delta^2/24);
    Gk = Gkgaussian .* Gkbox;
    Uf = ifft2(Gk.*fft2(U_DNS));

elseif strcmp(filterType,'boxPhysical') || strcmp(filterType,'boxPhysical+gaussian') || ...
        strcmp(filterType,'gaussian+boxPhysical') || strcmp(filterType,'spectral+boxPhysical')

    if strcmp(filterType,'boxPhysical+gaussian') || ...
        strcmp(filterType,'gaussian+boxPhysical')
        Gk = exp(-Ksq*Delta^2/24);
        U_DNS = ifft2(Gk.*fft2(U_DNS));

    elseif strcmp(filterType,'spectral+boxPhysical')
        U_DNS = spectralFilter_same_size(fft2(U_DNS),N_LES(1));
    end

    boxKernelSize = round(Delta/dx);

    U_FDNS = zeros(N_DNS(1)+boxKernelSize-1, N_DNS(2)+boxKernelSize-1);

    for countX = 1:N_DNS(1)+boxKernelSize-1
    
        xgrid = [countX-boxKernelSize:countX];
        for countGrid = 1:length(xgrid)
            if xgrid(countGrid) <= 0
                xgrid(countGrid) = xgrid(countGrid) + N_DNS(1);
            elseif xgrid(countGrid) > N_DNS(1)
                xgrid(countGrid) = xgrid(countGrid) - N_DNS(1);
            end
        end
    
        for countY = 1:N_DNS+boxKernelSize-1
            ygrid = [countY-boxKernelSize:countY];
            for countGrid = 1:length(ygrid)
                if ygrid(countGrid) <= 0
                    ygrid(countGrid) = ygrid(countGrid) + N_DNS(2);
                elseif ygrid(countGrid) > N_DNS(2)
                    ygrid(countGrid) = ygrid(countGrid) - N_DNS(2);
                end
            end
            U_FDNS(countX,countY) = mean(U_DNS(xgrid,ygrid), 'all');
        end     
    end

    temp1 = floor(boxKernelSize/2)+1;
    temp2 = boxKernelSize-temp1;
    sizeFDNS = size(U_FDNS);
    Uf = U_FDNS(temp1:(sizeFDNS(1)-temp2), temp1:(sizeFDNS(2)-temp2));

elseif strcmp(filterType,'spectral')
    Uf = spectralFilter_same_size(fft2(U_DNS),N_LES(1));

elseif strcmp(filterType,'none')
    Uf = U_DNS;
end


%% Coarse Graining
if strcmp(coarseGrainingType, 'spectral')

    Uf_c = spectralFilter(fft2(Uf),N_LES(1));

elseif strcmp(coarseGrainingType, 'physical')
% Subsampling in physical space

    xgrid = (1:N_LES(1))*(N_DNS(1)/N_LES(1)) - (N_DNS(1)/N_LES(1))/2;
    ygrid = (1:N_LES(2))*(N_DNS(2)/N_LES(2)) - (N_DNS(2)/N_LES(2))/2;

    Uf_c = Uf(xgrid,ygrid);

% elseif strcmp(coarseGrainingType, 'box')
% % Spacial averaging (Box filter + Subsampling)
% 
%     Uf_c = zeros(N_LES);
%     ratio = N_DNS/Nfilter;
%     
%     for countx = 1:N_LES(1)
%         for county = 1:N_LES(2)
%             xGrid = [(countx*ratio-(ratio-1)):(countx*ratio)];
%             yGrid = [(county*ratio-(ratio-1)):(county*ratio)];
%             Uf_c(countx,county) = mean(U_DNS(xGrid,yGrid),'all');
%         end
%     end

elseif strcmp(coarseGrainingType, 'none')
    Uf_c = Uf;
end

end