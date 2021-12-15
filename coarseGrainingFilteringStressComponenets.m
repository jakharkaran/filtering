close all;
clear;

addpath('../eqsdiscovery/');

dataType = 'Re20kNX1024nx4ny0r0p1';
% Re20kNX1024nx4ny0r0p1 Re20kNX1024nx4ny0r0p01 
% Re20kNX1024nx25ny25r0p1 Re100kNX2048nx4ny0r0p1

Lx = 2*pi;
Ngrid = 64;
Nfilter = 'default';
countSnap = 1;

if strcmp(dataType,'Re20kNX1024nx4ny0r0p1') || strcmp(dataType,'Re20kNX1024nx4ny0r0p01') || ...
        strcmp(dataType,'Re20kNX1024nx25ny25r0p1')
    N_DNS = 1024;
elseif strcmp(dataType,'Re100kNX2048nx4ny0r0p1')
    N_DNS = 2048;
end

load(['../data/2DTurbulence/' dataType '/DNS/train1/DNS1.mat'])

N_LES = [64,64];
Delta = 2*Lx/N_LES(1);

%% DNS grid
Lx = 2*pi;
dx = Lx/N_DNS;
x = linspace(0,Lx-dx,N_DNS);
kx = (2*pi/Lx)*[0:(N_DNS/2) (-N_DNS/2+1):-1];

[Ky,Kx] = meshgrid(kx,kx);

U_DNS = real(ifft2(1i*Ky.*fft2(slnPsiDNS(:,:,countSnap))));
V_DNS = -real(ifft2(1i*Kx.*fft2(slnPsiDNS(:,:,countSnap))));

plotStress = ['S12'];

filterTypes = ["gaussian", "boxSpectral", "boxPhysical", "gaussian+boxSpectral", ...
    "gaussian+boxPhysical", "spectral+boxPhysical", "spectral"];
coarseGrainingType = ["spectral", "spectral", "physical", "spectral", ...
    "physical", "physical", "spectral"];

% h1 = figure(1);python 2DTurbulenceUV.py S Re20kNX1024nx4ny0r0p1 train1 boxSpectral 1024 128 4 50 1e-6 1 1 0 ./
h2 = figure(2);

counterFilter = 0;
for plotCount = 1:28

    if plotCount == 1 || plotCount == 5 || plotCount == 9 || plotCount == 13 || ...
        plotCount == 17 || plotCount == 21 || plotCount == 25 || plotCount == 29

    counterFilter = counterFilter + 1;
%     counterFilter = 1;

    [S11_residual, S12_residual, S22_residual, ...
        S11_leonard, S12_leonard, S22_leonard, ...
        S11_cross, S12_cross, S22_cross,...
        S11_reynolds, S12_reynolds, S22_reynolds,...
        S11_reynolds_1, S12_reynolds_1, S22_reynolds_1,...
        S11_reynolds_2, S12_reynolds_2, S22_reynolds_2] = residualStressComponents2D( ...
            U_DNS,V_DNS, filterTypes(counterFilter), ... 
            coarseGrainingType(counterFilter), Delta, N_LES);

    [filterTypes(counterFilter) coarseGrainingType(counterFilter) 'Residual' min(S11_residual(:)) max(S11_residual(:))]
    [filterTypes(counterFilter) coarseGrainingType(counterFilter) 'Leonard' min(S11_leonard(:)) max(S11_leonard(:))]
    [filterTypes(counterFilter) coarseGrainingType(counterFilter) 'Cross' min(S11_cross(:)) max(S11_cross(:))]
    [filterTypes(counterFilter) coarseGrainingType(counterFilter) 'Reynolds' min(S11_reynolds(:)) max(S11_reynolds(:))]
    [filterTypes(counterFilter) coarseGrainingType(counterFilter) 'Reynolds' min(S11_reynolds_1(:)) max(S11_reynolds_1(:))]
    [filterTypes(counterFilter) coarseGrainingType(counterFilter) 'Reynolds' min(S11_reynolds_2(:)) max(S11_reynolds_2(:))]
    
    if strcmp(plotStress,'S11')
        S_residual = S11_residual;
        S_leonard = S11_leonard;
        S_cross = S11_cross;
        S_reynolds = S11_reynolds;
        S_reynolds_1 = S11_reynolds_1;
        S_reynolds_2 = S11_reynolds_2;
        titlesS = ["\tau_{11} Residual", "\tau_{11} Leonard", "\tau_{11} Cross", ...
            "\tau_{11} Reynolds"];
        titlesS2 = ["\tau_{11} = \bar{u'_{i}u'_{j}}", "\tau_{11} = \bar{u'_{i}}\bar{u'_{j}}"];
    
    elseif strcmp(plotStress,'S12')
        S_residual = S12_residual;
        S_leonard = S12_leonard;
        S_cross = S12_cross;
        S_reynolds = S12_reynolds;
        S_reynolds_1 = S12_reynolds_1;
        S_reynolds_2 = S12_reynolds_2;
        titlesS = ["\tau_{12} Residual", "\tau_{12} Leonard", "\tau_{12} Cross", ...
            "\tau_{12} Reynolds"];
    
    elseif strcmp(plotStress,'S22')
        S_residual = S22_residual;
        S_leonard = S22_leonard;
        S_cross = S22_cross;
        S_reynolds = S22_reynolds;
        S_reynolds_1 = S22_reynolds_1;
        S_reynolds_2 = S22_reynolds_2;
        titlesS = ["\tau_{22} Residual", "\tau_{22} Leonard", "\tau_{22} Cross", ...
            "\tau_{22} Reynolds"];
    end

    [corr2(S_residual,S_leonard+S_cross+S_reynolds), ...
        RMSE2(S_residual,S_leonard+S_cross+S_reynolds)];
    end

    if mod(plotCount,4) == 1
        S = S_residual;
        countTitle = 1;
        subplotCount = counterFilter;

        if N_LES(1) == 128
            limS = 0.2;
        elseif N_LES(1) == 64
            limS = 0.6;
        end
    
    if strcmp(filterTypes(counterFilter),"spectral") && strcmp(coarseGrainingType(counterFilter),"spectral")
        if N_LES(1) == 128
            limS = 0.02;
        elseif N_LES(1) == 64
            limS = 0.06;
        end
    end

    elseif mod(plotCount,4) == 2
        S = S_leonard;
        countTitle = 2;
        subplotCount = 7 + counterFilter;

        if N_LES(1) == 128
            limS = 0.2;
        elseif N_LES(1) == 64
            limS = 0.6;
        end

    if strcmp(filterTypes(counterFilter),"spectral") && strcmp(coarseGrainingType(counterFilter),"spectral")
        if N_LES(1) == 128
            limS = 0.02;
        elseif N_LES(1) == 64
            limS = 0.06;
        end
    end

    elseif mod(plotCount,4) == 3
        S = S_cross;
        countTitle = 3;
        subplotCount = 14 + counterFilter;

        if N_LES(1) == 128
            limS = 0.02;
        elseif N_LES(1) == 64
            limS = 0.15;
        end

    if strcmp(filterTypes(counterFilter),"spectral") && strcmp(coarseGrainingType(counterFilter),"spectral")
        if N_LES(1) == 128
            limS = 0.02;
        elseif N_LES(1) == 64
            limS = 0.06;
        end
    end

    elseif mod(plotCount,4) == 0
        S = S_reynolds;
        countTitle = 4;
        subplotCount = 21 + counterFilter;

        if N_LES(1) == 128
            limS = 0.005;
        elseif N_LES(1) == 64
            limS = 0.02;
        end

    if strcmp(filterTypes(counterFilter),"spectral") && strcmp(coarseGrainingType(counterFilter),"spectral")
        if N_LES(1) == 128
            limS = 0.0006;
        elseif N_LES(1) == 64
            limS = 0.005;
        end
    end
    end
    
    hold on;

%     
%     subplot(4,7,subplotCount);
%     
%     s = pcolor(S);
%     mapS1 =  b2r(-limS,limS);
%     colormap(mapS1); colorbar; caxis([-limS limS]);
%     axis equal; s.EdgeColor = 'none'; 
%     axis([0 N_LES(1) 0 N_LES(2)]);
%     title(titlesS(countTitle));
%     set(gca,'XTickLabel', [], 'YTickLabel', []);
%     hold off;

    if mod(plotCount,4) == 1 || mod(plotCount,4) == 2

        if mod(plotCount,4) == 1 

            S = S_reynolds_1;
            countTitle = 4;
            subplotCount = counterFilter;
    
            if N_LES(1) == 128
                limS = 0.005;
            elseif N_LES(1) == 64
                limS = 0.02;
            end
    
            if strcmp(filterTypes(counterFilter),"spectral") && strcmp(coarseGrainingType(counterFilter),"spectral")
                if N_LES(1) == 128
                    limS = 0.0006;
                elseif N_LES(1) == 64
                    limS = 0.005;
                end
            end

        elseif mod(plotCount,4) == 2

            S = S_reynolds_2;
            countTitle = 4;
            subplotCount = 7 + counterFilter;
    
            if N_LES(1) == 128
                limS = 0.005;
            elseif N_LES(1) == 64
                limS = 0.02;
            end
    
            if strcmp(filterTypes(counterFilter),"spectral") && strcmp(coarseGrainingType(counterFilter),"spectral")
                if N_LES(1) == 128
                    limS = 0.0006;
                elseif N_LES(1) == 64
                    limS = 0.005;
                end
            end
        end
        
    subplot(2,7,subplotCount);
    
    s = pcolor(S);
    mapS1 =  b2r(-limS,limS);
    colormap(mapS1); colorbar; caxis([-limS limS]);
    axis equal; s.EdgeColor = 'none'; 
    axis([0 N_LES(1) 0 N_LES(2)]);
%     title(titlesS2(countTitle));
    set(gca,'XTickLabel', [], 'YTickLabel', []);
    hold off;

    end


end

saveas(h2, '2.png')

% [min(S11_residual(:)) max(S11_residual(:))]
% [min(S11_leonard(:)) max(S11_leonard(:))]
% [min(S11_cross(:)) max(S11_cross(:))]
% [min(S11_reynolds(:)) max(S11_reynolds(:))]
% 
% [min(S12_residual(:)) max(S12_residual(:))]
% [min(S12_leonard(:)) max(S12_leonard(:))]
% [min(S12_cross(:)) max(S12_cross(:))]
% [min(S12_reynolds(:)) max(S12_reynolds(:))]
% 
% [min(S22_residual(:)) max(S22_residual(:))]
% [min(S22_leonard(:)) max(S22_leonard(:))]
% [min(S22_cross(:)) max(S22_cross(:))]
% [min(S22_reynolds(:)) max(S22_reynolds(:))]
