% Raster scan ISAM
addpath(genpath('I:\Intercambio_Informacion_EAFIT\matlab'));
cMap = (gray(512));
% cMap = viridis(512);
cMap = cMap(32:end, :);

LATEX_DEF = {'Interpreter', 'latex'};
set(0, 'defaultTextInterpreter', 'LaTex')
set(0, 'defaultAxesTickLabelInterpreter', 'LaTex')
set(groot,'defaultLegendInterpreter','latex');
set(0, 'DefaultAxesFontSize', 20)
% set(gca, 'TickLabelInterpreter', 'latex')
%% Parameters
% Tomogram number of points
nZ = 512; % axial, number of pixels per A-line, accounting for zero-padding
nX = 256; % fast scan axis, number of A-lines per B-scan
nY = 1; % slow scan axis, number of B-scans per tomogram
nK = 512; % Number of samples, <= nZ, difference is zero padding
xNyquistOversampling = 2;
nXOversampling = nX; % Number of Alines for oversampling PSF, <= nX, difference is zero padding

curFig = 10;

% Definition of sampling is backwards: one has a light source with a given
% bandwidth, one then finds which nK provides you with the zRange desired.
% zScale is just zRange / nK.

useGPU = true;

% Spectral parameters
wavelength = 1.310e-6; % Source central wavelength
axialRes = 5e-6; 8e-6; % Axial resolution. Typical system 8e-6, it was tested with 3e-6

% Confocal parameters
numAper = 0.25; % clcNumerical aperture

if useGPU
  varType = {'single', 'gpuArray'};
  ThisLinspace = @gpuArray.linspace;
  ToSingle = @(x) gpuArray(single(x));
else
  varType = {'single'};
  ThisLinspace = @linspace;
  ToSingle = @single;
end


% Wavelength spectral width of the source
wavelengthWidthSource = 2 * log(2) / pi * wavelength ^ 2 / axialRes;
% Wavenumber spectral width of the source
wavenumberWidthSource = 2 * pi / (wavelength - (wavelengthWidthSource / 2)) - ...
  2 * pi / (wavelength + (wavelengthWidthSource / 2)); 

% Wavelength spectral width, increase 1.0 for oversampling
wavelengthWidthFull = 1.5 * pi / 2 / log(2) * wavelengthWidthSource;
% Wavenumber spectral range
wavenumberRange = 2 * pi ./ (wavelength + ([wavelengthWidthFull -wavelengthWidthFull] / 2));
% Wavenumber sampling vector. Because we are simulating complex fringes, we
% need nZ and not 2 * nZ
zeroPadding = (nZ - nK) / 2;
kVect = ThisLinspace(wavenumberRange(1), wavenumberRange(2), nK)';
wavenumber = single((wavenumberRange(1) + wavenumberRange(2)) / 2);

% Linear in wavenumber spectrum of the source
wavenumberFWHMSource = wavenumberWidthSource / (2 * sqrt(2 * log(2)));
sourceSpec = exp(-(kVect - wavenumber) .^ 2 / 2 / wavenumberFWHMSource ^ 2);
% figure(1), plot(wavenumberVect, sourceSpec), axis tight

% Physical size of axial axis
zSize = 1 * pi * nK / diff(wavenumberRange);
axSampling = zSize / nZ; % Axial sampling

% Confocal parameters
alpha = pi / numAper; % Confocal constant
beamWaistDiam = 2 * alpha / wavenumber; % Beam waist diameter 
latSampling = beamWaistDiam / sqrt(2) / xNyquistOversampling; % Latereral sampling
confocalParm = pi * (beamWaistDiam / sqrt(2)) ^ 2 / wavelength; % Confocal paremeter (for info.)

% Zero-path delay. Changing this changes the focal plane in the HighNA
% model
zRef = 0; - zSize / 4;
% Distance from the top plane to the focal plane. This is independent from
% zRef only for the LowNA models and not for the HighNA model
focalPlane = zRef;

xSize = latSampling * nX; % Physical size of fast scan axis
ySize = latSampling * nY; % Physical size of slow scan axis

% Coordinates
% Cartesian coordinate
zVect = single(ThisLinspace(-zSize / 2, zSize / 2 - axSampling, nZ));
xVect = single(ThisLinspace(-xSize / 2, xSize / 2 - latSampling, nX));
yVect = single(0); ThisLinspace(- ySize / 2, ySize / 2 - latSampling, nY);

% Frequency coordinates
freqBWFac = 2; % Increase frequency bandwdith to avoid artifact in numerical FT
nFreqX = nX * freqBWFac;
freqXVect = single(ThisLinspace(-0.5, 0.5 - 1 / nFreqX, nFreqX)) /...
  (latSampling / freqBWFac) * 2 * pi;

% Show High-NA Gaussiam beam
gaussianBeam1 = 1 / (2 * pi) * fftshift(fft(fftshift(1 ./ ...
  ((alpha ./ wavenumber) .^ 2 + (1i * (zVect' - focalPlane) ./ wavenumber)) .* ...
  exp(2i * (zVect' - zRef) .* kVect) .* ...
  exp(-1i * (zVect' - focalPlane) .* freqXVect .^ 2 ./ wavenumber / 4) .* ...
  exp(- (freqXVect * alpha ./ wavenumber / 2) .^ 2), 2), [], 2), 2);
figH1 = figure(curFig); subplot(1, 2, 1), imagesc(xVect * 1e6, zVect * 1e6, abs(gaussianBeam1)),
xlabel('$x$ [$\mu$m]'), title('b) Gaussian beam'), %ylabel('$z$ [$\mu$m]'), 
axis image, colormap(cMap), set(gca, 'YTick', []);

%% Create object with point scatteres
% Number of point scatterers
nPointSource = 128;
rng('default')

% Range where point scatteres appear
objZRange = nZ - 100;
objXRange = nX - 100;
% Because of complex fringes, z range is defined as [-nZ / 2, nZ / 2] * axSampling
% objPosZ = permute(linspace(- objZRange / 2, objZRange / 2, nPointSource) * axSampling, [1 3 2]);
objPosZ = (objZRange * (rand(1, 1, nPointSource, varType{:})) - objZRange / 2) * axSampling;
objPosX = (objXRange * (rand(1, 1, nPointSource, varType{:})) - objXRange / 2) * latSampling;
% Amplitude of scatterers
objAmp = ones(1, 1, nPointSource, varType{:}); % All ones by default

% Just to visualize point scatterers, offset to nZ / 2 to index matrix
% correctly due to complex fringes. Coerce is important because sometimes
% we will get an object that rounds to 0 index.
objSuscep = zeros(nZ, nX, nY, varType{:});
for k = 1:nPointSource
  objSuscep(Coerce(round(nZ / 2 + objPosZ(k) / axSampling), 1, nZ),...
    Coerce(round(objPosX(k) / latSampling) + nX / 2, 1, nX), :) = objAmp(k);
end

%% Forward Model
tic
% High NA Model
fringes1 = zeros(nK, nX, varType{:});
for thisScan = 1:nX
  % Current beam position
  thisBeamPosX = xVect(thisScan);
  % Spectrum at this beam possition is the contribution of the Gaussian
  % beam at the location of the point sources
  fringes1(:, thisScan) =  ForwardModel_PointScatterers_HighNA(objAmp, objPosZ,...
    objPosX - thisBeamPosX, kVect, freqXVect,...
    alpha, focalPlane, zRef);
end
toc
% Calculate frignes with proper constants, including source spectrum
fringes1 = fringes1 .* 1i ./ (2 * pi) .* 1 .* sqrt(sourceSpec) ./ kVect;

tic
% High NA Model
fringes2 = zeros(nK, nX, varType{:});
for thisScan = 1:nX
  % Current beam position
  thisBeamPosX = xVect(thisScan);
  % Spectrum at this beam possition is the contribution of the Gaussian
  % beam at the location of the point sources
  fringes2(:, thisScan) =  ForwardModel_PointScatterers_FreqLowNA(objAmp, objPosZ,...
    objPosX - thisBeamPosX, kVect, wavenumber, freqXVect,...
    alpha, focalPlane, zRef);
end
toc
% Calculate frignes with proper constants, including source spectrum
fringes2 = fringes2 .* 1i ./ (2 * pi) .* 1 .* sqrt(sourceSpec) ./ kVect;

%%
xNyquistOversampling2 = 2.5 * xNyquistOversampling;

% Confocal parameters
numAper = 0.1; % clcNumerical aperture
% Confocal parameters
alpha = pi / numAper; % Confocal constant
beamWaistDiam = 2 * alpha / wavenumber; % Beam waist diameter 
latSampling = beamWaistDiam / sqrt(2) / xNyquistOversampling2; % Latereral sampling
confocalParm = pi * (beamWaistDiam / sqrt(2)) ^ 2 / wavelength; % Confocal paremeter (for info.)

xSize = latSampling * nX; % Physical size of fast scan axis
ySize = latSampling * nY; % Physical size of slow scan axis

% Coordinates
% Cartesian coordinate
zVect = single(ThisLinspace(-zSize / 2, zSize / 2 - axSampling, nZ));
xVect = single(ThisLinspace(-xSize / 2, xSize / 2 - latSampling, nX));
yVect = single(0); ThisLinspace(- ySize / 2, ySize / 2 - latSampling, nY);

% Frequency coordinates
nFreqX = nX * freqBWFac;
freqXVect = single(ThisLinspace(-0.5, 0.5 - 1 / nFreqX, nFreqX)) /...
  (latSampling / freqBWFac) * 2 * pi;

% zRef = -100e-6;
% focalPlane = 0;
% Show Low-NA Gaussiam beam
gaussianBeam2 = 1 / (2 * pi) * fftshift(fft(fftshift(1 ./ ...
  ((alpha ./ wavenumber) .^ 2 + (1i * (zVect' - focalPlane) ./ wavenumber)) .* ...
  exp(2i * (zVect' - zRef) .* kVect) .* ...
  exp(-1i * (zVect' - focalPlane) .* freqXVect .^ 2 ./ wavenumber / 4) .* ...
  exp(- (freqXVect * alpha ./ wavenumber / 2) .^ 2), 2), [], 2), 2);

figH1 = figure(curFig); subplot(1, 2, 2), imagesc(xVect * 1e6, zVect * 1e6, abs(gaussianBeam2)), 
xlabel('$x$ [$\mu$m]'), title('b) Gaussian beam'), %ylabel('$z$ [$\mu$m]'), 
axis image, colormap(cMap), set(gca, 'YTick', []);

%% Forward Model
tic
% High NA Model
fringes3 = zeros(nK, nX, varType{:});
for thisScan = 1:nX
  % Current beam position
  thisBeamPosX = xVect(thisScan);
  % Spectrum at this beam possition is the contribution of the Gaussian
  % beam at the location of the point sources
  fringes3(:, thisScan) =  ForwardModel_PointScatterers_HighNA(objAmp, objPosZ,...
    objPosX - thisBeamPosX, kVect, freqXVect,...
    alpha, focalPlane, zRef);
end
toc
% Calculate frignes with proper constants, including source spectrum
fringes3 = fringes3 .* 1i ./ (2 * pi) .* 1 .* sqrt(sourceSpec) ./ kVect;

%% Fringes processing
% Fourier transform fringes to get tomogram
% Apply Hanning to avoid artifacts
tom1 = fftshift(fft(fftshift(padarray(fringes1 .* hanning(nK), zeroPadding, 'both'), 1), [], 1), 1);
tom2 = fftshift(fft(fftshift(padarray(fringes2 .* hanning(nK), zeroPadding, 'both'), 1), [], 1), 1);
tom3 = fftshift(fft(fftshift(padarray(fringes3 .* hanning(nK), zeroPadding, 'both'), 1), [], 1), 1);

% Aline are shifted depending on the zero-path delay
refShift = round((zRef - zSize) / zSize * nZ);
tom1 = circshift(tom1, [refShift 0]);
tom2 = circshift(tom2, [refShift 0]);
tom3 = circshift(tom3, [refShift 0]);

%% Show OCT images
% Plot optioncs
logLim = [65 inf];

% Figure 1
figH1 = figure(curFig + 1); subplot(1,3,1), imagesc(squeeze(xVect) * 1e6, squeeze(flip(zVect))* 1e6,...
  objSuscep(:, :, 1)), axis image
xlabel('$x$ [$\mu$m]'), ylabel('$z$ [$\mu$m]'), title('a) Sample'),
set(gca,'YDir','normal'), colormap(cMap), drawnow,
hCB = colorbar; hCB.TickLabelInterpreter = 'latex'; hCB.Ticks = linspace(0, 1, 6);
hCB.Label.String = '[a.u.]'; hCB.Label.Interpreter = 'latex'; hCB.Label.FontSize = 20;

subplot(1,3,2), imagesc(xVect * 1e6, zVect * 1e6, abs(gaussianBeam1) / max(abs(gaussianBeam1(:)))), 
xlabel('$x$ [$\mu$m]'), title('b) Gaussian beam'), %ylabel('$z$ [$\mu$m]'), 
axis image, colormap(cMap), set(gca, 'YTick', []);
hCB = colorbar; hCB.TickLabelInterpreter = 'latex'; hCB.Ticks = linspace(0, 1, 6);
hCB.Label.String = '[a.u.]'; hCB.Label.Interpreter = 'latex'; hCB.Label.FontSize = 20;

subplot(1,3,3), imagesc(squeeze(xVect) * 1e6, squeeze(zVect)* 1e6, 10 * log10(abs(tom1) .^ 2), logLim),
xlabel('$x$ [$\mu$m]'), title('c) OCT image'), axis image % ylabel('$z$ [$\mu$m]'), 
set(gca, 'YTick', [])
hCB = colorbar; hCB.TickLabelInterpreter = 'latex'; hCB.Label.String = 'log. scale [a.u.]';
hCB.Label.Interpreter = 'latex';  hCB.Label.FontSize = 20;

logLim = [65 110];

% Figure 2
figH2 = figure(curFig + 2); subplot(1,3,1), imagesc(squeeze(xVect) * 1e6, squeeze(flip(zVect))* 1e6,...
  objSuscep(:, :, 1)), axis image
xlabel('$x$ [$\mu$m]'), ylabel('$z$ [$\mu$m]'), title('a) Sample'),
set(gca,'YDir','normal'), colormap(cMap), drawnow,
hCB = colorbar; hCB.TickLabelInterpreter = 'latex'; hCB.Ticks = linspace(0, 1, 6);
hCB.Label.String = '[a.u.]'; hCB.Label.Interpreter = 'latex'; hCB.Label.FontSize = 20;

subplot(1,3,2), imagesc(xVect * 1e6, zVect * 1e6, abs(gaussianBeam2) / max(abs(gaussianBeam2(:)))), 
xlabel('$x$ [$\mu$m]'), title('b) Gaussian beam'), %ylabel('$z$ [$\mu$m]'), 
axis image, colormap(cMap), set(gca, 'YTick', []);
hCB = colorbar; hCB.TickLabelInterpreter = 'latex'; hCB.Ticks = linspace(0, 1, 6);
hCB.Label.String = '[a.u.]'; hCB.Label.Interpreter = 'latex'; hCB.Label.FontSize = 20;

subplot(1,3,3), imagesc(squeeze(xVect) * 1e6, squeeze(zVect)* 1e6, 10 * log10(abs(tom3) .^ 2), logLim),
xlabel('$x$ [$\mu$m]'), title('c) OCT image'), axis image % ylabel('$z$ [$\mu$m]'), 
set(gca, 'YTick', [])
hCB = colorbar; hCB.TickLabelInterpreter = 'latex'; hCB.Label.String = 'log. scale [a.u.]';
hCB.Label.Interpreter = 'latex';  hCB.Label.FontSize = 20;

% Figure 3
figH3 = figure(curFig + 3); subplot(1,3,1), imagesc(squeeze(xVect) * 1e6, squeeze(flip(zVect))* 1e6,...
  objSuscep(:, :, 1)), axis image
xlabel('$x$ [$\mu$m]'), ylabel('$z$ [$\mu$m]'), title('a) Sample'),
set(gca,'YDir','normal'), colormap(cMap), drawnow,
hCB = colorbar; hCB.TickLabelInterpreter = 'latex'; hCB.Ticks = linspace(0, 1, 6);
hCB.Label.String = '[a.u.]'; hCB.Label.Interpreter = 'latex'; hCB.Label.FontSize = 20;

subplot(1,3,2), imagesc(squeeze(xVect) * 1e6, squeeze(zVect)* 1e6, 10 * log10(abs(tom3) .^ 2), logLim),
xlabel('$x$ [$\mu$m]'), title('c) OCT image'), axis image % ylabel('$z$ [$\mu$m]'), 
set(gca, 'YTick', [])
hCB = colorbar; hCB.TickLabelInterpreter = 'latex'; hCB.Label.String = 'log. scale [a.u.]';
hCB.Label.Interpreter = 'latex';  hCB.Label.FontSize = 20;

subplot(1,3,3), imagesc(squeeze(xVect) * 1e6, squeeze(zVect)* 1e6, 10 * log10(abs(tom3Refocused) .^ 2), logLim),
xlabel('$x$ [$\mu$m]'), title('c) Defocus correction'), axis image % ylabel('$z$ [$\mu$m]'), 
set(gca, 'YTick', [])
hCB = colorbar; hCB.TickLabelInterpreter = 'latex'; hCB.Label.String = 'log. scale [a.u.]';
hCB.Label.Interpreter = 'latex';  hCB.Label.FontSize = 20;

if 0
  %% Forward model with motion
  motionUm = 1 * smoothdata(rand(1, nX) - 0.5) * axSampling;
  tic
  % High NA Model
  fringes4 = zeros(nK, nX, varType{:});
  for thisScan = 1:nX
    % Current beam position
    thisBeamPosX = xVect(thisScan);
    % Spectrum at this beam possition is the contribution of the Gaussian
    % beam at the location of the point sources
    fringes4(:, thisScan) =  ForwardModel_PointScatterers_HighNA(objAmp, objPosZ + motionUm(thisScan),...
      objPosX - thisBeamPosX, kVect, freqXVect,...
      alpha, focalPlane, zRef);
  end
  toc
  % Calculate frignes with proper constants, including source spectrum
  fringes4 = fringes4 .* 1i ./ (2 * pi) .* 1 .* sqrt(sourceSpec) ./ kVect;
  tom4 = fftshift(fft(fftshift(padarray(fringes4 .* hanning(nK), zeroPadding, 'both'), 1), [], 1), 1);
  tom4 = circshift(tom4, [refShift 0]);
end

%% High magnitude motion
rng('default')
motionUm = 1.5 * smoothdata(mean(rand(10, nX), 1) - 0.5, 'gaussian', 7) * 12.3;
tom3Mot = ApplyAxialShift(tom3, motionUm);
tom3MotPhase = tom3Mot .* exp(2i * wavenumber * motionUm * axSampling);
% tom3Mot = tom4;
% tom3MotPhase = tom4;
phaseThresDb = 65;

% Figure 4a
figH4a = figure(curFig + 4); plot(xVect * 1e6, motionUm * 1e6 * axSampling, 'r', 'linewidth', 1), axis tight
ylabel('$\delta z_n$ [$\mu$m]'), xlabel('$x$ [$\mu$m]'), grid minor
title('d) Induce motion')

% Refocus
tom3MotionRefocused = RefocusTomogramWithModifiedAmpFilter1D(tom3Mot,...
  axSampling * 1e6, latSampling * 1e6 * nX / nXOversampling, wavelength * 1e6, z0Px, 2, ApplyOptFilt);
tom3MotionPhaseRefocused = RefocusTomogramWithModifiedAmpFilter1D(tom3MotPhase,...
  axSampling * 1e6, latSampling * 1e6 * nX / nXOversampling, wavelength * 1e6, z0Px, 2, ApplyOptFilt);

% Figure 4b
figH4b = figure(curFig + 5); subplot(1, 3, 1), imagesc(squeeze(xVect)* 1e6, squeeze(flip(zVect))* 1e6, ...
  10 * log10(abs(tom3) .^ 2), logLim), axis image, set(gca,'YDir','normal'), 
xlabel('$x$ [$\mu$m]'), ylabel('$z$ [$\mu$m]'), title('a) Intenisity with no motion'),
hCB = colorbar; hCB.TickLabelInterpreter = 'latex'; hCB.Label.String = 'log. scale [a.u.]';
hCB.Label.Interpreter = 'latex';  hCB.Label.FontSize = 20;

subplot(1, 3, 2), imagesc(squeeze(xVect)* 1e6, squeeze(zVect)* 1e6, ...
  10 * log10(abs(tom3Mot) .^ 2), logLim), axis image
xlabel('$x$ [$\mu$m]'), title('b) Intenisity with motion'), %ylabel('$z$ [$\mu$m]'), 
set(gca, 'YTick', [])
hCB = colorbar; hCB.TickLabelInterpreter = 'latex'; hCB.Label.String = 'log. scale [a.u.]';
hCB.Label.Interpreter = 'latex';  hCB.Label.FontSize = 20;

subplot(1, 3, 3), imagesc(squeeze(xVect)* 1e6, squeeze(zVect)* 1e6, ...
  10 * log10(abs(tom3MotionPhaseRefocused) .^ 2), logLim), axis image
xlabel('$x$ [$\mu$m]'), title('c) Refocused intensity with motion'), %ylabel('$z$ [$\mu$m]'), 
set(gca, 'YTick', []), colormap(cMap)
hCB = colorbar; hCB.TickLabelInterpreter = 'latex'; hCB.Label.String = 'log. scale [a.u.]';
hCB.Label.Interpreter = 'latex';  hCB.Label.FontSize = 20;

% Figure 4c
maskInt = single(20 * log10(abs(tom3)) > phaseThresDb);
maskInt(~maskInt) = nan;
maskIntMotion = single(20 * log10(abs(tom3Mot)) > phaseThresDb);
maskIntMotion(~maskIntMotion) = nan;
figH4c = figure(curFig + 6); subplot(1, 3, 1), imagescnan(squeeze(xVect)* 1e6, squeeze(flip(zVect))* 1e6, ...
  angle(tom3 .* maskInt), [-pi pi]), axis image, set(gca,'YDir','normal'), 
xlabel('$x$ [$\mu$m]'), ylabel('$z$ [$\mu$m]'), title('e) Phase with no motion')
hCB = colorbar; hCB.TickLabelInterpreter = 'latex'; hCB.Label.String = '[rad.]';
hCB.Label.Interpreter = 'latex';  hCB.Label.FontSize = 20; hCB.Ticks = linspace(-pi,pi,3);
hCB.TickLabels = {'$-\pi$', '$0$', '$\pi$'};

subplot(1, 3, 2), imagescnan(squeeze(xVect)* 1e6, squeeze(zVect)* 1e6, ...
  angle(tom3MotPhase .* maskIntMotion), [-pi pi]), axis image, set(gca, 'YTick', []),
xlabel('$x$ [$\mu$m]'), title('f) Phase with motion'), %ylabel('$z$ [$\mu$m]'), 
colormap(cmap('c3', 'shift', -0.25))
hCB = colorbar; hCB.TickLabelInterpreter = 'latex'; hCB.Label.String = '[rad.]';
hCB.Label.Interpreter = 'latex';  hCB.Label.FontSize = 20; hCB.Ticks = linspace(-pi,pi,3);
hCB.TickLabels = {'$-\pi$', '$0$', '$\pi$'};

subplot(1, 3, 3), imagescnan(squeeze(xVect)* 1e6, squeeze(zVect)* 1e6, ...
  angle(tom3Mot .* maskIntMotion), [-pi pi]), axis image, set(gca, 'YTick', []),
xlabel('$x$ [$\mu$m]'), title({'g) Phase with motion', 'ignoring phase offsets'}), %ylabel('$z$ [$\mu$m]'), 
colormap(cmap('c3', 'shift', -0.25))
hCB = colorbar; hCB.TickLabelInterpreter = 'latex'; hCB.Label.String = '[rad.]';
hCB.Label.Interpreter = 'latex';  hCB.Label.FontSize = 20; hCB.Ticks = linspace(-pi,pi,3);
hCB.TickLabels = {'$-\pi$', '$0$', '$\pi$'};

%%
saveFig = 1;
if saveFig
  saveas(figH1, 'FM-HighNA.svg')
  saveas(figH2, 'FM-LowNA.svg')
  saveas(figH3, 'IM-Refocus.svg')
end

saveFig = 1;
if saveFig
  saveas(figH4a, 'MotionMag.svg')
  saveas(figH4b, 'FM-OCTImagesMotion.svg')
  saveas(figH4c, 'FM-OCTPhaseMotion.svg')
end