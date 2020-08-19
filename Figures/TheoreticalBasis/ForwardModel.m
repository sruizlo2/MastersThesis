% Raster scan ISAM
addpath(genpath('../../../../matlab'));

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
wavelength = 1.285e-6; % Source central wavelength
axialRes = 6e-6; 8e-6; % Axial resolution. Typical system 8e-6, it was tested with 3e-6

% Confocal parameters
numAper = 0.3; % clcNumerical aperture

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
wavelengthWidthFull = 1 * pi / 2 / log(2) * wavelengthWidthSource;
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
% figure2(1), plot(wavenumberVect, sourceSpec), axis tight

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
zRef = 0; - 1 * zSize;
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
freqBWFac = 1; % Increase frequency bandwdith to avoid artifact in numerical FT
nFreqX = nX * freqBWFac;
freqXVect = single(ThisLinspace(-0.5, 0.5 - 1 / nFreqX, nFreqX)) /...
  (latSampling / freqBWFac) * 2 * pi;

%% Show Gaussiam beam
gaussianBeam = 1 / (2 * pi) * fftshift(fft(fftshift( ...
  exp(1i * (zVect' - zRef) .* sqrt((2 * kVect) .^ 2 - freqXVect .^ 2)) .* ...
  exp(- (freqXVect * alpha ./ kVect / 2) .^ 2) * 0.5 ./ ...
  ((alpha ./ kVect) .^ 2 + (1i * (zVect' - zRef) ./ kVect)), 2), [], 2), 2); %
figure2(curFig), imagesc(abs(gaussianBeam)),

%% Create object with point scatteres
% Number of point scatterers
nPointSource = 25;

% Range where point scatteres appear
objZRange = nZ - 40;
objXRange = 1; nX - 40;
% Because of complex fringes, z range is defined as [-nZ / 2, nZ / 2] * axSampling
objPosZ = permute(linspace(- objZRange / 2, objZRange / 2, nPointSource) * axSampling, [1 3 2]);
% objPosZ = (objZRange * (rand(1, 1, nPointSource, varType{:})) - objZRange / 2) * axSampling;
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

% Show point scatterers
figure2(curFig + 1), imagesc(objSuscep(:, :, 1)), axis image
xlabel('x'), ylabel('z'), title('Sample'), colorbar

cMap = viridis(256);
colormap(cMap), drawnow

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

%%
% Low NA Model in frequency domain
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

% Low NA Model in Spatial domain
fringes3 = zeros(nK, nX, varType{:});
for thisScan = 1:nX
  % Current beam position
  thisBeamPosX = xVect(thisScan);
  % Spectrum at this beam possition is the contribution of the Gaussian
  % beam at the location of the point sources
  fringes3(:, thisScan) =  ForwardModel_PointScatterers_LowNA(objAmp, objPosZ,...
    objPosX - thisBeamPosX, kVect, wavenumber,...
    alpha, focalPlane, zRef);
end

%% Fringes processing
% Calculate frignes with proper constants, including source spectrum
fringes1 = fringes1 .* 1i ./ (2 * pi) .* 1 .* sqrt(sourceSpec) ./ kVect;
fringes2 = fringes2 .* 1i ./ (2 * pi) .* 1 .* sqrt(sourceSpec) ./ kVect;
fringes3 = fringes3 .* 1i ./ (2 * pi) .* 1 .* sqrt(sourceSpec) ./ kVect;

% Normalize LowNA2 fringes because its scale is very different
fringes3 = fringes3 ./ max(abs(fringes3(:))) * max(abs(fringes2(:)));

% Fourier transform fringes to get tomogram
% Apply Hanning to avoid artifacts
tom1 = fftshift(fft(fftshift(padarray(fringes1 .* hanning(nK), zeroPadding, 'both'), 1), [], 1), 1);
tom2 = fftshift(fft(fftshift(padarray(fringes2 .* hanning(nK), zeroPadding, 'both'), 1), [], 1), 1);
tom3 = fftshift(fft(fftshift(padarray(fringes3 .* hanning(nK), zeroPadding, 'both'), 1), [], 1), 1);

% Add lateral oversampling
tom1 = fftshift(ifft(fftshift(tom1, 2), [], 2), 2);
tom1 = fftshift(fft(fftshift(padarray(tom1, [0 (nXOversampling - nX) / 2], 'both'), 2), [], 2), 2);

% Aline are shifted depending on the zero-path delay
refShift = round((zRef - zSize) / zSize * nZ);
tom1 = circshift(tom1, [refShift 0]);
tom2 = circshift(tom2, [refShift 0]);
tom3 = circshift(tom3, [refShift 0]);

%% Show OCT images
% Plot optioncs
logLimX = [50 inf];
logLim = [50 inf];

figure2(curFig + 3), imagesc(squeeze(xVect) * 1e6, squeeze(zVect)* 1e6,...
  10 * log10(abs(tom1) .^ 2), logLimX), axis image
colorbar, xlabel('x'), ylabel('z'), title('Beam with low/high-NA definition in \xi')

figure2(curFig + 4), imagesc(squeeze(xVect) * 1e6, squeeze(zVect)* 1e6,...
  10 * log10(abs(tom2) .^ 2), logLim), axis image
colorbar, xlabel('x'), ylabel('z'), title('Beam with low-NA definition in \xi')

figure2(curFig + 5), imagesc(squeeze(xVect) * 1e6, squeeze(zVect)* 1e6,...
  10 * log10(abs(tom3) .^ 2), logLim), axis image
colorbar, xlabel('x'), ylabel('z'), title('Beam with low-NA definition in X')

figure2(curFig + 6), subplot(1,3,1),
imagesc(squeeze(xVect) * 1e6, squeeze(zVect)* 1e6, 10 * log10(abs(tom1) .^ 2), logLimX),
colorbar, xlabel('x [\mum]'), ylabel('z [\mum]'), title('Beam with low/high-NA definition in \xi')

subplot(1,3,2),
imagesc(squeeze(xVect)* 1e6, squeeze(zVect)* 1e6, 10 * log10(abs(tom2) .^ 2), logLim),
colorbar, xlabel('x [\mum]'), ylabel('z [\mum]'), title('Beam with low-NA definition in \xi')

subplot(1,3,3),
imagesc(squeeze(xVect)* 1e6, squeeze(zVect)* 1e6, 10 * log10(abs(tom3) .^ 2), logLim),
colorbar, xlabel('x [\mum]'), ylabel('z [\mum]'), title('Beam with low-NA definition in X')

%% Apply CAO
z0Px = nZ / 2 + focalPlane / axSampling;
ApplyOptFilt = 0;
logLim = [70 140];
tom1Refocused = RefocusTomogramWithModifiedAmpFilter1D(flip(tom1, 3),...
  axSampling * 1e6, latSampling * 1e6 * nX / nXOversampling, wavelength * 1e6, z0Px, 2, ApplyOptFilt);

tom2Refocused = RefocusTomogramWithModifiedAmpFilter1D(flip(tom2, 3),...
  axSampling * 1e6, latSampling * 1e6 * nX / nXOversampling, wavelength * 1e6, z0Px, 2, ApplyOptFilt);

figure2(curFig + 7), imagesc(squeeze(xVect)* 1e6, squeeze(zVect)* 1e6, ...
  10 * log10(abs(tom1Refocused) .^ 2), logLim), axis image
colorbar, xlabel('x'), ylabel('z'), title('Refocused Beam with low/high-NA definition in \xi')

figure2(curFig + 8), hold off,
plot(squeeze(zVect)* 1e6, 10 * log10(abs(tom1Refocused(:, end/2)) .^ 2), 'k')
hold on, plot(squeeze(zVect)* 1e6, 10 * log10(abs(tom2Refocused(:, end/2)) .^ 2), 'r')
xlabel('x'), title('Refocused Beam with low/high-NA definition in \xi')

%% ISAM
% fringes1FT = fftshift(ifft(ifftshift(...
%   flip(gather(fringes1), 1), 2), [],  2), 2);
fringes1FT = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(...
  flip(gather(tom1), 1), 1), 2), [], 1), [], 2), 1), 2);
figure(100), imagesc(squeeze(20 * log10(abs(fringes1FT(:, :, 1)))))

kVect2 = gather(kVect(1:1:end));

freqZVect = gather(sqrt(4 * kVect2 .^ 2 - freqXVect .^ 2)) / 2;
figure2(100), plot( freqZVect(:,end/2), 'r'), hold on
plot(kVect2, 'k--'), hold off
figure2(101), imagesc(freqZVect)

tom1FT = zeros(nZ, nX, nY);
for thisFreqX = 1:nX
  for thisFreqY = 1:nY
    tom1FT(:, thisFreqX, thisFreqY) = interp1(kVect2, ...
      fringes1FT(:, thisFreqX, thisFreqY), freqZVect(:, thisFreqX, thisFreqY), 'pchip');
  end
end

% tom1FTRe = interp2(repmat(kVect2, [1 nX]), ...
%   real(fringes1FT), freqZVect, 'spline');
% tom1FTIm = interp2(repmat(kVect2, [1 nX]), ...
%   imag(fringes1FT), freqZVect, 'spline');
% tom1FT = complex(tom1FTRe, tom1FTIm);

tom1ISAM = ifftshift(ifftshift(fft(fft(fftshift(fftshift(...
  tom1FT, 1), 2), [], 1), [], 2), 1), 2);

figure(102), imagesc(squeeze(xVect)* 1e6, squeeze(zVect)* 1e6, ...
  squeeze(20 * log10(abs(tom1ISAM))))
axis image, colorbar

figure2(curFig + 9), hold off,
plot(squeeze(zVect)* 1e6, 10 * log10(abs(tom1Refocused(:, end/2)) .^ 2), 'k'), hold on
plot(squeeze(zVect)* 1e6, 10 * log10(abs(tom1ISAM(:, end/2)) .^ 2), 'r')
xlabel('x'), title('ISAM vs Refocused Beam with low/high-NA definition in \xi')
