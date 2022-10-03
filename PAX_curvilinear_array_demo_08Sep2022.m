% PAX_curvilinear_array_demo_08Sep2022.m
% Matlab code for demonstrating the PAX method (Gray and Coussios 2022, IEEE UFFC) of sound speed estimation for 
% coregistration of B-Mode and PAM (passive acoustic mapping) images.
% The demonstration uses simulated backscatter data from a point source on the central axis of a C52v curvilinear array
% Signals were propagated through a two-layer model consisting of a fast (eg muscle/liver) ambient medium and a peripheral
% layer of slow (eg fat) material.

% Some implementation notes:
% 1. The PAM imaging domain (L50) is kept relatively narrow here to save time. Similarly, the Bmode calculations 
%    are only made on a single central A-line. Both of these simplifications were made for purpose of illustrating PAX with
%    an on-axis target. When operating on experiment data, a pre-processing step is performed to identify a region of interest 
%    This is follwed by calculating the PAM and B-mode images within suitably narrow bounds. The temporal window (L58) used
%    here suffices to capture wave front curvature from the target when received by a curvilinear array.
%
% 2. The axial step size in the PAM image domain (L51) is smaller than suggested by diffraction theory to give more precision
%    to the PAX calculation. An alternative turns out to be more time efficient is to use a coarser grid and interpolate the
%    image to a finer depth-grid.
%
% 3. The speed-depth survey is implemented here with four speed estimates (L47) spanning 1480 - 1570 m/s. This was adequate
%    but certainly using more, a non-uniformly distributed set of speeds, or an iterative approach may slightly improve the 
%    final PAX estimate.
%
% 4. The Bmode processing used here (L64, 87) is vastly simplified, but serves the general purpose of this demonstration. 
%    More rigorous reconstructions should be possible since the full receive time series data set is provided.

% This top level code calls the following additional functions:
%   teaImg_08Sep2022.m: creates a Time Exposure Acoustics image using conventional passive beamforming
%   delayfm.m:          applys explicit delays/leads to a matrix of time domain signal data
%   intersections.m:    finds the intersection of two curves in 2D space (written by Douglas M. Schwarz, dmschwarz@ieee.org)
%                       This function was obtained from the matlab file exchange:
%                       https://uk.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections?s_tid=ta_fx_results

% Michael Gray, University of Oxford, michael.gray@eng.ox.ac.uk, 08September2022

% NOTE: Make sure you update the datapath (line 41) and put this code and its support functions in the matlab path.
 

%% 1. Initial Setup

% Load the data set
datapath = 'C:\Users\engs1486\OneDrive - Nexus365\PAMC\Simulation\C52_3layer\';
disp('   Loading data set...')
load([datapath,'PAX_curvilinear_array_demo_data_08Sep2022.mat']);
Ne = length(aPos);                                          % number of array elements  
Nt = length(t);                                             % length of time vector

% Define a few sound speeds to try (same as in PAX IEEE paper)
cest = 1480:30:1570; Nc = length(cest);

% Define PAM imaging domain and presize storage
vy = (-0.25e-3 : 0.25e-3 : 0.25e-3);    Nvy = length(vy);   % lateral span, keep narrow in this example
vz = (-6e-3 : 0.10e-3 : 6e-3) + zs;     Nvz = length(vz);   % axial span

imT = zeros(Nvy,Nvz,Nc);                                    % image
imTmax = zeros(Nc,1);                                       % image peak value for each cest 
zmaxpam = zeros(Nc,1);                                      % axial location of image peak  

% Define 4cm temporal ROI based on near-center element response - use to accelerate PAM calculations
[~,ind] = max(sigs(Ne/2,:));                                % time index of signal peak
Nwin = ceil(0.04/1540*fs);                                  % window length
wsft = 80;                                                  % window start offset 
win = ind + (-(Nwin/2-wsft):(Nwin/2+wsft-1));               % window indices 

% A-line intensity from near-center array element  
xc = xcorr(sigs(Ne/2,:)',s0);                               % cross correlation of array signal with reference signal
[upper,lower] = envelope(xc);                               % envelope of correlation output
imB = upper(Nt:end)/max(upper(Nt:end));                     % backscatter intensity for positive lags
[~,Bind] = max(imB);                                        % index of maximum intensity
zmaxbmode = zeros(Nc,1);
Bz = zeros(Nt,Nc);



%% 2. Loop over speed estimates to calculate image peaks for PAM and Bmode
disp('   Running beamformers...')
for qc = 1:Nc
    
    % PAM
        % Form image
        imT(:,:,qc) = teaImg_08Sep2022(sigs(:,win) ,aPos, fs, cest(qc), vy, vz);
        
        % Find value and location of image maximum         
        [imTmax(qc), ind] = max(imT(:,:,qc),[],'all','linear');
        [~,v] = ind2sub([Nvy Nvz],ind);
        zmaxpam(qc) = vz(v)*1e3;                        % depth of peak, mm
    
        
    % Bmode
        % Form (approximate) A-line coordinate, and find peak
        Bz(:,qc) = t*cest(qc)/2*1e3;                    % depth vector (one way), mm
        zmaxbmode(qc) = Bz(Bind,qc);                    % depth of peak, mm
    
end


%% 3. Find intersection of Bmode and PAM (cest vs zmax) curves to estimate co-registration speed
disp('   Optimizing sound speed...')
[cpax, zpax] = intersections(cest, zmaxpam, cest, zmaxbmode);

% Recalculate PAM with optimized speed
imTpax = teaImg_08Sep2022(sigs(:,win), aPos, fs, cpax, vy, vz);
[imTmaxpax, ind] = max(imTpax,[],'all','linear');
[u,v] = ind2sub([Nvy Nvz],ind);
zmaxpampax = vz(v)*1e3;                        % depth of peak, mm

% Recalculate Bmode depth vector with optimized speed
Bzpax = t*cpax/2*1e3;  
zmaxbmodepax = Bzpax(Bind);                    % depth of peak, mm



%% 4. Plot it
figure
subplot(2,2,1)      % on-axis images using 1540 m/s or thereabouts
    cu = find(cest>=1540,1);
    plot(imT(vy==0,:,cu)/imTmax(cu),vz*1e3, 'r', imB, Bz(:,cu),'b') 
    ylim([vz(1)*1e3 vz(end)*1e3]), set(gca,'ydir','reverse')
    xlabel('Normalized Intensity'), ylabel('Depth (mm)')
    title(['On-axis Intensities, c_{est} = ',num2str(cest(cu))],'interpreter','tex')
    set(gca,'fontsize',14), legend({'PAM','B-mode'},'location','se')
    
subplot(2,2,2)      % on-axis images using cpax
    plot(imTpax(vy==0,:)/imTmaxpax,vz*1e3, 'r', imB, Bzpax,'b') 
    ylim([vz(1)*1e3 vz(end)*1e3]), set(gca,'ydir','reverse')
    xlabel('Normalized Intensity'), ylabel('Depth (mm)')
    title(['On-axis Intensities, c_{pax} = ',num2str(cpax,'%4.1f')],'interpreter','tex')
    set(gca,'fontsize',14)
    
subplot(2,2,3:4)    % speed/depth curves    
    plot(cest, zmaxpam,'-ro',cest, zmaxbmode,'-bs', cpax, zmaxpampax,'gd')
    set(gca,'ydir','reverse')
    ch = get(gca,'children'); for q=1:length(ch), set(ch(q),'MarkerSize',10); end   
    axis([(min(cest)-10) max(cest+10), 48, 60])
    xlabel('Sound Speed (m/s)'), ylabel('Image Peak Depth (mm)')
    title('Image Peaks')
    legend({'PAM','B-mode','PAX'},'location','s','orientation','horizontal')
    set(gca,'fontsize',14)
    
disp(['sound speeds: c_pax = ',num2str(cpax,'%4.1f'),', c_meanpath = ',num2str(cmean,'%4.1f'),' m/s'])
disp(['target depths: z_pax = ',num2str(zmaxpampax,'%2.1f'),', z_actual = ',num2str(zs*1e3,'%2.1f'),' mm'])


