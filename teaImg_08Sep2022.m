function img = teaImg_08Sep2022(rf,arrayPos,fs,c,vy,vz)
% Generate passively beamformed map from rf data using 
% the Time Exposure Acoustics (TEA) algorithm (Gyongy et al 2008, DOI 10.1109/ULTSYM.2008.0210)
% M Gray, updated 08September 2022, based on original code by M Gyongy

%% Configure
ray.rho = 1000;                                             % assumed density
ray.c = c;                                                  % specified sound speed
ray.r = [arrayPos(:,1), arrayPos(:,2), arrayPos(:,3)];      % array element positions: [elevation, lateral, axial]
ray.N = size(arrayPos,1);                                   % number of array elements
ray.NRF = size(rf,2);                                       % number of time points
ray.fs = fs;                                                % data sample rate
t = (0:(ray.NRF-1))/fs;                                     % time vector


%% Define Imaging Domain
vox.y = vy;
vox.z = vz;

vox.x = 0;
vox.Nz = length(vox.z);
vox.Ny = length(vox.y);
vox.N = vox.Nz*vox.Ny;
vox.r = zeros(vox.N,3);

for Ivox = 1:vox.N
    Iz = ceil(Ivox/vox.Ny);
    Iy = mod(Ivox-1,vox.Ny)+1;
    vox.r(Ivox,:) = [vox.x vox.y(Iy) vox.z(Iz)]; % voxels go down columns
end

img = zeros(vox.Ny,vox.Nz);


%% Form Image
Escale = (4*pi*ray.NRF/fs)/(ray.rho*ray.c);                         % scaling term for energy calculation

for Ivox = 1:vox.N     

    dist = sqrt(sum(((ones(ray.N,1)*vox.r(Ivox,:))-ray.r).^2,2));   % slant range to voxel location
    dels = dist/ray.c;                                              % steering delays
    dels = dels - min(dels);                                        % delays with minimum value offset removed

    delrf = delayfm(t,rf,-dels).*repmat(dist(:),[1 ray.NRF]);       % Remove delays and scale by slant range
    dterm = mean(delrf,1);                                          % Average over elements
    
    % Final energy estimate: units are Joules if rf data are in Pascals
    img(Ivox) = Escale*( mean((dterm).^2) - mean(dterm).^2 );  
    
end

