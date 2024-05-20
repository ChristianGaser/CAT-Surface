function make_sphere_phantom

simu_thickness = 3;
csf_val = 1;
gm_val  = 2;
wm_val  = 3;

name = sprintf('sphere_dilated%gmm.nii',simu_thickness);

vx = [0.5 0.5 0.5];
N = nifti;
N.dat = file_array(name,[320 320 320],[spm_type('float32'); 0],0,1,0);;
N.mat = eye(4);
for i=1:3, N.mat(i,i) = vx(i); end
N.mat0 = N.mat;
N.descrip  = 'Phantom';
N.mat_intent = 'Scanner';
N.mat0_intent = 'Scanner';
create(N);

fwhm   = [30];
coord  = round(N.dat.dim/2);

vol = zeros(N.dat.dim(1:3));

for j=1:size(coord,1)
  % compute parameters for spm_conv_vol
  %-----------------------------------------------------------------------
  s  = fwhm(j,:)./vx;       % voxel anisotropy
  s  = max(s,ones(size(s)));      % lower bound on FWHM
  s  = s/sqrt(8*log(2));        % FWHM -> Gaussian parameter    
  x  = round(6*s(1)); x = [-x:x];
  y  = round(6*s(2)); y = [-y:y];
  z  = round(6*s(3)); z = [-z:z];
  x  = exp(-(x).^2/(2*(s(1)).^2)); 
  y  = exp(-(y).^2/(2*(s(2)).^2)); 
  z  = exp(-(z).^2/(2*(s(3)).^2));

  xi  = (length(x) - 1)/2;
  xj  = (length(y) - 1)/2;
  xk  = (length(z) - 1)/2;

  % 3D gaussian
  xy = x'*y;
  xyz  = zeros([size(xy) length(z)]);
  for i=1:length(z)
    xyz(:,:,i) = xy'*z(i);
  end

  % range of image
  zx = coord(j,1)-xi:coord(j,1)+xi;
  zy = coord(j,2)-xj:coord(j,2)+xj;
  zz = coord(j,3)-xk:coord(j,3)+xk;
  
  vol(zx,zy,zz) = xyz;

end

% create wm and csf
ind_csf = vol>0.00001;
ind_gm  = vol>0.0001;
ind_wm  = vol>0.001;
vol(ind_csf) = csf_val;
vol(ind_gm)  = gm_val;
vol(ind_wm)  = wm_val;
ind = vol < 1E-5;

vol(1:75,:,:) = gm_val;
vol(245:end,:,:) = gm_val;
vol(:,1:75,:) = gm_val;
vol(:,245:end,:) = gm_val;
vol(:,:,1:75) = gm_val;
vol(:,:,245:end) = gm_val;
vol(ind) = 0;

% euclidean distance to wm
wm  = vol==wm_val;
Dwm = (bwdist(wm) - 0.5) * mean(vx);

% set gm and csf to csf
vol(vol<3 & vol>1) = csf_val;

% set dilated gm to gm
vol(vol < wm_val & Dwm <= simu_thickness) = gm_val;

% check distance inside gm
D = (bwdist(vol==csf_val) + bwdist(vol==wm_val) - 1)*mean(vx);

D(vol ~= gm_val) = 0;

figure
hist(D(find(D>0)),100); %xlim([1 3])

N.dat(:,:,:) = vol;

