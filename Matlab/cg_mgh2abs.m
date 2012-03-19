function cg_mgh2abs

P = spm_select(Inf,'.mgh','Select mgh-files');

n = size(P,1);
for i=1:n
  inname = deblank(P(i,:));
  [pth, name, ext] = spm_fileparts(inname);
  [vol, M] = load_mgh(inname);
  outname = fullfile(pth,[name '_abs' ext]);
  save_mgh(abs(vol),outname,M);
end