function x = cg_read_curv_txt(name)
% Read ascii or freesurfer file
% 
%-Get a filename if none was passed
%-----------------------------------------------------------------------
if nargin==0
	name = spm_get(1,'any','Select ASCII data file');
end

[pth, nam, ext] = spm_fileparts(name);

% check for ascii, gifti or freesurfer format
if strcmp(ext,'.txt')

  fid = fopen(name, 'r');

  if fid == -1
    error(['Could not open file' name]);
  end

  x = fscanf(fid,'%f');

  fclose(fid);
elseif strcmp(ext,'.gii')
  V = spm_data_hdr_read(name);
  x = spm_data_read(V);
else
  x = read_curv(name);
end

return