function x = cg_read_curv_txt(name)
% Read ascii or freesurfer file
% 
%-Get a filename if none was passed
%-----------------------------------------------------------------------
if nargin==0
	name = spm_get(1,'*','Select ASCII data file');
end

[pth, nam, ext] = fileparts(name);

% check for ascii or freesurfer format
if strcmp(ext,'.txt')

  fid = fopen(name, 'r');

  if fid == -1
    error(['Could not open file' name]);
  end

  x = fscanf(fid,'%f');

  fclose(fid);
else
  x = read_curv(name);
end
return