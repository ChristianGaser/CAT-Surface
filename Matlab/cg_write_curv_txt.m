function cg_write_curv_txt(name, x)
% Read ascii or freesurfer file
% 
%-Get a filename if none was passed
%-----------------------------------------------------------------------
if nargin > 2
  error('Wrong number of arguments');
end

[pth, basename, ext] = spm_fileparts(name);

% check for ascii, gifti or freesurfer format
if strcmp(ext,'.txt')
  eval(['save -ascii ' name  ' x']);
  
elseif strcmp(ext,'.gii')
  save(gifti(struct('cdata',x)),name);
else
  write_curv(name, x);
end

return

