function cg_write_curv_txt(name, x)
% Read ascii or freesurfer file
% 
%-Get a filename if none was passed
%-----------------------------------------------------------------------
if nargin > 2
  error('Wrong number of arguments');
end

[pth, basename, ext] = fileparts(name);

% check for ascii or freesurfer format
if strcmp(ext,'.txt')
  eval(['save -ascii ' name  ' x']);
  
else
  write_curv(name, x);
end

return

