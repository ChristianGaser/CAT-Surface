function r = cg_invert_values_curv_txt
% invert scalar values
%
%_______________________________________________________________________
% @(#)cg_invert_values_curv_txt.m   1.00 Christian Gaser 2010/01/06

n = 1;

P = spm_select(Inf,'.*','Select inult image(s)');

for i=1:size(P,1)

    in = deblank(P(i,:));
    [pth,nm,xt,vr] = fileparts(in);
    fprintf('%s: \n',nm);

    values = cg_read_curv_txt(in);
    values = -values;
    
    out = fullfile(pth,[nm '_inv' xt]);
    cg_write_curv_txt(out, values);

end
