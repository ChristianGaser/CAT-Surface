% a powerful way of creating average values for each twin pair
function create_mean_image_twinpair

for i=1:1000,
    P = spm_select([0 2],'*',['Select data of one hemisphere for one twin pair ' num2str(i)]);
    n = size(P,1);
    if n < 2
      break;
    end;
    
    V{i} = P;
end

for i = 1:length(V)
    [pth,nm,xt,vr] = fileparts(V{i}(1,:));
    name = fullfile(pth,['avg_' nm xt]);
    avg = cg_read_curv_txt(deblank(V{i}(1,:))) + cg_read_curv_txt(deblank(V{i}(2,:)));
    avg = avg/2;
    cg_write_curv_txt(name,avg);
end

return

