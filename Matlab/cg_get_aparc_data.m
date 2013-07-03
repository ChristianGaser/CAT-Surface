function cg_get_aparc_fsaverage

atlas_sel = 'aparc.a2009s|aparc';
atlas_names = char('aparc.a2009s','aparc');
meas_sel = 'thickness|area|amc_central|fd';
meas_names = char('thickness','area','amc_central','fd');
hemi_names = char('lh','rh');

P = spm_select(Inf,'dir','select subject folders');
n = size(P,1);

atlas = deblank(atlas_names(spm_input('Atlas',1,'m',atlas_sel),:));
meas  = deblank(meas_names(spm_input('Measure',2,'m',meas_sel),:));

for k=1:2
  hemi = hemi_names(k,:);

  output = [hemi '.' atlas '.' meas '.table'];
  
  for i=1:n
    name = deblank(P(i,:));
  
    [vertices,label,colortable] = read_annotation(fullfile(name,['label/' hemi '.' atlas '.annot']));
    curv = read_curv(fullfile(name,['surf/' hemi '.' meas]));

    label_val = colortable.table(:,5);
  
    if i==1
      annot_val = zeros(n,length(label_val)-1);
    end
  
    % skip first value for unknown areas
    for j=2:length(label_val)
      ind = find(label==label_val(j));
      annot_val(i,j-1) = mean(curv(ind));
    end
  
  end

  fid = fopen(output,'w');

  fprintf(fid,'Subject\t');
  for j=2:length(label_val)
    fprintf(fid,'%s\t',char(colortable.struct_names(j,:)));
  end

  fprintf(fid,'\n');

  P2 = spm_str_manip(P,'c');
  for i=1:n
    fprintf(fid,'%s\t',deblank(P2(i,:)));
    for j=2:length(label_val)
      fprintf(fid,'%3.3f\t',annot_val(i,j-1));  
    end
    fprintf(fid,'\n');
  end

  fclose(fid);
end
