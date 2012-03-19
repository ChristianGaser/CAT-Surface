function cg_cov_curv_txt
%function to calculate cc between surface measures

P = spm_select(Inf,'*','Select data');
n = size(P,1);

data = cg_read_curv_txt(deblank(P(1,:)));
m = length(data);
data_array = zeros(m,n);
data_array(:,1) = data;

for i=2:n
  data = cg_read_curv_txt(deblank(P(i,:)));
  if (length(data) ~= m)
    error('Data have to be same size and resampled.');
  end
  data_array(:,i) = data;
end

cc = corrcoef(data_array)

figure(11)
imagesc(cc);
set(gca,'XTick',1:n,'YTick',1:n)
colorbar
axis image