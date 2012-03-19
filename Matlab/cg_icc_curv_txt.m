function cg_icc_curv_txt
%function to work out ICCs according to shrout & fleiss' schema (Shrout PE,
%Fleiss JL. Intraclass correlations: uses in assessing rater reliability.
%Psychol Bull. 1979;86:420-428) at each vertex of a surface mesh.
%
% 'Raters' is either 1,2,3. 'Raters' is: 1 if each target is measured by a
% different set of raters from a population of raters, 2 if each target is
% measured by the same raters, but that these raters are sampled from a
% population of raters, 3 if each target is measured by the same raters and
% these raters are the only raters of interest.
%
% 'Type of measurement' is either 'single' or 'k' & denotes whether the ICC is based on a
% single measurement or on an average of k measurements, where k = the
% number of ratings/raters.

n = Inf;
for i = 1:1000,
  P = spm_select([0 n],'*',['Select data of one hemisphere for rater/scanner ' num2str(i)]);
  n = size(P,1);
  if n < 2
    break;
  end;
  V{i} = P;
end

icc_typ = menu('Type of measurement','Single measure (1)','Average measure (k)');
icc_cse = menu('Raters','Different set of raters (1)','Same raters, random measures (2)','Same raters, mixed measures (3)');
icc_typ_str = str2mat('1','k');
icc_typ_str2 = str2mat('single','k');

pvalue  = spm_input('SB Prophecy: Level of reliability?',1,'r',0.9,1);

close all
[tmp, name]=spm_str_manip(spm_str_manip(str2mat(V{1},V{2}),'t'),'C');
if strcmp(name.s(end),'.')
  name.s = name.s(1:end-1);
end

fprintf('Model according to Shrout & Fleiss 1979: ICC(%d,%s)\n',icc_cse,icc_typ_str(icc_typ));

% number of subjects/targets
n = size(V{1},1);

% number of raters
m = i - 1;

data = cg_read_curv_txt(deblank(V{1}(1,:)));
data_array = zeros([n m length(data)]);

for i=1:m
  for j=1:n
    data_array(j,i,:) = cg_read_curv_txt(deblank(V{i}(j,:)));
  end
end

icc_array = zeros(length(data),1);
sb_prophecy_array = zeros(length(data),1);

for i=1:length(data)
  data = data_array(:,:,i);
  icc_array(i) = icc(icc_cse,deblank(icc_typ_str2(icc_typ,:)),data);
  sb_prophecy_array(i) = pvalue*(1-icc_array(i))/icc_array(i)*(1-pvalue);
end

name1  = spm_input('Name of output file?',1,'s',['ICC' num2str(icc_cse) icc_typ_str(icc_typ) '_' name.s]);
cg_write_curv_txt(name1, icc_array);
name2  = spm_input('name of output file?',1,'s',['SBprophecy_ICC' num2str(icc_cse) icc_typ_str(icc_typ) '_' name.s]);
cg_write_curv_txt(name, sb_prophecy_array);

function out = icc(cse,typ,dat)
%function to work out ICCs according to shrout & fleiss' schema (Shrout PE,
%Fleiss JL. Intraclass correlations: uses in assessing rater reliability.
%Psychol Bull. 1979;86:420-428).
%
% 'dat' is data whose rows represent different ratings or raters & whose
% columns represent different cases or targets being measured. Each target
% is assumed to be a random sample from a population of targets.
%
% 'cse' is either 1,2,3. 'cse' is: 1 if each target is measured by a
% different set of raters from a population of raters, 2 if each target is
% measured by the same raters, but that these raters are sampled from a
% population of raters, 3 if each target is measured by the same raters and
% these raters are the only raters of interest.
%
% 'typ' is either 'single' or 'k' & denotes whether the ICC is based on a
% single measurement or on an average of k measurements, where k = the
% number of ratings/raters.
%
% This has been tested using the example data in the paper by shrout & fleiss.
% 
% Example: out = ICC(3,'k',S_Fdata)
% returns ICC(3,k) of data 'S_Fdata' to double 'out'.
%
% Kevin Brownhill, Imaging Sciences, KCL, London kevin.brownhill@kcl.ac.uk
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%number of raters/ratings
k = size(dat,2);
%number of targets
n = size(dat,1);
%mean per target
mpt = mean(dat,2);
%mean per rater/rating
mpr = mean(dat);
%get total mean
tm = mean(mpt);
%within target sum sqrs
WSS = sum(sum(bsxfun(@minus,dat,mpt).^2));
%within target mean sqrs
WMS = WSS / (n * (k - 1));
%between rater sum sqrs
RSS = sum((mpr - tm).^2) * n;
%between rater mean sqrs
RMS = RSS / (k - 1);
% %get total sum sqrs
% TSS = sum(sum((dat - tm).^2));
%between target sum sqrs
BSS = sum((mpt - tm).^2) * k;
%between targets mean squares
BMS = BSS / (n - 1);
%residual sum of squares
ESS = WSS - RSS;
%residual mean sqrs
EMS = ESS / ((k - 1) * (n - 1));

switch cse
    case 1
        switch typ
            case 'single'
                out = (BMS - WMS) / (BMS + (k - 1) * WMS);
            case 'k'
                out = (BMS - WMS) / BMS;
            otherwise
               error('Wrong value for input typ') 
        end
    case 2
        switch typ
            case 'single'
                out = (BMS - EMS) / (BMS + (k - 1) * EMS + k * (RMS - EMS) / n);
            case 'k'
                out = (BMS - EMS) / (BMS + (RMS - EMS) / n);
            otherwise
               error('Wrong value for input typ') 
        end
    case 3
        switch typ
            case 'single'
                out = (BMS - EMS) / (BMS + (k - 1) * EMS);
            case 'k'
                out = (BMS - EMS) / BMS;
            otherwise
               error('Wrong value for input typ') 
        end
    otherwise
        error('Wrong value for input cse')
end