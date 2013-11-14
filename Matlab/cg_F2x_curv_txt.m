function cg_F2x_curv_txt
%CG_F2X transformation of F-maps to P, -log(P), r or d-maps
%
% The following formulas are used:
%
% --------------------------------
% coefficient of determination R2
% --------------------------------
%
%           F*(n-1)
% R2 = ------------------
%        n-p + F*(n-1)
%
% --------------------------------
% p-value
% --------------------------------
%
% p = 1-spm_Fcdf
%
% --------------------------------
% log p-value
% --------------------------------
%
% -log10(1-P) = -log(1-spm_Fcdf)
%
% For the last case of log transformation this means that a p-value of p=0.99 (0.01)
% is transformed to a value of 2
%
% Examples:
% p-value	-log10(1-P)
% 0.1		1
% 0.05		1.3
% 0.01		2
% 0.001		3
% 0.0001	4
%
%_______________________________________________________________________
% @(#)cg_f2x.m   1.22 Christian Gaser 2005/11/09


P = spm_select(Inf,'F.*','Select result image(s)');
sel = spm_input('Convert F value to?',1,'m',...
	'1-p|-log(1-p)|coefficient of determination R^2',[1 2 3], 2);

correction = spm_input('Thresholding?',1,'m',...
    'FDR|uncorrected p-value|no threshold',[1 2 3], 1);

if correction < 3
    u  = spm_input('p value (1 for no threshold)','+1','r',0.05,1,[0,1]);
end

for i=1:size(P,1)
	Res = deblank(P(i,:));
	[pth,nm,xt,vr] = fileparts(Res);

	SPM_name = fullfile(pth, ['SPM.mat' vr]);
	
	% SPM.mat exist?
	try
		load(SPM_name);
		spm_mat = 1;
	catch
		spm_mat = 0;		
	end

	if spm_mat
		xCon = SPM.xCon;
		if strcmp(nm(1:3),'F_0')	
			num = str2num(nm(4:5));
			df = [xCon(num).eidf SPM.xX.erdf];
		else
			df = spm_input('df ?',1,'r',[],2)';		
		end
	else
		df = spm_input('df ?',1,'r',[],2)';		
	end

  F = cg_read_curv_txt(Res);
  F(isnan(F)) = 0;
  PF = (1-spm_Fcdf(F,df));

  switch correction
  case 1
      Fs = flipud(sort(F(:)));
      Ps = (1-spm_Fcdf(Fs,df));

      S = length(Ps);
      Fi  = (1:S)'/S*u/1;
      I = max(find(Ps<=Fi));

      if isempty(I)
        F_threshold = Inf;
      else
        F_threshold = Fs(I);
      end
      PF(find(F < F_threshold)) = 1;
      F(find(F < F_threshold)) = 0;
      fprintf('Threshold F = %3.3f\n',F_threshold);
    case 2
      F_threshold = spm_invFcdf(1 - u, df);
      PF(find(F < F_threshold)) = 1;
      F(find(F < F_threshold)) = 0;
      fprintf('Threshold F = %3.3f\n',F_threshold);
  case 3
      F_threshold = 0;
      fprintf('\n');  
  end

	switch sel
	  case 1
	  	f2x = 1-spm_Fcdf(F,df);
		nm2 = 'P';
	  case 2
	  	f2x = -log10(max(eps,1-spm_Fcdf(F,df)));
		nm2 = 'logP';
	  case 3
	  	f2x = (df(2)-1)*F./(df(2) - df(1)+F*(df(2) -1));
		nm2 = 'R2';
	end

    switch correction
      case 1
        tmp_str=['FDR' num2str(u*100)];
      case 2
        tmp_str= [num2str(u*100)];
      case 3
        tmp_str='';
    end

	% name should follow convention spm?_0*
	out = fullfile(pth,['spm' nm2 '_' tmp_str nm xt vr]);

  if isfinite(F_threshold)  cg_write_curv_txt(out, f2x); end

end
