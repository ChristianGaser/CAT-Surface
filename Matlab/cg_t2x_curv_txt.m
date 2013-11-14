function r = cg_t2x_curv_txt
%CG_T2X transformation of t-maps to P, -log(P), r or d-maps
%
% The following formulas are used:
%
% --------------------------------
% correlation coefficient:
% --------------------------------
%             1
% r = ------------------
%           n-2
%     sqrt(------ + 1)
%          sqr(t)
%
% --------------------------------
% effect-size
% --------------------------------
%            2r
% d = ----------------
%     sqrt(1-sqr(r))
%
% --------------------------------
% p-value
% --------------------------------
%
% p = 1-spm_Tcdf
%
% --------------------------------
% log p-value
% --------------------------------
%
% -log10(1-P) = -log(1-spm_Tcdf)
%
% For the last case of log transformation this means that a p-value of p=0.99 (0.01)
% is transformed to a value of 2
%
% Examples:
% p-value   -log10(1-P)
% 0.1       1
% 0.05      1.3
% 0.01      2
% 0.001     3
% 0.0001    4
%
%_______________________________________________________________________
% @(#)cg_t2x.m   1.21 Christian Gaser 2006/04/27

n = 1;
STAT = 'T';

P = spm_select(Inf,'T.','Select result image(s)');

convert = spm_input('Convert t value to?',1,'m',...
    '1-p|-log(1-p)|correlation coefficient cc|effect size d|no conversion',[1 2 3 4 5], 2);

correction = spm_input('Thresholding?',1,'m',...
    'FDR|uncorrected p-value|no threshold',[1 2 3], 1);

if correction < 3
    u  = spm_input('p value (1 for no threshold)','+1','r',0.05,1,[0,1]);
end

neg_results = spm_input('Show also inverse effects (e.g. neg. values)','+1','b','yes|no',[1 0],2);

for i=1:size(P,1)

    T_threshold = 0;
    Res = deblank(P(i,:));
    [pth,nm,xt,vr] = spm_fileparts(Res);
    fprintf('%s: ',nm);

    SPM_name = fullfile(pth, ['SPM.mat' vr]);
    
    % SPM.mat exist?
    try
        load(SPM_name);
        spm_mat = 1;
    catch
        spm_mat = 0;        
    end

    if spm_mat
        df = SPM.xX.erdf;
    else
        df = spm_input('df ?',1,'e');       
    end

    T = cg_read_curv_txt(Res);
    T(isnan(T)) = 0;
    PZ = (1-spm_Tcdf(T,df)).^n;

    switch correction
    case 1

        Ts = flipud(sort(T(:)));
        Ps = (1-spm_Tcdf(Ts,df)).^n;

        S = length(Ps);
        Fi  = (1:S)'/S*u/1;
        I = max(find(Ps<=Fi));

        if isempty(I)
          T_threshold = Inf;
        else
          T_threshold = Ts(I);
        end

        if neg_results
          ind_T = find((T < T_threshold) & (T > -T_threshold));
        else
          ind_T = find(T < T_threshold);
        end

        PZ(ind_T) = 1;
        T(ind_T) = 0;
        fprintf('Threshold T = %3.3f\n',T_threshold);
    case 2
        T_threshold = spm_invTcdf(1 - u, df);
        if neg_results
          ind_T = find((T < T_threshold) & (T > -T_threshold));
        else
          ind_T = find(T < T_threshold);
        end
        PZ(ind_T) = 1;
        T(ind_T) = 0;
        fprintf('Threshold T = %3.3f\n',T_threshold);
    case 3
        fprintf('\n');  
    end
    
    switch convert
      case 1
        t2x = PZ;
        nm2 = 'P';
      case 2
        t2x = -log10(max(eps,PZ));
        ind = find(T < 0);
        if ~isempty(ind)
            t2x(ind) = log10(max(eps,1-PZ(ind)));
        end
        ind = find(T == 0);
        if ~isempty(ind)
            t2x(ind) = 0;
        end
        nm2 = 'logP';
      case 3
        t2x = sign(T).*(1./((df./((T.*T)+eps))+1)).^0.5;
        nm2 = 'R';
      case 4
        tmp = ((df./((T.*T)+eps))+1);
        t2x = 2./((1-(1./tmp)).*tmp).^0.5;
        nm2 = 'D';
      case 5
        t2x = T;
        nm2 = 'T';
    end

    switch correction
      case 1
        tmp_str=['FDR' num2str(u*100)];
      case 2
        tmp_str= [num2str(u*100)];
      case 3
        tmp_str='';
    end
    
    if neg_results
      neg_str = '_bi'; 
    else
      neg_str = '';
    end

    if strcmp(nm(1),'T')    
        out = fullfile(pth,[nm2 '_' tmp_str nm(2:end) neg_str xt vr]);
    else
        out = fullfile(pth,[nm2 '_'  tmp_str nm neg_str xt vr]);
    end
    
    if isfinite(T_threshold) cg_write_curv_txt(out, t2x); end

end
