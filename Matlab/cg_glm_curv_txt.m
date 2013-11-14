function [Ic,xCon] = cg_glm_curv_txt(Ic,xCon)
%FORMAT [Ic,xCon] = cg_glm_curv_txt(Ic,xCon)
%
% Design structure fields and option definitions
% ==============================================
%
% D.Desname  - a string naming the design
%
% In general, cg_glm_curv_txt.m accomodates four factors. Usually these are
% 'group', 'subject', 'condition' & 'replication', but to allow for a
% flexible interface these are dynamically named for different designs,
% and are referred to as Factor4, Factor3, Factor2, and Factor1
% respectively. The first part of the D definition dictates the names
% and number of factor levels (i.e. number of subjects etc.) relevant
% for this design, and also how the H (condition) and B (block)
% partitions of the design matrix should be constructed.
%
% D.n        - a 1x4 vector, indicating the number of levels. D.n(i)
%              for i in [1:4] is the number of levels for factor i.
%              Specify D.n(i) as 1 to ignore this factor level,
%              otherwise the number of levels can be pre-specified as a
%              given number, or given as Inf to allow the user to
%              choose the number of levels.
%
% D.sF       - a 1x4 cellstr containing the names of the four
%              factors. D.sF{i} is the name of factor i.
%
% D.b.aTime  - a binary indicator specifying whether images within F3
%              level (subject) are selected in time order. For time
%              order (D.b.aTime=1), F2 levels are indicated by a user
%              input "condition" string (input handled by spm_input's
%              'c' type). When (D.b.aTime=0), images for each F3 are
%              selected by F2 (condition). The latter was the mode of
%              SPM95 and SPM96. (SPM94 and SPMclassic didn't do
%              replications of conditions.)
%
% Once the user has entered the images and indicated the factor levels,
% a nScan x 4 matrix, I, of indicator variables is constructed
% specifying for each scan the relevant level of each of the four
% factors. I(n,i) is the level of factor i corresponding to image n.
% This I matrix of factor indicators is then used to construct the H
% and B forms of the design matrix according to the prescripton in the
% design definition D:
%
% D.Hform    - a string specifying the form of the H partition of the
%              design matrix. The string is evaluated as an argument
%              string for spm_DesMtx, which builds design matrix
%              partitions from indicator vectors.
%             (eval(['[H,Hnames] = spm_DesMtx(',D.Hform,');']))
%
% D.BForm   - a string specifying the form of the G partition.
%
% ( Note that a constant H partition is dropped if the B partition can   )
% ( model the constant effect.                                           )
%
% The next part of the design definition defines covariate options.
% Covariates are split into covariates (of interest) and nuisance
% variables. The covariates of interest and nuisance variables are put
% into the C & G partitions of the design matrox (the final design
% matrix is [H,C,B,G], where global nuisance covariates are appended to
% G). In SPM94/5/6 the design matrix was partitioned into effects of
% interest [H,C] and effects of no interest [B,G], with an F-test for
% no effects of interest and adjusted data (for effects of no interest)
% following from these partitions. SPM99 is more freestyle, with
% adjustments and F-tests specified by contrasts. However, the concept
% of effects of interest and of no interest has been maintained for
% continuity, and cg_glm_curv_txt.m computes an F-contrast to test for "no
% effects of interest".
%
% D.nC       - a 1x2 vector: D.nC(1) is the number of covariates,
%              D.nC(2) the number of nuisance variables. Specify zero
%              to skip covariate entry, the actual number of
%              covariates, or Inf to let the user specify the number of
%              covariates. As with earlier versions, blocks of design
%              matrix can be entered. However, these are now treated as
%              a single covariate entity, so the number of
%              covariates.nuisance variables is now the number of items
%              you are prompted for, regardless of their dimension. (In
%              SPM95-6 this number was the number of covariate vectors
%              that could be entered.)
%
% D.iCC      - a 1x2 cell array containing two vectors indicating the
%              allowable covariate centering options for this design.
%              These options are defined in the body of cg_glm_curv_txt.m,
%              in variables sCC & CFIforms. Use negative indices to
%              indicate the default, if any - the largest negative
%              wins.
%
% D.iCFI     - a 1x2 cell array containing two vectors indicating the
%              allowable covariate by factor interactions for this
%              design. Interactions are only offered with a factor if
%              it has multiple levels. The options are defined in the
%              body of cg_glm_curv_txt.m, in variables sCFI & CFIforms. Use
%              negative indicies to indicate a default.
%
% The next part defines global options:
%
% D.iGXcalc  - a vector of possible global calculation options for
%              this design, as listed in the body of cg_glm_curv_txt.m in
%              variable sGXcalc. (If other global options are chosen,
%              then the "omit" option is not offered.) Again, negative
%              values indicate a default.
%
% D.iGloNorm - a vector of possible global normalisation options for
%              this design, as described in the body of cg_glm_curv_txt.m in
%              variable sGloNorm.
%
% D.iGMsca   - a vector of possible grand mean scaling options, as
%              described in the body of cg_glm_curv_txt.m in variable
%              sGMsca. (Note that grand mean scaling is redundent when
%              using proportional scaling global flow normalisation.)
%
% D.iGC      - a vector of possible global covariate centering
%              options, corresponding to the descriptions in variable
%              iCC given in the body of cg_glm_curv_txt.m. This is only
%              relevant for AnCova type global normalisation, and even
%              then only if you're actually interested in constraining
%              the values of the parameters in some useful way.
%              Usually, one chooses option 10, "as implied by AnCova".
%
% The next component specifies masking options:
%
% D.M_.T     - a vector defining the analysis threshold: Specify
%              "-Inf" as an element to offer "None" as an option. If a
%              real element is found, then absolute thresholding is
%              offered, with the first real value proffered as default
%              threshold. If an imaginary element is found, then
%              proportional thresholding if offered (i.e. the threshold
%              is a proportion of the image global), with the (abs of)
%              the first imaginary element proffered as default.
%
% D.M_.I     - Implicit masking? 0-no, 1-yes, Inf-ask. (This is
%              irrelevant for image types with a representation of NaN,
%              since NaN is then the mask value, and NaN's are always
%              masked.)
%
% D.M.X      - Explicit masking? 0-no, 1-yes, Inf-ask.
% 
%                           ----------------
%
% To use a customised design structure D, type cg_glm_curv_txt('cfg',D) in the
% Matlab command window.
%
% The easiest way to generate a customised design definition structure
% is to tweak one of the pre-defined definitions. The following code
% will prompt you to select one of the pre-defined designs, and return
% the design definition structure for you to work on:
%


%-Option definitions
%-----------------------------------------------------------------------
%-Generic factor names
sF = {'sF1','sF2','sF3','sF4'};

%-Covariate by factor interaction options
sCFI = {'<none>';...                                        %-1
      'with sF1';'with sF2';'with sF3';'with sF4';...             %-2:5
      'with sF2 (within sF4)';'with sF3 (within sF4)'};           %-6,7

%-DesMtx argument components for covariate by factor interaction options
% (Used for CFI's Covariate Centering (CC), GMscale & Global normalisation)
CFIforms = {      '[]',       'C',  '{}';...                %-1
            'I(:,1)',   'FxC',      '{D.sF{1}}';...               %-2
            'I(:,2)',   'FxC',      '{D.sF{2}}';...               %-3
            'I(:,3)',   'FxC',      '{D.sF{3}}';...               %-4
            'I(:,4)',   'FxC',      '{D.sF{4}}';...               %-5
            'I(:,[4,2])',     'FxC',      '{D.sF{4},D.sF{2}}';...       %-6
            'I(:,[4,3])',     'FxC',      '{D.sF{4},D.sF{3}}'     };    %-7

%-Centre (mean correction) options for covariates & globals            (CC)
% (options 9-12 are for centering of global when using AnCova GloNorm) (GC)
sCC = {           'around overall mean';...                       %-1
            'around sF1 means';...                          %-2
            'around sF2 means';...                          %-3
            'around sF3 means';...                          %-4
            'around sF4 means';...                          %-5
            'around sF2 (within sF4) means';...             %-6
            'around sF3 (within sF4) means';...             %-7
            '<no centering>';...                            %-8
            'around user specified value';...               %-9
            '(as implied by AnCova)';...                    %-10
            'GM';...                                  %-11
            '(redundant: not doing AnCova)'}';              %-12
%-DesMtx I forms for covariate centering options
CCforms = {'ones(nScan,1)',CFIforms{2:end,1},''}';


%-Global normalization options (options 1-7 match CFIforms)       (GloNorm)
sGloNorm = {      'AnCova';...                                    %-1
            'AnCova by sF1';...                             %-2
            'AnCova by sF2';...                             %-3
            'AnCova by sF3';...                             %-4
            'AnCova by sF4';...                             %-5
            'AnCova by sF2 (within sF4)';...                %-6
            'AnCova by sF3 (within sF4)';...                %-7
            'proportional scaling';...                      %-8
            '<no global normalisation>'};                   %-9

%-Grand mean scaling options                                        (GMsca)
sGMsca = {  'scaling of overall grand mean';...             %-1
            'scaling of sF1 grand means';...                %-2
            'scaling of sF2 grand means';...                %-3
            'scaling of sF3 grand means';...                %-4
            'scaling of sF4 grand means';...                %-5
            'scaling of sF2 (within sF4) grand means';...         %-6
            'scaling of sF3 (within sF4) grand means';...         %-7
            '(implicit in PropSca global normalisation)';...      %-8
            '<no grand Mean scaling>'     };                %-9
%-NB: Grand mean scaling by subject is redundent for proportional scaling


%-Global calculation options                                       (GXcalc)
sGXcalc  = {      'omit';...                                %-1
            'user specified';...                            %-2
            'mean voxel value (within per image fullmean/8 mask)'};     %-3



%=======================================================================
%-D E S I G N   P A R A M E T E R S
%=======================================================================
%-Get design type
%-----------------------------------------------------------------------

sel = spm_input('Select design class...',1,'m',...
            'Cross-sectional|Cross-sectional different style|Longitudinal',...
            1:3);

switch sel
case 1

    D = struct(...
            'DesName','Cross-sectional design (1 scan per subject)',...
            'n',  [Inf Inf 1 1],    'sF',{{'repl','group','',''}},...
            'Hform',          'I(:,2),''-'',''group''',...
            'Bform',          'I(:,3),''-'',''\mu''',...
            'nC',[Inf,Inf],'iCC',{{[1,3,4,8],[1,4,8]}},'iCFI',{{[1,3,-4],[1,-4]}},...
            'iGXcalc',1,'iGMsca',9,'GM',[],...
            'iGloNorm',9,'iGC',12,...
            'M_',struct('T',0.1,'I',0,'X',Inf),...
            'b',struct('aTime',0));
            

case 2

    D = struct(...
            'DesName','Cross-sectional design (1 scan per subject)',...
            'n',  [Inf Inf 1 1],    'sF',{{'repl','group','',''}},...
            'Hform',          'I(:,2),''-'',''group''',...
            'Bform',          'I(:,3),''-'',''\mu''',...
            'nC',[Inf,Inf],'iCC',{{[1,3,4,8],[1,4,8]}},'iCFI',{{[1,3,-4],[1,-4]}},...
            'iGXcalc',1,'iGMsca',9,'GM',[],...
            'iGloNorm',9,'iGC',12,...
            'M_',struct('T',0.1,'I',0,'X',Inf),...
            'b',struct('aTime',1));

            
case 3
            
    D = struct(...
            'DesName','Longitudinal design: (more than 1 time point per subject)',...
            'n',[Inf Inf Inf Inf],  'sF',{{'repl','cond','subject','group'}},...
            'Hform',          'I(:,[4,2]),''-'',{''group'',''repl''}',...
            'Bform',          'I(:,[4,3]),''-'',{''group'',''subj''}',...
            'nC',[Inf,Inf],'iCC',{{[5:8],[5,7,8]}},'iCFI',{{[1,5,6,-7],[1,5,-7]}},...
            'iGXcalc',1,'iGMsca',9,'GM',[],...
            'iGloNorm',9,'iGC',12,...
            'M_',struct('T',0.1,'I',0,'X',Inf),...
            'b',struct('aTime',1));

end

%-Set factor names for this design
%-----------------------------------------------------------------------
sCC      = sf_estrrep(sCC,[sF',D.sF']);
sCFI     = sf_estrrep(sCFI,[sF',D.sF']);
sGloNorm = sf_estrrep(sGloNorm,[sF',D.sF']);
sGMsca   = sf_estrrep(sGMsca,[sF',D.sF']);

%-Get filenames & factor indicies
%-----------------------------------------------------------------------
[P,I]    = cg_get_files(D.sF,D.n,D.b.aTime,'any');
%[P,I]    = spm_spm_ui('Files&Indices',D.sF,D.n,D.b.aTime);
nScan    = size(I,1);                                 %-#obs

%-Additional design parameters
%-----------------------------------------------------------------------
bL       = any(diff(I,1),1);  %-Multiple factor levels?
          % NB: bL(2) might be thrown by user specified f1 levels
          %     (D.b.aTime & D.n(2)>1) - assumme user is consistent?
bFI      = [bL(1),bL(2:3)&~bL(4),bL(4),bL([2,3])&bL(4)];
          %-Allowable interactions for covariates
          %-Only offer interactions with multi-level factors, and
          % don't offer by F2|F3 if bL(4)!

%-Build Condition (H) and Block (B) partitions
%=======================================================================
eval(['[H,Hnames] = spm_DesMtx(',D.Hform,');'])
if rank(H)==nScan, error('unestimable condition effects'), end
eval(['[B,Bnames] = spm_DesMtx(',D.Bform,');'])
if rank(B)==nScan, error('unestimable block effects'), end

%-Drop a constant H partition if B partition can model constant
if size(H,2)>0 & all(H(:)==1) & (rank([H B])==rank(B))
      H = []; Hnames = {};
      warning('Dropping redundant constant H partition')
end


%-Covariate partition(s): interest (C) & nuisance (G) excluding global
%=======================================================================
nC = D.nC;              %-Default #covariates
C  = {[],[]}; Cnames = {{},{}};     %-Covariate DesMtx partitions & names
xC = [];                %-Struct array to hold raw covariates


dcname = {'CovInt','NusCov'}; %-Default root names for covariates
dstr   = {'covariate','nuisance variable'};

GUIpos = spm_input('!NextPos');
nc     = [0,0];
for i  = 1:2                  % 1:covariates of interest, 2:nuisance variables

    if isinf(nC(i)), nC(i)=spm_input(['# ',dstr{i},'s'],GUIpos,'w1'); end

    while nc(i) < nC(i)

      %-Create prompt, get covariate, get covariate name
        %---------------------------------------------------------------
      if nC(i)==1, str=dstr{i}; else, str=sprintf('%s %d',dstr{i},nc(i)+1); end
        c = spm_input(str,GUIpos,'r',[],[nScan,Inf]);
        if any(isnan(c(:))), break, end         %-NaN is dummy value to exit
      nc(i)  = nc(i)+1;             %-#Covariates (so far)
      if nC(i)>1, tstr = sprintf('%s^{%d}',dcname{i},nc(i));
      else,       tstr = dcname{i}; end
            cname  = spm_input([str,' name?'],'+1','s',tstr);
            rc     = c;                   %-Save covariate value
      rcname = cname;                     %-Save covariate name

        %-Interaction option? (if single covariate vector entered)?
        %---------------------------------------------------------------
        if size(c,2) == 1
                if length(D.iCFI{i})>1
                  %-User choice of interaction options, default is negative
                  %-Only offer interactions for appropriate factor combinations
            iCFI = intersect(abs(D.iCFI{i}),find([1,bFI]));
            dCFI = max([1,intersect(iCFI,-D.iCFI{i}(D.iCFI{i}<0))]);
            iCFI = spm_input([str,': interaction?'],'+1','m',...
                  sCFI(iCFI),iCFI,find(iCFI==dCFI));
          else
            iCFI = abs(D.iCFI{i});        %-AutoSelect default option
          end
      else
          iCFI = 1;
      end

        %-Centre covariate(s)? (Default centring to correspond to CFI)
        % Always offer "no centering" as default for design matrix blocks
        %---------------------------------------------------------------
      DiCC = D.iCC{i};
      if size(c,2)>1, DiCC = union(DiCC,-8); end
        if length(DiCC)>1
            %-User has a choice of centering options
            %-Only offer factor specific for appropriate factor combinations
            iCC = intersect(abs(DiCC),find([1,bFI,1]) );
            %-Default is max -ve option in D, overridden by iCFI if CFI
            if iCFI == 1, dCC = -DiCC(DiCC<0); else, dCC = iCFI; end
            dCC = max([1,intersect(iCC,dCC)]);
            iCC = spm_input([str,': centre?'],'+1','m',...
                  sCC(iCC),iCC,find(iCC==dCC));
        else
            iCC = abs(DiCC);  %-AutoSelect default option
        end
      %-Centre within factor levels as appropriate
        if any(iCC == [1:7]), c = c - spm_meanby(c,eval(CCforms{iCC})); end

        %-Do any interaction (only for single covariate vectors)
        %---------------------------------------------------------------
        if iCFI > 1                       %-(NB:iCFI=1 if size(c,2)>1)
                  tI        = [eval(CFIforms{iCFI,1}),c];
            tConst    = CFIforms{iCFI,2};
            tFnames   = [eval(CFIforms{iCFI,3}),{cname}];
            [c,cname] = spm_DesMtx(tI,tConst,tFnames);
      elseif size(c,2)>1                  %-Design matrix block
            [null,cname] = spm_DesMtx(c,'X',cname);
      else
            cname = {cname};
      end

      %-Store raw covariate details in xC struct for reference
      %-Pack c into appropriate DesMtx partition
        %---------------------------------------------------------------
      %-Construct description string for covariate
      str = {sprintf('%s: %s',str,rcname)};
      if size(rc,2)>1, str = {sprintf('%s (block of %d covariates)',...
            str{:},size(rc,2))}; end
      if iCC < 8, str=[str;{['used centered ',sCC{iCC}]}]; end
      if iCFI> 1, str=[str;{['fitted as interaction ',sCFI{iCFI}]}]; end

      tmp       = struct(     'rc',rc,    'rcname',rcname,...
                        'c',c,            'cname',{cname},...
                        'iCC',iCC,  'iCFI',iCFI,...
                        'type',i,...
                        'cols',[1:size(c,2)] + ...
                                    size([H,C{1}],2) +  ...
                                    size([B,C{2}],2)*(i-1),...
                        'descrip',{str}                     );
      if isempty(xC), xC = tmp; else, xC = [xC,tmp]; end
      C{i}      = [C{i},c];
      Cnames{i} = [Cnames{i}; cname];

    end     % (while)

end % (for)
clear c tI tConst tFnames
spm_input('!SetNextPos',GUIpos);

%-Unpack into C & G design matrix sub-partitions
G = C{2}; Gnames = Cnames{2};
C = C{1}; Cnames = Cnames{1};

%-Construct full design matrix (X), parameter names and structure (xX)
%=======================================================================
X      = [H C B G];
tmp    = cumsum([size(H,2), size(C,2), size(B,2), size(G,2)]);
xX     = struct(  'X',        X,...
                  'iH',       [1:size(H,2)],...
                  'iC',       [1:size(C,2)] + tmp(1),...
                  'iB',       [1:size(B,2)] + tmp(2),...
                  'iG',       [1:size(G,2)] + tmp(3),...
                  'name',           {[Hnames; Cnames; Bnames; Gnames]},...
                  'I',        I,...
                  'sF',       {D.sF});
                  
[nScan nBeta] = size(xX.X);

%-Design space and projector matrix [pseudoinverse] for WLS
%=======================================================================
xX.xKXs = spm_sp('Set',xX.X);                   % KWX
xX.pKX  = spm_sp('x-',xX.xKXs);                       % projector
xX.V = eye(nScan);
xX.pKXV  = xX.pKX*xX.V;                   %-for contrast variance weight
xX.Bcov  = xX.pKXV*xX.pKX';               %-Variance of est. param.
[xX.trRV,xX.trRVRV] ...                   %-Variance expectations
         = spm_SpUtil('trRV',xX.xKXs,xX.V);
xX.erdf  = xX.trRV^2/xX.trRVRV;                 %-Effective residual d.f.

if nargin == 0
  Fcname = 'effects of interest';
  iX0    = [xX.iB xX.iG];
  xCon   = spm_FcUtil('Set',Fcname,'F','iX0',iX0,xX.xKXs);
end

SPM.xX = xX;
SPM.xCon = xCon;
SPM.xY.VY = cell2struct(P,'fname',2);
tmp = {	sprintf('%d condition, +%d covariate, +%d block, +%d nuisance',...
		size(H,2),size(C,2),size(B,2),size(G,2));...
	sprintf('%d total, having %d degrees of freedom',...
		size(X,2),rank(X));...
	sprintf('leaving %d degrees of freedom from %d images',...
		size(X,1)-rank(X),size(X,1))				};
xsDes = struct(	'Design',			{D.DesName},...
		'Global_calculation',		{sGXcalc{1}},...
		'Grand_mean_scaling',		{sGMsca{9}},...
		'Global_normalisation',		{sGloNorm{9}},...
		'Parameters',			{tmp}			);
SPM.xsDes = xsDes;

if nargin == 0
  [Ic,xCon] = spm_conman(SPM,'T|F',1,...
            '     Select contrasts...',' for conjunction',1);
end

spm_progress_bar('Init',nScan,'reading data...','subjects completed');

Y = cg_read_curv_txt(deblank(char(P(1,:))));
[pth, name, ext] = fileparts(deblank(char(P(1,:))));

spm_progress_bar('Set',1)
for i=2:nScan
      tmp = cg_read_curv_txt(deblank(char(P(i,:))));
      if (size(tmp,1)==size(Y(:,1),1) & size(tmp,2)==size(Y(:,1),2))
            Y = [Y tmp];
      else
            error('Size of data differs');
      end
      spm_progress_bar('Set',i)
end
spm_progress_bar('Clear');

KY = Y';

beta = xX.pKX * KY;                 % parameter estimates
res   = spm_sp('r',xX.xKXs,KY);           %-Residuals
ResSS = sum(res.^2)';               %-Residual SSQ
ResMS = ResSS/xX.trRV;

% ResMS
Q = ['ResMS' ext];
cg_write_curv_txt(Q, ResMS)

% betas
for i=1:nBeta
      dat = beta(i,:)';
      Q = sprintf(['beta0%d' ext],i);
      cg_write_curv_txt(Q, dat)
end

% T-tests
for j=1:size(xCon,2)
      if xCon(j).STAT == 'T'
        % T-test
        con = (xCon(j).c'*beta)';     
        t   = con./(sqrt(ResMS*(xCon(j).c'*xX.Bcov*xCon(j).c)));

        name     = replace_char(xCon(j).name);
        Q = ['con_' name ext];
        cg_write_curv_txt(Q, con)
        Q = ['T_' name ext];
        cg_write_curv_txt(Q, t)
      elseif xCon(j).STAT == 'F'
        % F-test
        X1o           = spm_FcUtil('X1o',xCon(j),xX.xKXs);
        [trMV,trMVMV] = spm_SpUtil('trMV',X1o,xX.V);
        SPM.xCon(j).eidf = trMV^2/trMVMV;
        h       = spm_FcUtil('Hsqr',xCon(j),xX.xKXs);
        e = h*beta;
        ss = sum(e.^2,1);
        MVM = ss/trMV;
        RVR = ResMS;
        Z = MVM'./RVR;
        name     = replace_char(xCon(j).name);
        Q = sprintf(['F_%03d_%s' ext],j,name);
        cg_write_curv_txt(Q, Z)
        fprintf('Degrees of freedom for %s: %g,%g\n',name,[SPM.xCon(j).eidf xX.erdf]);
      else
        error('Only F and T-tests are supported.');
      end
end

save SPM SPM

return

%=======================================================================

function str = sf_estrrep(str,srstr)
%=======================================================================
for i = 1:size(srstr,1)
      str = strrep(str,srstr{i,1},srstr{i,2});
end

return

function varargout = cg_get_files(varargin)
% Setting up the general linear model for independent data
%_______________________________________________________________________
% based on @(#)spm_spm_ui.m   2.49 Andrew Holmes 03/03/20
SCCSid  = '2.49';

%=======================================================================
% - FORMAT specifications for programers
%=======================================================================
%( This is a multi function function, the first argument is an action  )
%( string, specifying the particular action function to take.          )
%
%
% FORMAT [P,I] = cg_get_files(DsF,Dn,DbaTime)
% PET/SPECT file & factor level input
% DsF     - 1x4 cellstr of factor names (ie D.sF)
% Dn      - 1x4 vector indicating the number of levels (ie D.n)
% DbaTime - ask for F3 images in time order, with F2 levels input by user?
% P       - nScan x 1 cellsrt of image filenames
% I       - nScan x 4 matrix of factor level indices
%


%-Condition arguments
%-----------------------------------------------------------------------

%=======================================================================
% - Get files and factor indices
%=======================================================================
% DbaTime=D.b.aTime; Dn=D.n; DsF=D.sF;
if nargin<4, ext = 'any'; else, ext = varargin{4}; end
if nargin<3, DbaTime = 1; else, DbaTime = varargin{3}; end
if nargin<2, Dn  = [Inf,Inf,Inf,Inf]; else, Dn=varargin{2}; end
if nargin<1, DsF = {'Fac1','Fac2','Fac3','Fac4'}; else, DsF=varargin{1}; end

%-Initialise variables
%-----------------------------------------------------------------------
i4 = [];          % factor 4 index (usually group)
i3 = [];          % factor 3 index (usually subject), per f4
i2 = [];          % factor 2 index (usually condition), per f3/f4
i1 = [];          % factor 1 index (usually replication), per f2/f3/f4
P  = {};          % cell array of string filenames

nV = 1;

%-Accrue filenames and factor level indicator vectors
%-----------------------------------------------------------------------
bMV = nV>1;
if isinf(Dn(4)), n4 = spm_input(['#',DsF{4},'''s'],'+1','n1');
      else, n4 = Dn(4); end
bL4 = n4>1;

ti2 = '';
GUIpos = spm_input('!NextPos');
for j4  = 1:n4
    spm_input('!SetNextPos',GUIpos);
    sF4P=''; if bL4, sF4P=[DsF{4},' ',int2str(j4),': ']; end
    if isinf(Dn(3)), n3=spm_input([sF4P,'#',DsF{3},'''s'],'+1','n1');
          else, n3 = Dn(3); end
    bL3 = n3>1;
    
    if DbaTime & Dn(2)>1
      %disp('NB:selecting in time order - manually specify conditions')
      %-NB: This means f2 levels might not be 1:n2
      GUIpos2 = spm_input('!NextPos');
      for j3 = 1:n3
          sF3P=''; if bL3, sF3P=[DsF{3},' ',int2str(j3),': ']; end
          str = [sF4P,sF3P];
          tP  = {};
          n21 = Dn(2)*Dn(1);
          for v=1:nV
            vstr=''; if bMV, vstr=sprintf(' (var-%d)',v); end
            ttP = cellstr(spm_select(n21,ext,[str,'select images',vstr]));
            n21 = length(ttP);
            tP  = [tP,ttP];
          end
          ti2 = spm_input([str,' ',DsF{2},'?'],GUIpos2,'c',ti2',n21,Dn(2));
          %-Work out i1 & check
          [tl2,null,j] = unique(ti2);
          tn1 = zeros(size(tl2)); ti1 = zeros(size(ti2));
          for i=1:length(tl2)
                tn1(i)=sum(j==i); ti1(ti2==tl2(i))=1:tn1(i); end
          if isfinite(Dn(1)) & any(tn1~=Dn(1))
            %-#i1 levels mismatches specification in Dn(1)
            error(sprintf('#%s not %d as pre-specified',DsF{1},Dn(1)))
          end
          P  = [P;tP];
          i4 = [i4; j4*ones(n21,1)];
          i3 = [i3; j3*ones(n21,1)];
          i2 = [i2; ti2];
          i1 = [i1; ti1];
      end

    else

      if isinf(Dn(2))
          n2 = spm_input([sF4P,'#',DsF{2},'''s'],'+1','n1');
      else
          n2 = Dn(2);
      end
      bL2 = n2>1;

      if n2==1 & Dn(1)==1 %-single scan per f3 (subj)
          %disp('NB:single scan per f3')
          str = [sF4P,'select images, ',DsF{3},' 1-',int2str(n3)];
          tP = {};
          for v=1:nV
            vstr=''; if bMV, vstr=sprintf(' (var-%d)',v); end
            ttP = cellstr(spm_select(n3,ext,[str,vstr]));
            tP = [tP,ttP];
          end
          P   = [P;tP];
          i4  = [i4; j4*ones(n3,1)];
          i3  = [i3; [1:n3]'];
          i2  = [i2; ones(n3,1)];
          i1  = [i1; ones(n3,1)];
      else
          %-multi scan per f3 (subj) case
          %disp('NB:multi scan per f3')
          for j3 = 1:n3
            sF3P=''; if bL3, sF3P=[DsF{3},' ',int2str(j3),': ']; end
            if Dn(1)==1
                  %-No f1 (repl) within f2 (cond)
                  %disp('NB:no f1 within f2')
                  str = [sF4P,sF3P,'select images: ',DsF{2},...
                         ' 1-',int2str(n2)];
                  tP = {};
                  for v=1:nV
                        vstr=''; if bMV, vstr=sprintf(' (var-%d)',v); end
                        ttP = cellstr(spm_select(n2,ext,[str,vstr]));
                        tP = [tP,ttP];
                  end
                  P   = [P;tP];
                  i4  = [i4; j4*ones(n2,1)];
                  i3  = [i3; j3*ones(n2,1)];
                  i2  = [i2; [1:n2]'];
                  i1  = [i1; ones(n2,1)];
            else
                %-multi f1 (repl) within f2 (cond)
                %disp('NB:f1 within f2')
                for j2  = 1:n2
                  sF2P='';
                  if bL2, sF2P=[DsF{2},' ',int2str(j2),': ']; end
                  str = [sF4P,sF3P,sF2P,' select images...'];
                  tP  = {};
                  n1  = Dn(1);
                  for v=1:nV
                        vstr=''; if bMV, vstr=sprintf(' (var-%d)',v); end
                        ttP = cellstr(spm_select(n1,ext,[str,vstr]));
                        n1  = length(ttP);
                        tP  = [tP,ttP];
                  end
                  P   = [P;tP];
                  i4  = [i4; j4*ones(n1,1)];
                  i3  = [i3; j3*ones(n1,1)];
                  i2  = [i2; j2*ones(n1,1)];
                  i1  = [i1; [1:n1]'];
                end                         % (for j2)
            end                             % (if Dn(1)==1)
          end                                 % (for j3)
      end                                     % (if  n2==1 &...)
    end                                         % (if DbaTime & Dn(2)>1)
end                                             % (for j4)
varargout = {P,[i1,i2,i3,i4]};

% -------------------------------------------------------------------------------------
function out_str = replace_char(in_str)
% replace spaces with "_" and characters like "<" or ">" with "gt" or "lt"

str_num = deblank(in_str);

str_num(findstr(str_num,' ')) = '_';
strpos = findstr(str_num,' > ');
if ~isempty(strpos), str_num = [str_num(1:strpos-1) '_gt_' str_num(strpos+1:end)]; end
strpos = findstr(str_num,' < ');
if ~isempty(strpos), str_num = [str_num(1:strpos-1) '_lt_' str_num(strpos+1:end)]; end
strpos = findstr(str_num,'>');
if ~isempty(strpos), str_num = [str_num(1:strpos-1) 'gt' str_num(strpos+1:end)]; end
strpos = findstr(str_num,'<');
if ~isempty(strpos), str_num = [str_num(1:strpos-1) 'lt' str_num(strpos+1:end)]; end
str_num = spm_str_manip(str_num,'v');

out_str = str_num;
return