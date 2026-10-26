function make_reference(outdir)
% Reference fixtures for the CAT-Surface port of SPM's shooting primitives
% FORMAT make_reference(outdir)
% outdir - directory for the fixtures (created if needed);
%          default: fullfile(tempdir,'spm_shoot_ref')
%
% Runs spm_diffeo, spm_bsplinc and spm_bsplins of SPM26 on seeded, smooth
% random inputs at a small anisotropic size (31x27x23, 1.5x1.2x1.8 mm) and
% at the size CAT12's shooting uses (113x137x113, 1.5 mm), with the
% regularisation CAT12 uses (prm = [vx, rparam*sched*prod(vx)]), and writes
% every input and output as a raw little-endian array:
%   <outdir>/<name>.bin   column-major data, float32 (f32) or float64 (f64)
%   <outdir>/manifest.txt one line per array: name type ndims dim1 ... dimN
%   <outdir>/timing.txt   run time of the mex calls on the real-size data
% tests/spm_shoot_reference/compare_reference.c reads them back.
%
% Example (SPM26 must be on the path):
%   matlab -batch "addpath('/Users/gaser/spm/spm12'); ...
%       addpath('tests/spm_shoot_reference'); make_reference('/tmp/ref')"
%__________________________________________________________________________

if nargin < 1, outdir = fullfile(tempdir, 'spm_shoot_ref'); end
if ~exist(outdir, 'dir'), mkdir(outdir); end
if ~exist('spm_diffeo', 'file')
    error('make_reference:spm', 'SPM26 (spm_diffeo) is not on the path.');
end
if ~exist('spm_shoot3d', 'file')
    addpath(fullfile(spm('Dir'), 'toolbox', 'Shoot'));
end

mf = fopen(fullfile(outdir, 'manifest.txt'), 'w');
tf = fopen(fullfile(outdir, 'timing.txt'), 'w');
c1 = onCleanup(@() fclose(mf));
c2 = onCleanup(@() fclose(tf));
put = @(name, x) save_arr(mf, outdir, name, x);

fprintf(tf, '# %s, %s, maxNumCompThreads=%d\n', spm('Ver'), version, maxNumCompThreads);
rng(20261026, 'twister');
spm_diffeo('boundary', 0);

rparam = [1e-4 0.001 0.2 0.05 0.2];   % spm_shoot_defaults: rparam
bs     = [2 2 2  1 1 1];              % spm_shoot_defaults: bs_args
cases  = {'small', [31 27 23],   [1.5 1.2 1.8], [32 4 1]
          'real',  [113 137 113], [1.5 1.5 1.5], [32 4 1]};

for ic = 1:size(cases, 1)
    pf = cases{ic, 1};
    d  = cases{ic, 2};
    vx = cases{ic, 3};
    sc = cases{ic, 4};
    timed = strcmp(pf, 'real');
    p  = @(s) [pf '_' s];
    fprintf('%s: %dx%dx%d\n', pf, d);

    % ---------------------------------------------------------------------
    % Inputs: velocity (|v| up to 3 voxels), PSD Hessian, gradient, tissues
    % ---------------------------------------------------------------------
    v = zeros([d 3], 'single');
    for k = 1:3, v(:,:,:,k) = single(3*smooth_unit(d, 4)); end
    G = zeros([d 9]);
    for k = 1:9, G(:,:,:,k) = 0.5*smooth_unit(d, 2); end
    A = zeros([d 6], 'single');
    A(:,:,:,1) = G(:,:,:,1).^2 + G(:,:,:,2).^2 + G(:,:,:,3).^2;
    A(:,:,:,2) = G(:,:,:,4).^2 + G(:,:,:,5).^2 + G(:,:,:,6).^2;
    A(:,:,:,3) = G(:,:,:,7).^2 + G(:,:,:,8).^2 + G(:,:,:,9).^2;
    A(:,:,:,4) = G(:,:,:,1).*G(:,:,:,4) + G(:,:,:,2).*G(:,:,:,5) + G(:,:,:,3).*G(:,:,:,6);
    A(:,:,:,5) = G(:,:,:,1).*G(:,:,:,7) + G(:,:,:,2).*G(:,:,:,8) + G(:,:,:,3).*G(:,:,:,9);
    A(:,:,:,6) = G(:,:,:,4).*G(:,:,:,7) + G(:,:,:,5).*G(:,:,:,8) + G(:,:,:,6).*G(:,:,:,9);
    clear G
    b = zeros([d 3], 'single');
    for k = 1:3, b(:,:,:,k) = single(0.1*smooth_unit(d, 1) + 0.01*randn(d)); end
    f = zeros([d 3], 'single');
    for k = 1:3, f(:,:,:,k) = single(0.5 + 0.5*smooth_unit(d, 2)); end
    put(p('v'), v); put(p('A'), A); put(p('b'), b); put(p('f'), f);

    % ---------------------------------------------------------------------
    % vel2mom, fmg (FMG path, cyc 3 nit 2 as in spm_shoot_update), trapprox
    % ---------------------------------------------------------------------
    for s = sc
        prm = [vx, rparam*s*prod(vx)];
        put(p(sprintf('prm_s%d', s)), prm);
        t = tic; m = spm_diffeo('vel2mom', v, prm); tv = toc(t);
        put(p(sprintf('vel2mom_s%d', s)), m);
        t = tic; u = spm_diffeo('fmg', A, b, [prm 3 2]); tf_ = toc(t);
        put(p(sprintf('fmg_s%d', s)), u);
        put(p(sprintf('trapprox_s%d', s)), spm_diffeo('trapprox', A, prm));
        if timed
            fprintf(tf, 'vel2mom_s%d %.6f\nfmg_s%d %.6f\n', s, tv, s, tf_);
        end
    end
    prm = [vx, rparam*4*prod(vx)];

    % fmg with a starting estimate (V-cycle path), and the other relaxers
    u0 = spm_diffeo('fmg', A, b, [prm 1 1]);
    put(p('fmg_u0'), u0);
    put(p('fmg_vcyc'), spm_diffeo('fmg', A, b, [prm 3 2], u0));
    if strcmp(pf, 'small')
        variants = {'me',  [1e-4 0.001 0   0    0  ]
                    'be',  [1e-4 0.001 0.2 0    0  ]
                    'le',  [1e-4 0     0   0.05 0.2]};
        for iv = 1:size(variants, 1)
            pv = [vx, variants{iv, 2}*4*prod(vx)];
            put(p(['prm_' variants{iv, 1}]), pv);
            put(p(['fmg_' variants{iv, 1}]), spm_diffeo('fmg', A, b, [pv 3 2]));
            put(p(['mom2vel_' variants{iv, 1}]), spm_diffeo('mom2vel', b, [pv 3 2]));
            put(p(['vel2mom_' variants{iv, 1}]), spm_diffeo('vel2mom', v, pv));
        end
        put(p('mom2vel_all'), spm_diffeo('mom2vel', b, [prm 3 2]));
        put(p('kernel_be'), spm_diffeo('kernel', d, [vx, [1e-4 0.001 0.2 0 0]*4*prod(vx)]));
    end
    put(p('kernel'), spm_diffeo('kernel', d, prm));

    % ---------------------------------------------------------------------
    % Shooting (spm_shoot3d, 8 steps) for a realistic phi / Jphi
    % ---------------------------------------------------------------------
    [phi, Jphi] = spm_shoot3d(v, prm, [8 2 2]);
    put(p('phi'), phi); put(p('Jphi'), Jphi);

    t = tic; [dp, dJ] = spm_diffeo('smalldef', v, 1/8); ts = toc(t);
    put(p('smalldef_y'), dp); put(p('smalldef_J'), dJ);
    [dpm, dJm] = spm_diffeo('smalldef', v, -1/8);
    put(p('smalldefm_y'), dpm); put(p('smalldefm_J'), dJm);
    put(p('smalldef1_y'), spm_diffeo('smalldef', v, 0.37));

    t = tic; [phi2, Jphi2] = spm_diffeo('comp', dp, phi, dJ, Jphi); tc = toc(t);
    put(p('comp_y'), phi2); put(p('comp_J'), Jphi2);
    put(p('comp2_y'), spm_diffeo('comp', phi, dpm));

    m0 = spm_diffeo('vel2mom', v, prm);
    put(p('m0'), m0);
    t = tic; mt = spm_diffeo('pushg', m0, phi, Jphi); tg = toc(t);
    put(p('pushg'), mt);
    put(p('pushg_noJ'), spm_diffeo('pushg', m0, phi));

    t = tic; f1 = spm_diffeo('pullc', f, phi); tp = toc(t);
    put(p('pullc'), f1);
    put(p('pull'), spm_diffeo('pull', f, phi));
    [f0, c0] = spm_diffeo('pushc', f1, phi);
    put(p('pushc_f'), f0); put(p('pushc_c'), c0);
    [f0, c0] = spm_diffeo('push', f1, phi);
    put(p('push_f'), f0); put(p('push_c'), c0);

    t = tic; dt = spm_diffeo('det', Jphi); td = toc(t);
    put(p('det'), dt);
    t = tic; jd = spm_diffeo('def2det', phi); tj = toc(t);
    put(p('def2det'), jd);
    if strcmp(pf, 'small'), put(p('def2jac'), spm_diffeo('def2jac', phi)); end

    % invdef as CAT12 calls it: output grid of another size, affine maps
    M1 = [1.02 0.03 -0.01 -1.5; -0.02 0.97 0.02 2.0; 0.01 -0.03 1.04 -0.7; 0 0 0 1];
    M2 = [0.99 -0.01 0.02 0.8; 0.015 1.03 -0.02 -1.1; -0.01 0.02 0.98 0.4; 0 0 0 1];
    di = d + [2 -3 1];
    put(p('invdef_M1'), M1); put(p('invdef_M2'), M2); put(p('invdef_dim'), di);
    t = tic; iy = spm_diffeo('invdef', phi, di, M1, M2); ti = toc(t);
    put(p('invdef'), iy);
    % (not spm_diffeo('invdef',phi): with one argument the gateway leaves its
    % two matrices uninitialised)
    put(p('invdef_id'), spm_diffeo('invdef', phi, d, eye(4), eye(4)));

    % multigrid transfer operators
    put(p('restrict'), spm_diffeo('restrict', f(:,:,:,1)));
    put(p('resize_dn'), spm_diffeo('resize', f(:,:,:,1), ceil(d/2)));
    put(p('resize_up'), spm_diffeo('resize', spm_diffeo('restrict', f(:,:,:,1)), d));

    % NaN branches of pullc / comp / pushg (small only)
    if strcmp(pf, 'small')
        phin = phi;
        phin(5,6,7,1) = NaN; phin(10,3,2,2) = Inf; phin(1,1,1,3) = -Inf;
        put(p('phin'), phin);
        put(p('pullc_nan'), spm_diffeo('pullc', f, phin));
        [cy, cJ] = spm_diffeo('comp', dp, phin, dJ, Jphi);
        put(p('comp_nan_y'), cy); put(p('comp_nan_J'), cJ);
        put(p('pushg_nan'), spm_diffeo('pushg', m0, phin, Jphi));
    end

    % ---------------------------------------------------------------------
    % B-splines as CAT12 uses them: coefficients of log template, degree 2,
    % wrapped; sampled at the identity grid with derivatives
    % ---------------------------------------------------------------------
    g  = single(0.02 + 0.96*(0.5 + 0.5*smooth_unit(d, 2)));
    lg = log(g);            % saved: C's logf need not match MATLAB's log
    put(p('bs_lg'), lg);
    t = tic; c = spm_bsplinc(lg, bs); tbc = toc(t);
    put(p('bs_c'), c);
    [x0, x1, x2] = ndgrid(1:d(1), 1:d(2), 1:d(3));
    t = tic; [s0, s1, s2, s3] = spm_bsplins(c, x0, x1, x2, bs); tbs = toc(t);
    put(p('bs_id_f'), s0); put(p('bs_id_d0'), s1); put(p('bs_id_d1'), s2); put(p('bs_id_d2'), s3);
    clear x0 x1 x2

    % the 3-slice trick of spm_shoot_update/mnom_derivs gives the same values
    [id{1}, id{2}, id{3}] = ndgrid(1:d(1), 1:d(2), ceil(bs(3)/2)+1);
    zz = [1 2 round(d(3)/2) d(3)];
    maxd = 0;
    mx = @(a, b) max(abs(a(:) - b(:)));
    for z = zz
        ind = z + (-ceil(bs(3)/2):ceil(bs(3)/2));
        ind = rem(ind-1+d(3), d(3)) + 1;
        [g1, e1, e2, e3] = spm_bsplins(double(c(:,:,ind)), id{:}, bs);
        maxd = max([maxd, mx(g1, s0(:,:,z)), mx(e1, s1(:,:,z)), ...
                  mx(e2, s2(:,:,z)), mx(e3, s3(:,:,z))]);
    end
    fprintf(tf, '# mnom_derivs 3-slice sampling vs full-volume sampling: max abs diff %g\n', maxd);

    % arbitrary points, outside the volume, and non-finite coordinates
    n  = 4000;
    xr = [rand(n,1)*(d(1)+12)-6, rand(n,1)*(d(2)+12)-6, rand(n,1)*(d(3)+12)-6];
    xr(1:3,:) = [NaN 2 2; 3 Inf 3; 4 4 -Inf];
    xr(4,:)   = [1-0.05 d(2)+0.05 d(3)];    % at the TINY margin
    xr(5,:)   = [1-0.051 2 2];              % just outside it
    put(p('bs_xr'), xr);
    [r0, r1, r2, r3] = spm_bsplins(c, xr(:,1), xr(:,2), xr(:,3), bs);
    put(p('bs_xr_f'), r0); put(p('bs_xr_d0'), r1); put(p('bs_xr_d1'), r2); put(p('bs_xr_d2'), r3);
    put(p('bs_xr_fonly'), spm_bsplins(c, xr(:,1), xr(:,2), xr(:,3), bs));

    if strcmp(pf, 'small')
        % other settings: mirrored degree 3, mixed degree/boundaries, NaN voxels
        gn = log(g);
        gn(3,4,5) = NaN; gn(7,7,7) = Inf;
        put(p('bs_gn'), gn);
        bsm = [3 3 3 0 0 0];
        cm  = spm_bsplinc(gn, bsm);
        put(p('bs_cm'), cm);
        [r0, r1, r2, r3] = spm_bsplins(cm, xr(:,1), xr(:,2), xr(:,3), bsm);
        put(p('bs_m_f'), r0); put(p('bs_m_d0'), r1); put(p('bs_m_d1'), r2); put(p('bs_m_d2'), r3);
        bsx = [4 1 2 1 0 1];
        cx  = spm_bsplinc(double(gn), bsx);
        put(p('bs_cx'), cx);
        [r0, r1, r2, r3] = spm_bsplins(cx, xr(:,1), xr(:,2), xr(:,3), bsx);
        put(p('bs_x_f'), r0); put(p('bs_x_d0'), r1); put(p('bs_x_d1'), r2); put(p('bs_x_d2'), r3);
    end

    if timed
        fprintf(tf, 'smalldef_jac1 %.6f\ncomposition_jacobian %.6f\npushc_grads %.6f\n', ts, tc, tg);
        fprintf(tf, 'pullc %.6f\ndeterminant %.6f\ndef2det %.6f\ninvdef %.6f\n', tp, td, tj, ti);
        fprintf(tf, 'bsplinc %.6f\nbsplins_deriv %.6f\n', tbc, tbs);
    end
end
fprintf('done: %s\n', outdir);
end

%==========================================================================
function x = smooth_unit(d, sig)
% Gaussian-smoothed white noise scaled to max |x| = 1
r = ceil(3*sig);
k = exp(-(-r:r).^2/(2*sig^2));
k = k/sum(k);
x = randn(d);
x = convn(x, reshape(k, [], 1, 1), 'same');
x = convn(x, reshape(k, 1, [], 1), 'same');
x = convn(x, reshape(k, 1, 1, []), 'same');
x = x/max(abs(x(:)));
end

%==========================================================================
function save_arr(mf, outdir, name, x)
% Write x as raw little-endian data and add a manifest line
if isa(x, 'single')
    t = 'f32'; prec = 'float32';
elseif isa(x, 'double')
    t = 'f64'; prec = 'float64';
else
    error('make_reference:class', 'Unsupported class %s for %s.', class(x), name);
end
fid = fopen(fullfile(outdir, [name '.bin']), 'w', 'ieee-le');
fwrite(fid, x, prec);
fclose(fid);
d = size(x);
fprintf(mf, '%s %s %d', name, t, numel(d));
fprintf(mf, ' %d', d);
fprintf(mf, '\n');
end
