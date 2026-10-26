function make_reference(outdir)
%MAKE_REFERENCE Reference data for the C ports in Lib/CAT_Cat12Utils.c.
%
%  make_reference(outdir)
%
%  Runs the CAT12/CAT26 helpers that the safeguards of run_Shooting
%  (cat_main_registration.m) call -- cat_vol_smooth3X, cat_vol_div,
%  cat_vol_morph ('l', 'd', 'o', 'ldc', 'ldo'), cat_vol_localstat,
%  cat_vol_laplace3R and cat_vbdist -- on synthetic inputs at the five
%  registration grid sizes and writes raw inputs and outputs (float32 or
%  uint8, MATLAB order) plus manifest.txt to outdir.  compare_reference.c
%  reads the manifest, runs the C ports and reports the differences.
%
%  Inputs: the shooting template (GM/WM) resampled to each grid, with and
%  without an enlarged ventricle (hydrocephalus) and a tunnel from it to the
%  outside; the steps of run_Shooting's ventricle detection on those; smooth
%  random fields, velocity- and Jacobian-like fields with hot spots, NaN and
%  Inf; binary blobs with holes, several components, objects touching the
%  border, ties, empty and full volumes.
%
%  Needs SPM (spm12/spm26) and CAT12/CAT26 on the MATLAB path, e.g.
%    matlab -batch "addpath('~/spm/spm12'); addpath('~/spm/spm12/toolbox/CAT'); \
%                   cd tests/cat12_utils_reference; make_reference('/tmp/ref')"
%
%  Manifest line: name func nx ny nz p1 p2 p3 p4 input mask output [output2]
%  ('-' for no mask).
% ______________________________________________________________________
% Christian Gaser
% Structural Brain Mapping Group, Jena University Hospital
% ______________________________________________________________________

  if nargin < 1, outdir = fullfile(tempdir, 'cat12_utils_ref'); end
  if ~exist(outdir, 'dir'), mkdir(outdir); end
  fid = fopen(fullfile(outdir, 'manifest.txt'), 'w');
  ctx = struct('fid', fid, 'dir', outdir);

  G        = [67 82 67; 75 91 75; 84 102 84; 96 117 96; 113 137 113];
  tempres2 = [2.5 2.25 2.0 1.75 1.5];
  rred     = [0.6 0.667 0.75 0.857 1.0];
  sched    = [0.2 0.45 0.5 0.8 1.3 2.0];   % max(0.2,min(2,sd.sched(it)/8))

  catdir = fileparts(which('cat_vol_morph'));
  Vt  = spm_vol(fullfile(catdir, 'templates_MNI152NLin2009cAsym', 'Template_0_GS.nii'));
  TGM = spm_read_vols(Vt(1));
  TWM = spm_read_vols(Vt(2));

  for ti = 1:size(G, 1)
    d  = G(ti, :);
    t  = tempres2(ti);
    vr = 3 / t;
    rng(100 + ti, 'twister');

    % ------------------------------------------------------------------
    % brain-like inputs and the ventricle detection of run_Shooting (n1 == 2)
    % ------------------------------------------------------------------
    for variant = 1:3
      [f1, f2] = brain(TGM, TWM, d, variant);
      tag = sprintf('L%d_b%d', ti, variant);

      in  = single(f2 + f1);
      msk = (f2 + f1 + 1/2) > 0;
      Ygwm = cat_vol_localstat(in, msk, 1, 1);
      addcase(ctx, [tag '_localstat'], 'localstat', d, 1, in, msk, Ygwm);

      Yb = cat_vol_morph(Ygwm, 'ldc', 12/t);
      addcase(ctx, [tag '_ldc'], 'ldc', d, [12/t 1 1 1], uint8(Ygwm > 0.5), [], uint8(Yb > 0.5));

      Ybd = cat_vol_smooth3X(Yb, 6/t);
      addcase(ctx, [tag '_smooth_Yb'], 'smooth3X', d, 6/t, single(Yb), [], Ybd);
      Ybd = Ybd / max(Ybd(:));

      Yven0 = (Yb - Ygwm) > 0.9 & Ybd > 0.9;
      Yven1 = cat_vol_morph(Yven0, 'ldo', 3, t);
      addcase(ctx, [tag '_ldo'], 'ldo', d, [3 t t t], uint8(Yven0), [], uint8(Yven1 > 0.5));

      ndil = max(1, 1/rred(ti));
      Yven2 = cat_vol_morph(Yven1, 'd', ndil);
      addcase(ctx, [tag '_dilate'], 'dilate', d, [ndil 1 1 1], uint8(Yven1 > 0.5), [], uint8(Yven2));

      Yven3 = cat_vol_morph(Yven2, 'l');
      addcase(ctx, [tag '_lab'], 'lab', d, [], uint8(Yven2), [], uint8(Yven3));

      Ys = cat_vol_smooth3X(Yven3, 4/t);
      addcase(ctx, [tag '_smooth_Yven'], 'smooth3X', d, 4/t, single(Yven3), [], Ys);

      Ys = cat_vol_smooth3X(single(min(1, Ys * 1.5)), 2/t);
      addcase(ctx, [tag '_smooth_Yven2'], 'smooth3X', d, 2/t, single(min(1, cat_vol_smooth3X(Yven3, 4/t) * 1.5)), [], Ys);
    end

    % ------------------------------------------------------------------
    % hot-spot correction: dt, u, mxu
    % ------------------------------------------------------------------
    dt = exp(0.3 * smoothfield(d, 6));
    spots = rand(d) > 0.9995;
    dt(spots) = dt(spots) * 40;
    dt(rand(d) > 0.9998) = 1/30;
    dt = single(dt);
    Sdt = cat_vol_smooth3X(dt, vr);
    addcase(ctx, sprintf('L%d_smooth_dt', ti), 'smooth3X', d, vr, dt, [], Sdt);
    mxu = dt ./ Sdt > 1.5;
    Smxu = cat_vol_smooth3X(mxu, 2*vr);
    addcase(ctx, sprintf('L%d_smooth_mxu', ti), 'smooth3X', d, 2*vr, single(mxu), [], Smxu);

    u = 3 * smoothfield(d, 8);
    Ul = cat_vol_laplace3R(u/10, Smxu > 0.01, 0.1) * 10;
    addcase(ctx, sprintf('L%d_laplace_0.1', ti), 'laplace', d, 0.1, u, Smxu > 0.01, Ul);
    mxu2 = cat_vol_smooth3X(dt ./ cat_vol_smooth3X(dt, 1*vr) > 1.5, 2*vr) > 0.01;
    Ul2 = cat_vol_laplace3R(u/10, mxu2, 0.02) * 10;
    addcase(ctx, sprintf('L%d_laplace_0.02', ti), 'laplace', d, 0.02, u, mxu2, Ul2);

    addcase(ctx, sprintf('L%d_smooth_u_1', ti), 'smooth3X', d, 1, u, [], cat_vol_smooth3X(u, 1));
    for s = sched
      if ti == 1 || ti == 5 || any(abs(s - sched(mod(ti, numel(sched)) + 1)) < eps)
        addcase(ctx, sprintf('L%d_smooth_u_%g', ti, s), 'smooth3X', d, s, u, [], cat_vol_smooth3X(u, s));
      end
    end
    for s = [2/t 3/t 6/t]
      addcase(ctx, sprintf('L%d_smooth_noise_%0.4f', ti, s), 'smooth3X', d, s, ...
        single(randn(d)), [], []);
    end

    % divergence of a velocity component, as in the hot-spot loop
    addcase(ctx, sprintf('L%d_div', ti), 'div', d, [], u, [], cat_vol_div(u));
    [f1, f2] = brain(TGM, TWM, d, 1);
    ub = u .* single(cat_vol_smooth3X(single(f1 + f2), 2) > 0.1);
    addcase(ctx, sprintf('L%d_div_brain', ti), 'div', d, [], ub, [], cat_vol_div(ub));
    un = u; un(smoothfield(d, 4) > 0.55) = NaN;
    addcase(ctx, sprintf('L%d_div_nan', ti), 'div', d, [], un, [], cat_vol_div(un));

    % ------------------------------------------------------------------
    % binary blobs: several components, holes, objects at the border
    % ------------------------------------------------------------------
    B1 = smoothfield(d, 5) > 0.25;
    B2 = blobs(d);
    for b = 1:2
      if b == 1, B = B1; else, B = B2; end
      tag = sprintf('L%d_B%d', ti, b);
      addcase(ctx, [tag '_lab'], 'lab', d, [], uint8(B), [], uint8(cat_vol_morph(B, 'l')));
      for ndil = unique([1 max(1, 1/rred(ti)) 1.5 2])
        addcase(ctx, sprintf('%s_dilate_%0.4f', tag, ndil), 'dilate', d, [ndil 1 1 1], uint8(B), [], ...
          uint8(cat_vol_morph(B, 'd', ndil)));
      end
      addcase(ctx, [tag '_open'], 'open', d, [1 1 1 1], uint8(B), [], uint8(cat_vol_morph(B, 'o', 1)));
      addcase(ctx, [tag '_ldc'], 'ldc', d, [12/t 1 1 1], uint8(B), [], uint8(cat_vol_morph(B, 'ldc', 12/t) > 0.5));
      addcase(ctx, [tag '_ldo'], 'ldo', d, [3 t t t], uint8(B), [], uint8(cat_vol_morph(B, 'ldo', 3, t) > 0.5));
      addcase(ctx, [tag '_labmulti'], 'labmulti', d, [10 0.1], uint8(B), [], ...
        uint8(cat_vol_morph(B, 'l', [10 0.1]) > 0.5));
      [D, I] = cat_vbdist(single(B), true(d), [t t t]);
      addcase(ctx, [tag '_vbdist'], 'vbdist', d, [t t t], single(B), [], D, uint32(I));
    end
  end

  % --------------------------------------------------------------------
  % every ldc/ldo distance at every grid, on the blobs of the smallest grid
  % --------------------------------------------------------------------
  d  = G(1, :);
  rng(7, 'twister');
  B1 = smoothfield(d, 5) > 0.25;
  B2 = blobs(d);
  for ti = 1:numel(tempres2)
    t = tempres2(ti);
    for b = 1:2
      if b == 1, B = B1; else, B = B2; end
      addcase(ctx, sprintf('X_B%d_ldc_%0.4f', b, 12/t), 'ldc', d, [12/t 1 1 1], uint8(B), [], ...
        uint8(cat_vol_morph(B, 'ldc', 12/t) > 0.5));
      addcase(ctx, sprintf('X_B%d_ldo_%0.2f', b, t), 'ldo', d, [3 t t t], uint8(B), [], ...
        uint8(cat_vol_morph(B, 'ldo', 3, t) > 0.5));
    end
  end
  % other open sizes and anisotropic voxels
  for ndil = [1 1.4 2 3]
    addcase(ctx, sprintf('X_B2_open_%g', ndil), 'open', d, [ndil 1 1 1], uint8(B2), [], uint8(cat_vol_morph(B2, 'o', ndil)));
  end
  addcase(ctx, 'X_B1_dilate_aniso', 'dilate', d, [2 1 1.5 2.5], uint8(B1), [], uint8(cat_vol_morph(B1, 'd', 2, [1 1.5 2.5])));
  addcase(ctx, 'X_B1_open_aniso', 'open', d, [2 1 1.5 0.8], uint8(B1), [], uint8(cat_vol_morph(B1, 'o', 2, [1 1.5 0.8])));
  addcase(ctx, 'X_B2_ldc_aniso', 'ldc', d, [4 1 1.25 2], uint8(B2), [], uint8(cat_vol_morph(B2, 'ldc', 4, [1 1.25 2]) > 0.5));
  addcase(ctx, 'X_B2_ldo_aniso', 'ldo', d, [3 1.3 1.6 2.2], uint8(B2), [], uint8(cat_vol_morph(B2, 'ldo', 3, [1.3 1.6 2.2]) > 0.5));
  addcase(ctx, 'X_B2_ldc_7.3', 'ldc', d, [7.3 1 1 1], uint8(B2), [], uint8(cat_vol_morph(B2, 'ldc', 7.3) > 0.5));
  addcase(ctx, 'X_B1_labmulti_abs', 'labmulti', d, [4 -50], uint8(B1), [], uint8(cat_vol_morph(B1, 'l', [4 -50]) > 0.5));
  addcase(ctx, 'X_B1_labmulti_abs2', 'labmulti', d, [6 200], uint8(B1), [], uint8(cat_vol_morph(B1, 'l', [6 200]) > 0.5));

  % distance transform: non-dyadic and anisotropic voxels, masks
  for vx = {[1 1 1], [1.5 1.5 1.5], [1 1.3 2], [0.7 0.7 0.7]}
    v = vx{1};
    [D, I] = cat_vbdist(single(B2), true(d), v);
    addcase(ctx, sprintf('X_vbdist_%g_%g_%g', v), 'vbdist', d, v, single(B2), [], D, uint32(I));
  end
  M = smoothfield(d, 6) > -0.2;
  [D, I] = cat_vbdist(single(B1), M, [1 1 1]);
  addcase(ctx, 'X_vbdist_mask', 'vbdist', d, [1 1 1], single(B1), M, D, uint32(I));
  [D, I] = cat_vbdist(single(B1), M, [1.75 1.75 1.75]);
  addcase(ctx, 'X_vbdist_mask_1.75', 'vbdist', d, [1.75 1.75 1.75], single(B1), M, D, uint32(I));
  Vk = zeros(11, 11, 11, 'single'); Vk(6, 6, 6) = 1;
  [D, I] = cat_vbdist(Vk, true(size(Vk)), [1.5 1.5 1.5]);
  addcase(ctx, 'X_vbdist_kernel', 'vbdist', size(Vk), [1.5 1.5 1.5], Vk, [], D, uint32(I));

  % thin volumes: the neighbour test of cat_vbdist/cat_vol_gradient3 wraps
  for sz = {[1 7 5], [2 6 5], [6 1 5], [5 2 4], [3 3 3]}
    d = sz{1};
    B = rand(d) > 0.6; B(1) = true;
    for v = {[1 1 1], [1.5 1.5 1.5]}
      [D, I] = cat_vbdist(single(B), true(d), v{1});
      addcase(ctx, sprintf('T_vbdist_%dx%dx%d_%g', d, v{1}(1)), 'vbdist', d, v{1}, single(B), [], D, uint32(I));
    end
    X = single(1 + 2 * rand(d));
    addcase(ctx, sprintf('T_div_%dx%dx%d', d), 'div', d, [], X, [], cat_vol_div(X));
    addcase(ctx, sprintf('T_lab_%dx%dx%d', d), 'lab', d, [], uint8(B), [], uint8(cat_vol_morph(B, 'l')));
  end
  d = G(1, :);

  % localstat with zeros, NaN and a partial mask; div with Inf
  Yl = single(max(0, smoothfield(d, 4)));
  Yl(rand(d) > 0.999) = NaN;
  Yl(rand(d) > 0.9995) = Inf;
  M = smoothfield(d, 6) > -0.3;
  addcase(ctx, 'X_localstat_mask', 'localstat', d, 1, Yl, M, cat_vol_localstat(Yl, M, 1, 1));
  addcase(ctx, 'X_localstat_nb2', 'localstat', d, 2, Yl, M, cat_vol_localstat(Yl, M, 2, 1));
  ui = 3 * smoothfield(d, 8);
  ui(30:31, 40, 33) = Inf; ui(5, 6, 7) = -Inf; ui(1, 1, 1) = Inf;
  addcase(ctx, 'X_div_inf', 'div', d, [], ui, [], cat_vol_div(ui));
  un = 3 * smoothfield(d, 8); un(:, :, 1:20) = NaN;
  addcase(ctx, 'X_div_nanslab', 'div', d, [], un, [], cat_vol_div(un));

  % special volumes
  E = false(d);
  F = true(d);
  C = false(d); C(1, 1, 1) = true;
  P = false(d); P(34, 40, 30) = true;
  T = false(d); T(10:14, 10:14, 10:14) = true; T(40:44, 50:54, 40:44) = true;   % tie
  W = false(d); W(:, :, 1) = true; W(20:30, 20:30, 20:30) = true;                % plane at the border
  specials = {'empty', E; 'full', F; 'corner', C; 'point', P; 'tie', T; 'wall', W};
  for k = 1:size(specials, 1)
    B = specials{k, 2}; nm = specials{k, 1};
    addcase(ctx, ['S_' nm '_lab'], 'lab', d, [], uint8(B), [], uint8(cat_vol_morph(B, 'l')));
    addcase(ctx, ['S_' nm '_dilate'], 'dilate', d, [1.5 1 1 1], uint8(B), [], uint8(cat_vol_morph(B, 'd', 1.5)));
    addcase(ctx, ['S_' nm '_open'], 'open', d, [1 1 1 1], uint8(B), [], uint8(cat_vol_morph(B, 'o', 1)));
    addcase(ctx, ['S_' nm '_ldc'], 'ldc', d, [4.8 1 1 1], uint8(B), [], uint8(cat_vol_morph(B, 'ldc', 4.8) > 0.5));
    addcase(ctx, ['S_' nm '_ldc8'], 'ldc', d, [8 1 1 1], uint8(B), [], uint8(cat_vol_morph(B, 'ldc', 8) > 0.5));
    addcase(ctx, ['S_' nm '_ldo'], 'ldo', d, [3 1.5 1.5 1.5], uint8(B), [], uint8(cat_vol_morph(B, 'ldo', 3, 1.5) > 0.5));
  end

  % smoothing of small volumes and special filter sizes
  for sz = {[5 7 9], [10 9 8], [13 14 15], [7 7 7], [6 20 20]}
    d = sz{1};
    X = single(randn(d));
    for s = [0 -1 0.3 0.75 1 1.6 2 4 4.8 8]
      addcase(ctx, sprintf('S_smooth_%dx%dx%d_%g', d, s), 'smooth3X', d, s, X, [], cat_vol_smooth3X(X, s));
    end
  end
  d = G(1, :);
  X = single(randn(d)); X(rand(d) > 0.999) = NaN; X(rand(d) > 0.999) = Inf; X(rand(d) > 0.999) = -Inf;
  for s = [0.3 1 2.4]
    addcase(ctx, sprintf('S_smooth_nan_%g', s), 'smooth3X', d, s, X, [], cat_vol_smooth3X(X, s));
  end

  fclose(fid);
  fprintf('Reference data written to %s\n', outdir);
end

% ======================================================================
function addcase(ctx, name, func, d, p, in, mask, out, out2)
%ADDCASE Write input, mask and output(s) of one test case and its manifest line.
  if isempty(out) && strcmp(func, 'smooth3X')
    out = cat_vol_smooth3X(in, p(1));
  end
  p(end+1:4) = 0;
  fin = [name '_in.raw'];
  writeraw(fullfile(ctx.dir, fin), in);
  if isempty(mask)
    fmask = '-';
  else
    fmask = [name '_mask.raw'];
    writeraw(fullfile(ctx.dir, fmask), uint8(mask));
  end
  fout = [name '_out.raw'];
  writeraw(fullfile(ctx.dir, fout), out);
  line = sprintf('%s %s %d %d %d %.17g %.17g %.17g %.17g %s %s %s', name, func, d, p(1:4), fin, fmask, fout);
  if nargin > 8
    fout2 = [name '_out2.raw'];
    writeraw(fullfile(ctx.dir, fout2), out2);
    line = [line ' ' fout2];
  end
  fprintf(ctx.fid, '%s\n', line);
end

function writeraw(fn, X)
%WRITERAW Raw dump in MATLAB order: float32 for floating point, else the class.
  if islogical(X), X = uint8(X); end
  if isfloat(X), prec = 'float32'; else, prec = class(X); end
  f = fopen(fn, 'w');
  fwrite(f, X, prec);
  fclose(f);
end

function Y = smoothfield(d, fwhm)
%SMOOTHFIELD Gaussian noise smoothed with spm_smooth, scaled to max |Y| = 1.
  Y = randn(d);
  spm_smooth(Y, Y, [fwhm fwhm fwhm]);
  Y = single(Y / max(abs(Y(:))));
end

function [f1, f2] = brain(TGM, TWM, d, variant)
%BRAIN Template GM/WM resampled to d; 2: enlarged ventricle; 3: plus a tunnel to the outside.
  [X, Y, Z] = ndgrid(linspace(1, size(TGM, 1), d(1)), linspace(1, size(TGM, 2), d(2)), ...
    linspace(1, size(TGM, 3), d(3)));
  f1 = single(interpn(TGM, X, Y, Z, 'linear'));
  f2 = single(interpn(TWM, X, Y, Z, 'linear'));
  if variant > 1
    [x, y, z] = ndgrid(1:d(1), 1:d(2), 1:d(3));
    c = d .* [0.5 0.53 0.55];
    r = d .* [0.15 0.2 0.12];
    V = ((x - c(1)) / r(1)).^2 + ((y - c(2)) / r(2)).^2 + ((z - c(3)) / r(3)).^2 <= 1;
    if variant > 2
      V = V | (abs(x - c(1)) <= 1 & abs(z - c(3)) <= 1 & y >= c(2));
    end
    f1(V) = 0;
    f2(V) = 0;
  end
end

function B = blobs(d)
%BLOBS Ellipsoid with holes and a cavity, small separate objects, objects at the border.
  [x, y, z] = ndgrid(1:d(1), 1:d(2), 1:d(3));
  c = d / 2;
  B = ((x - c(1)) / (0.3 * d(1))).^2 + ((y - c(2)) / (0.3 * d(2))).^2 + ((z - c(3)) / (0.3 * d(3))).^2 <= 1;
  for k = 1:12
    p = c + (rand(1, 3) - 0.5) .* d * 0.4;
    B(((x - p(1)).^2 + (y - p(2)).^2 + (z - p(3)).^2) <= (1 + 3 * rand)^2) = false;
  end
  B(((x - c(1)).^2 + (y - c(2)).^2 + (z - c(3)).^2) <= 36) = false;                     % cavity
  B(((x - c(1)).^2 + (y - c(2)).^2 + (z - c(3)).^2) <= 4) = true;                      % island in it
  B(x <= 4 & y > d(2) * 0.2 & y < d(2) * 0.5 & z > d(3) * 0.3 & z < d(3) * 0.6) = true;  % slab at x = 1
  B(((x - d(1)).^2 + (y - d(2)).^2 + (z - 1).^2) <= 64) = true;                         % corner ball
  for k = 1:25
    p = 1 + rand(1, 3) .* (d - 1);
    B(((x - p(1)).^2 + (y - p(2)).^2 + (z - p(3)).^2) <= (0.5 + 2 * rand)^2) = true;
  end
  B(5:7, 5:7, 5:7) = true; B(5:7, 5:7, 9:11) = true;                                   % equal pair
end
