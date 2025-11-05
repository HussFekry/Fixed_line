function DEM_plate_pull_3x3()
% 3x3 grid pulled by a right rigid plate; left plate fixed (matches sketch)
% Physics (unchanged):
%   - Bonds: linear spring (k_b) + dashpot (c_b) along bond axis
%   - Failure: breaks when tensile strain > epsCrit
%   - Integrator: semi-implicit Euler
%   - Optional global (air) damping

clc; close all;

%% --------- PHYSICS ---------
S.k_b     = 1.0e4;     % bond stiffness (N/m)
S.c_b     = 80;        % bond damping
S.epsCrit = 0.18;      % tensile failure strain
S.m       = 1.0;       % bead mass (kg)
S.R       = 0.35;      % bead radius (m)
S.g       = 0.0;       % gravity off
S.c_air   = 0.05;      % small global damping

S.dt      = 0.001;     % time step (s)
S.T       = 2.0;       % total time (s)

%% --------- GEOMETRY to match the sketch ---------
nx = 3; ny = 3;                 % 3 columns × 3 rows
dx = 2*S.R; dy = 2*S.R;         % spacing
ofs = [0.5, 0.5];               % little margin around

% node coordinates (column-major)
N  = nx*ny;
id = @(ix,iy) (iy-1)*nx + ix;

xy0 = zeros(N,2);
for iy = 1:ny
  for ix = 1:nx
    xy0(id(ix,iy),:) = [ (ix-1)*dx, (iy-1)*dy ] + ofs;
  end
end
xy = xy0;                      % positions
v  = zeros(N,2);
f  = zeros(N,2);
C  = lines(N);

% neighbor bonds: horizontal + vertical (no diagonals)
pairs = [];
L0    = [];
% horizontals
for iy = 1:ny
  for ix = 1:nx-1
    i = id(ix,iy); j = id(ix+1,iy);
    pairs = [pairs; i j]; %#ok<AGROW>
    L0    = [L0; dx];     %#ok<AGROW>
  end
end
% verticals
for ix = 1:nx
  for iy = 1:ny-1
    i = id(ix,iy); j = id(ix,iy+1);
    pairs = [pairs; i j]; %#ok<AGROW>
    L0    = [L0; dy];     %#ok<AGROW>
  end
end

Nb    = size(pairs,1);
alive = true(Nb,1);

% columns (for BCs)
leftCol  = arrayfun(@(iy) id(1,iy),      1:ny).';
rightCol = arrayfun(@(iy) id(nx,iy),     1:ny).';

% fixed left plate (wall) keeps left column at their reference
anchorPos_fixed = xy0(leftCol,:);

% right moving plate (ram) imposes the SAME x displacement to the whole right column
Umax  = 0.40;     % imposed displacement (m)
Tramp = 0.20;     % ramp time (s)

%% --------- HISTORY & FIGURES ---------
steps = round(S.T/S.dt);
H.bonds = zeros(steps,1);
H.reac  = zeros(steps,1);  % reaction at the fixed wall
H.uimp  = zeros(steps,1);  % imposed displacement of the ram

figure('Color','k','Position',[60 80 1120 470]);
fg = [.85 .85 .85];

% precompute wall/ram geometry for nicer drawing
x_leftWall  = min(xy0(:,1)) - 0.35;  % a bit left of left column
x_rightFace0= max(xy0(:,1)) + S.R;   % initial position (through bead centers)
plateW      = 0.15;                  % drawn thickness
ymin = min(xy0(:,2)) - 1.0;
ymax = max(xy0(:,2)) + 1.0;

%% ===================== TIME LOOP =====================
for s = 1:steps
  t = s*S.dt;

  % reset forces
  f(:) = 0;
  if S.c_air > 0, f = f - S.c_air*v; end
  f(:,2) = f(:,2) - S.m*S.g;

  % bonds: spring + dashpot + failure
  for b = 1:Nb
    if ~alive(b), continue; end
    i = pairs(b,1); j = pairs(b,2);
    rij = xy(j,:) - xy(i,:);
    L   = norm(rij); if L < 1e-12, continue; end
    n   = rij / L;
    vij = v(j,:) - v(i,:);
    vrel= dot(vij,n);

    dL  = L - L0(b);
    eps = dL / L0(b);
    Fb  = S.k_b*dL + S.c_b*vrel;

    f(i,:) = f(i,:) + Fb*n;
    f(j,:) = f(j,:) - Fb*n;

    if eps > S.epsCrit, alive(b) = false; end
  end

  % imposed displacement on the right plate (same u for all beads in rightCol)
  if t <= Tramp
    uimp = Umax * (t/Tramp);
  else
    uimp = Umax;
  end
  xy(rightCol,1) = xy0(rightCol,1) + uimp;
  v(rightCol,1)  = 0;        % keep them moving with the plate (no relative x)

  % reaction on the fixed wall = internal force that would act on the clamps
  Rx = sum(f(leftCol,1));
  H.reac(s) = Rx;
  H.uimp(s) = uimp;

  % hard clamp at the left plate
  xy(leftCol,:) = anchorPos_fixed;
  v(leftCol,:)  = 0;
  f(leftCol,:)  = 0;

  % integrate free nodes
  a  = f / S.m;
  v  = v + a*S.dt;
  xy = xy + v*S.dt;

  % history
  H.bonds(s) = nnz(alive);

  %% -------------------- DRAW --------------------
  clf;

  % left: world view (exactly like the sketch)
  subplot(1,3,1); hold on; axis equal
  set(gca,'Color','k','XColor',fg,'YColor',fg);
  xlim([x_leftWall-0.5, x_rightFace0+Umax+2]); ylim([ymin, ymax]);
  title(sprintf('3×3 | left wall fixed | right plate pulled  |  t=%.2fs  |  u=%.3f m', t, uimp),'Color','w');

  % draw left wall (fixed)
  patch([x_leftWall x_leftWall+plateW x_leftWall+plateW x_leftWall], ...
        [ymin ymax ymax ymin], [0.8 0.8 0.9], 'EdgeColor','none');

  % draw right moving ram (plate)
  x_rightFace = x_rightFace0 + uimp;
  patch([x_rightFace x_rightFace+plateW x_rightFace+plateW x_rightFace], ...
        [ymin ymax ymax ymin], [0.9 0.8 0.4], 'EdgeColor','none');

  % bonds (cyan alive, gray broken)
  for b = 1:Nb
    i = pairs(b,1); j = pairs(b,2);
    col = alive(b)*[0 .8 1] + (~alive(b))*[.4 .4 .4];
    plot([xy(i,1) xy(j,1)], [xy(i,2) xy(j,2)], '-', 'Color', col, 'LineWidth', 2.5);
  end

  % beads
  th = linspace(0,2*pi,40);
  for i = 1:N
    fill(xy(i,1)+S.R*cos(th), xy(i,2)+S.R*sin(th), C(i,:), 'EdgeColor','none');
  end

  % markers: show which nodes are constrained / imposed
  plot(xy(leftCol,1), xy(leftCol,2), 'ws','MarkerFaceColor','w','MarkerSize',9);          % fixed
  quiver(xy(rightCol,1)-0.6, xy(rightCol,2), 0.6, 0, 'y','LineWidth',2,'MaxHeadSize',2); % plate push arrows

  % middle: bonds alive (should start at 12 for a 3×3 orthogonal net)
  subplot(1,3,2);
  plot((1:s)*S.dt, H.bonds(1:s), 'c-','LineWidth',2); grid on;
  set(gca,'Color','k','XColor',fg,'YColor',fg);
  xlabel('Time (s)'); ylabel('Bonds Alive'); title('Failure Progress','Color','w');

  % right: reaction (fixed wall) vs imposed displacement (ram)
  subplot(1,3,3);
  plot(H.uimp(1:s), H.reac(1:s), 'y-','LineWidth',2); grid on;
  set(gca,'Color','k','XColor',fg,'YColor',fg);
  xlabel('Imposed displacement u (m) at right plate');
  ylabel('Reaction force at fixed wall (N)');
  title('Reaction–Displacement','Color','w');

  drawnow;
end

disp('Done: 3×3 grid with left fixed plate and right moving plate (displacement control).');
end
