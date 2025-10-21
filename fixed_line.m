function dem_main()
% DEM main — same physics as your scripts:
% - bonds: linear spring (k_b) + dashpot (c_b) along the bond axis
% - failure: breaks when tensile strain > epsCrit
% - integration: semi-implicit Euler
% - optional air drag (c_air)
%
% Scenario 1 (default): fixed-left line, pull right-most bead with a ramp force.
% Change the SCENARIO section to do other cases without touching the physics.

clc; clear; close all;

%% ===================== USER KNOBS (scenario) =====================
SCENARIO.name = 'line_pull_right';        % 'line_pull_right' (default)
SCENARIO.N     = 8;                       % beads in the line
SCENARIO.R     = 0.35;                    % radius (m), spacing = 2R
SCENARIO.y0    = 0.0;                     % vertical position of the line

% Boundary conditions (indices in [1..N])
SCENARIO.fixed = 1;                       % fix these beads (vector allowed)
SCENARIO.pull  = SCENARIO.N;              % bead that receives the actuator

% Actuator: ramped force Fx(t) to bead "pull"
SCENARIO.Fmax  = 1500;                    % N
SCENARIO.Tramp = 1.0;                     % s (linear ramp to Fmax)

% Gravity OFF for tensile test; set g>0 for other scenes
S.g       = 0.0;                          % m/s^2
S.T       = 4.0;                          % total time (s)
S.dt      = 0.001;                        % time step (s)
S.c_air   = 0.05;                         % light air damping (0=off)

% Bond model (exactly your previous values/meaning)
S.k_b     = 1.0e4;                        % bond stiffness (N/m)
S.c_b     = 80;                           % bond damping
S.epsCrit = 0.15;                         % failure strain (tension)
S.m       = 1.0;                          % mass per bead (kg)

% Plot window
domainX = [-1, (SCENARIO.N-1)*(2*SCENARIO.R)+4];
domainY = [-2, 2];

%% ===================== GEOMETRY & BONDS ==========================
N  = SCENARIO.N; R = SCENARIO.R; dx = 2*R;
xy0 = [ (0:N-1)'*dx , SCENARIO.y0*ones(N,1) ];
xy  = xy0;                      % positions
v   = zeros(N,2);               % velocities
f   = zeros(N,2);               % forces
C   = lines(N);                 % colors

% Neighbor bonds only (1-2, 2-3, ..., N-1--N)
pairs = [(1:N-1)' (2:N)'];
Nb    = size(pairs,1);
L0    = repmat(dx, Nb, 1);      % reference lengths
alive = true(Nb,1);

fixedIdx = SCENARIO.fixed(:).';
pullIdx  = SCENARIO.pull;

anchorPos = xy0(fixedIdx,:);    % exact anchors for fixed beads

%% ===================== HISTORY / FIGURE ==========================
steps   = round(S.T/S.dt);
H.bonds = zeros(steps,1);
H.Fapp  = zeros(steps,1);
H.ext   = zeros(steps,1);

figure('Color','k','Position',[80 80 1120 460]);
fg = [.85 .85 .85];

%% ======================= TIME LOOP ===============================
for s = 1:steps
  t = s*S.dt;

  % ---------- reset forces ----------
  f(:) = 0;
  if S.c_air > 0, f = f - S.c_air * v; end
  % gravity (kept for completeness; here g=0)
  f(:,2) = f(:,2) - S.m*S.g;

  % ---------- bonds: spring + dashpot + failure ----------
  for b = 1:Nb
    if ~alive(b), continue; end
    i = pairs(b,1); j = pairs(b,2);
    rij = xy(j,:) - xy(i,:);
    L   = norm(rij); if L < 1e-12, continue; end
    n   = rij / L;
    vij = v(j,:) - v(i,:);
    vrel = dot(vij,n);

    dL  = L - L0(b);
    eps = dL / L0(b);             % tensile strain
    Fb  = S.k_b*dL + S.c_b*vrel;  % axial force (+tension)

    f(i,:) = f(i,:) + Fb*n;
    f(j,:) = f(j,:) - Fb*n;

    if eps > S.epsCrit, alive(b) = false; end
  end

  % ---------- actuator: ramped force on pull bead ----------
  if t <= SCENARIO.Tramp
    Fapp = SCENARIO.Fmax * (t/SCENARIO.Tramp);
  else
    Fapp = SCENARIO.Fmax;
  end
  f(pullIdx,1) = f(pullIdx,1) + Fapp;   % ONLY the right bead is pulled

  % ---------- enforce fixed BC before integrate ----------
  v(fixedIdx,:)  = 0;
  f(fixedIdx,:)  = 0;                   % remove any residual force

  % ---------- integrate (semi-implicit Euler) ----------
  a  = f / S.m;
  v  = v + a*S.dt;
  xy = xy + v*S.dt;

  % ---------- enforce fixed BC after integrate (hard clamp) ----------
  xy(fixedIdx,:) = anchorPos;
  v(fixedIdx,:)  = 0;

  % ---------- history ----------
  H.bonds(s) = sum(alive);
  H.Fapp(s)  = Fapp;
  H.ext(s)   = xy(pullIdx,1) - xy0(pullIdx,1);

  % ===================== DRAW ======================
  clf;

  % left: chain view
  subplot(1,3,1); hold on; axis equal
  set(gca,'Color','k','XColor',fg,'YColor',fg);
  xlim(domainX); ylim(domainY);
  title(sprintf('Fixed-Left Line  |  t = %.2fs  |  ext = %.3f m', t, H.ext(s)), 'Color','w');

  % bonds (cyan alive / gray broken)
  for b = 1:Nb
    i = pairs(b,1); j = pairs(b,2);
    col = alive(b)*[0 .8 1] + (~alive(b))*[.4 .4 .4];
    plot([xy(i,1) xy(j,1)], [xy(i,2) xy(j,2)], '-', 'Color', col, 'LineWidth', 2);
  end
  % beads
  th = linspace(0,2*pi,40);
  for i = 1:N
    fill(xy(i,1)+R*cos(th), xy(i,2)+R*sin(th), C(i,:), 'EdgeColor','none');
  end
  % markers: fixed squares; pull arrow on the right bead
  plot(xy(fixedIdx,1),xy(fixedIdx,2),'ws','MarkerSize',9,'MarkerFaceColor','w');
  quiver(xy(pullIdx,1),xy(pullIdx,2), 0.9, 0, 'y','LineWidth',2,'MaxHeadSize',2);

  % middle: bonds alive
  subplot(1,3,2);
  plot((1:s)*S.dt, H.bonds(1:s), 'c-','LineWidth',2); grid on;
  set(gca,'Color','k','XColor',fg,'YColor',fg);
  xlabel('Time (s)'); ylabel('Bonds Alive'); title('Failure Progress','Color','w');

  % right: Force–Extension (applied force vs actuator extension)
  subplot(1,3,3);
  plot(H.ext(1:s), H.Fapp(1:s), 'y-','LineWidth',2); grid on;
  set(gca,'Color','k','XColor',fg,'YColor',fg);
  xlabel('Extension at right end (m)'); ylabel('Applied Force F (N)');
  title('Force–Extension','Color','w');

  drawnow;
end

disp('Done. Same DEM physics; scenario set by SCENARIO.* knobs.');
end
