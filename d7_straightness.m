%% Peaucellier-Lipkin linkage: D7 straightness (spec note 1)
% Exact mechanism: tracer path is straight to numerical precision
% over the full stroke.
% Lengthen ONE link by 1% and report e/S over the
% same input range.
%   S = extent of the path along its least-squares line
%   e = max perpendicular deviation from that line
%
% Naming as in peaucellier_fk_ik.m:
%   O = arm pivot (origin), E = crank pivot (r,0), B = inner rhombus vertex,
%   A / C = outer vertices, D = tracer.
% Links (8 independent lengths so any one can be perturbed):
%   ell(1) [0] ground OE   ell(2) [1] crank EB   ell(3) [2] arm OA   ell(4) [3] arm OC
%   ell(5) [4] AB          ell(6) [5] BC         ell(7) [6] CD       ell(8) [7] DA
% Unknowns (absolute link angles, rad, along the listed direction):
%   qa = theta (E->B);  qp = [phi2 phi3 phi4 phi5 phi6 phi7] for O->A, O->C,
%   A->B, B->C, C->D, D->A.
clear; clc; close all;

%% 1. Link lengths (mm)
r = 60;  L = 150;  s = 60;
ELL   = [r r L L s s s s];
names = {'[0] ground OE', '[1] crank EB', '[2] arm OA', '[3] arm OC', ...
         '[4] AB', '[5] BC', '[6] CD', '[7] DA'};
k      = L^2 - s^2;
x_line = k / (2*r);

%% 2. Input ranges
% Nominal full stroke ends at the toggle |OB| = L - s.
theta_full = 2*acosd((L - s) / (2*r));                 % 82.82 deg
% A 1%-long arm moves that toggle in to |OB| = 1.01L - s, so use the largest
% round range that every perturbed linkage can still assemble over.
theta_arm  = 2*acosd((1.01*L - s) / (2*r));            % 80.6 deg
theta_com  = floor(theta_arm);                         % 80 deg
fprintf('Full stroke |theta| <= %.2f deg; 1%%-arm toggle at %.2f deg -> common range |theta| <= %d deg\n\n', ...
        theta_full, theta_arm, theta_com);

th_full = deg2rad(-82.8:0.1:82.8);
th_com  = deg2rad(-theta_com:0.1:theta_com);

%% 3. Check Hp against finite differences (trust nothing unchecked)
th_t = deg2rad(50);
qp_t = seed_qp(th_t, r, L, s);
hstep = 1e-6;  Hfd = zeros(6);
for j = 1:6
    dq = zeros(1,6);  dq(j) = hstep;
    Hfd(:,j) = (loop_g(th_t, qp_t+dq, ELL) - loop_g(th_t, qp_t-dq, ELL)) / (2*hstep);
end
fprintf('Check: |g| at seed = %.1e,  max|Hp - FD| = %.1e\n\n', ...
        norm(loop_g(th_t, qp_t, ELL)), max(abs(loop_Hp(qp_t, ELL) - Hfd), [], 'all'));

%% 4. Nominal linkage
[~, D_full] = fk_sweep(th_full, ELL, r, L, s);
[S_f, e_f]  = line_fit(D_full);
[~, D_nom]  = fk_sweep(th_com, ELL, r, L, s);
[S_n, e_n, dev_nom] = line_fit(D_nom);
fprintf('NOMINAL  (D_x range %.12f to %.12f, predicted %.1f)\n', min(D_full(:,1)), max(D_full(:,1)), x_line);
fprintf('  %-22s S = %8.3f mm   e = %.2e mm   e/S = %.2e\n', 'full stroke +/-82.8', S_f, e_f, e_f/S_f);
fprintf('  %-22s S = %8.3f mm   e = %.2e mm   e/S = %.2e\n\n', sprintf('common +/-%d', theta_com), S_n, e_n, e_n/S_n);

%% 5. One link 1% long, each in turn, over the common range
fprintf('ONE LINK 1%% LONG, |theta| <= %d deg\n', theta_com);
fprintf('  %-15s %9s %9s %10s\n', 'link', 'S (mm)', 'e (mm)', 'e/S');
devs = zeros(numel(th_com), 8);  eS = zeros(1, 8);
for i = 1:8
    ell    = ELL;
    ell(i) = 1.01 * ell(i);
    [~, D] = fk_sweep(th_com, ell, r, L, s);
    [S, e, dev] = line_fit(D);
    devs(:,i) = dev;  eS(i) = e/S;
    fprintf('  %-15s %9.3f %9.4f %10.2e\n', names{i}, S, e, e/S);
end

%% 6. Figure: deviation from the fitted line vs input
figure('Color', 'w', 'Position', [100 100 950 360]);
subplot(1,2,1);
plot(rad2deg(th_com), dev_nom*1e12, 'k'); grid on;
xlabel('\theta (deg)'); ylabel('deviation (\times10^{-12} mm)');
title('Nominal: deviation from fitted line');
subplot(1,2,2); hold on; grid on;
styles = {'-', '-', '--', '--', ':', ':', '-.', '-.'};
for i = 1:8
    plot(rad2deg(th_com), devs(:,i), styles{i}, 'LineWidth', 1.2, ...
         'DisplayName', sprintf('%s: e/S = %.1e', names{i}, eS(i)));
end
xlabel('\theta (deg)'); ylabel('deviation from fitted line (mm)');
title('One link 1% long'); legend('Location', 'southwest', 'FontSize', 7);
exportgraphics(gcf, 'fig_d7_straightness.png', 'Resolution', 200);
fprintf('\nSaved fig_d7_straightness.png\n');

%% ===================== Local functions =====================

function v = ev(a)        % unit vector at angle a: e(a) = [cos a; sin a]
    v = [cos(a); sin(a)];
end

function v = dev_(a)      % de/da = [-sin a; cos a]  (e rotated 90 deg)
    v = [-sin(a); cos(a)];
end

function G = loop_g(th, qp, ell)
% Loop closure g(q) = 0 (6x1), one vector equation per independent loop.
    OE = [ell(1); 0];
    g1 = OE + ell(2)*ev(th) - ell(5)*ev(qp(3)) - ell(3)*ev(qp(1));  % O-E-B-A
    g2 = OE + ell(2)*ev(th) + ell(6)*ev(qp(4)) - ell(4)*ev(qp(2));  % O-E-B-C
    g3 = ell(5)*ev(qp(3)) + ell(6)*ev(qp(4)) ...
       + ell(7)*ev(qp(5)) + ell(8)*ev(qp(6));                       % A-B-C-D
    G = [g1; g2; g3];
end

function H = loop_Hp(qp, ell)
% Hp = dg/dqp (6x6). Columns phi2..phi7, rows [loop1; loop2; loop3].
    H = zeros(6);
    H(1:2,1) = -ell(3)*dev_(qp(1));      % loop 1, phi2 (arm OA)
    H(1:2,3) = -ell(5)*dev_(qp(3));      % loop 1, phi4 (AB)
    H(3:4,2) = -ell(4)*dev_(qp(2));      % loop 2, phi3 (arm OC)
    H(3:4,4) =  ell(6)*dev_(qp(4));      % loop 2, phi5 (BC)
    for j = 3:6                          % loop 3, phi4..phi7 (rhombus)
        H(5:6,j) = ell(j+2)*dev_(qp(j));
    end
end

function Dp = tracer(qp, ell)
% Output: D = OA + AD = l2 e(phi2) - l7 e(phi7)   (1x2 row)
    Dp = (ell(3)*ev(qp(1)) - ell(8)*ev(qp(6)))';
end

function qp = seed_qp(th, r, L, s)
% Nominal geometric construction, used only as Newton's first guess.
% Picks the assembly mode: A on +n side of OB, C on -n side, rhombus open.
    B = [r 0] + r*[cos(th) sin(th)];
    a = norm(B);  u = B/a;  n = [-u(2) u(1)];
    OM = (a^2 + L^2 - s^2) / (2*a);
    h  = sqrt(max(L^2 - OM^2, 0));
    A = OM*u + h*n;  C = OM*u - h*n;  D = A + C - B;  O = [0 0];
    ang = @(P,Q) atan2(Q(2)-P(2), Q(1)-P(1));
    qp = [ang(O,A) ang(O,C) ang(A,B) ang(B,C) ang(C,D) ang(D,A)];
end

function qp = newton_qp(th, qp, ell)
% Newton-Raphson: repeat qp <- qp - Hp \ g until |g| < 1e-12 mm.
% The update happens before the test, so one extra step polishes to round-off.
    for it = 1:50
        G  = loop_g(th, qp, ell);
        qp = qp - (loop_Hp(qp, ell) \ G)';
        if norm(G) < 1e-12, return; end
    end
    G = loop_g(th, qp, ell);
    if norm(G) > 1e-8
        error('Newton failed at theta = %.2f deg (|g| = %.1e): not assemblable?', rad2deg(th), norm(G));
    end
end

function [QP, D] = fk_sweep(thetas, ell, r, L, s)
% FK over a theta grid containing 0. Start at theta = 0 from the geometric
% seed, step outward, and seed each solve with its neighbour so the solver
% stays on one assembly mode.
    n  = numel(thetas);
    [~, i0] = min(abs(thetas));
    QP = zeros(n, 6);
    QP(i0,:) = newton_qp(thetas(i0), seed_qp(thetas(i0), r, L, s), ell);
    for i = i0+1:n                      % 0 -> +max
        QP(i,:) = newton_qp(thetas(i), QP(i-1,:), ell);
    end
    for i = i0-1:-1:1                   % 0 -> -max
        QP(i,:) = newton_qp(thetas(i), QP(i+1,:), ell);
    end
    D = zeros(n, 2);
    for i = 1:n
        D(i,:) = tracer(QP(i,:), ell);
    end
end

function [S, e, dev] = line_fit(D)
% Total least-squares line through the points D (n x 2):
%   line passes through the centroid c, direction t = first right singular
%   vector of (D - c), normal nrm = second. Minimises perpendicular distances.
    c = mean(D, 1);
    M = D - c;
    [~, ~, V] = svd(M, 'econ');
    t = V(:,1);  nrm = V(:,2);
    along = M * t;
    dev   = M * nrm;
    S = max(along) - min(along);
    e = max(abs(dev));
end
