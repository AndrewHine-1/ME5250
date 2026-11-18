%% ME5250 Project 1 - Peaucellier-Lipkin: mechanical advantage vs input (D6)
%   q_a = theta, q_p = [phi2..phi7] = absolute angles of
%         O->A, O->C, A->B, B->C, C->D, D->A
%   J  = df/dqa - df/dqp * inv(Hp) * Ha   (MR 7.3),  tau = J' F
%   Mechanical advantage for a load F_y at D:  F_y / tau = 1 / J_y

clear; clc; close all;

r = 60;  L = 150;  s = 60;                 % mm
thMax = 2*acos((L - s)/(2*r));             % toggle OB = L - s  [rad]

% Sweep the open range; Hp is singular at the toggle ends (+-thMax)
th = linspace(-thMax, thMax, 501);
th = th(2:end-1);
invJ = zeros(size(th));
for i = 1:numel(th)
    q = fk(th(i), r, L, s);
    [Ha, Hp, dfdqp] = closure_jacobians(th(i), q, r, L, s);
    J = -dfdqp*(Hp\Ha);                    % df/dqa = 0
    invJ(i) = 1/J(2);
end

% Configuration-space singularities at theta = +-thMax (J_y finite there)
k = L^2 - s^2;
invJsing = 1/(k/(4*r)*sec(thMax/2)^2);

figure('Color', 'w'); hold on; grid on;
plot(rad2deg(th), invJ, 'LineWidth', 2);
plot(rad2deg(thMax)*[-1 1], invJsing*[1 1], '*', 'MarkerSize', 14);
xline(rad2deg(thMax)*[-1 1], 'k--');
xlabel('\theta (deg)'); ylabel('F_y / \tau = 1/J_y  (mm^{-1})');
legend('Mechanical advantage', 'C-space singularity', 'Location', 'best');
title('Mechanical advantage vs input');

%% ================= Local functions =================
function q = fk(th, r, L, s)
% Closed-form FK on the branch used in the report: A above OB, rhombus open.
    B = [r 0] + r*[cos(th) sin(th)];
    b = norm(B);  u = B/b;  nrm = [-u(2) u(1)];
    a = (L^2 - s^2 + b^2)/(2*b);
    hh = sqrt(max(L^2 - a^2, 0));
    A = a*u + hh*nrm;  C = a*u - hh*nrm;  D = A + C - B;
    ang = @(v) atan2(v(2), v(1));
    q = [ang(A), ang(C), ang(B - A), ang(C - B), ang(D - C), ang(A - D)];
end

function [Ha, Hp, dfdqp] = closure_jacobians(th, q, r, L, s)
% Ha = dg/dtheta (6x1), Hp = dg/dq_p (6x6), df/dq_p (2x6).
    ep = @(t) [-sin(t); cos(t)];
    z = zeros(2, 1);
    Ha = [r*ep(th); r*ep(th); z];
    Hp = [-L*ep(q(1)),  z,          -s*ep(q(3)),  z,           z,          z;
           z,          -L*ep(q(2)),  z,           s*ep(q(4)),  z,          z;
           z,           z,           s*ep(q(3)),  s*ep(q(4)),  s*ep(q(5)), s*ep(q(6))];
    dfdqp = [L*ep(q(1)), z, z, z, z, -s*ep(q(6))];
end
