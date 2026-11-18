%% Peaucellier-Lipkin linkage
% O = fixed pivot at origin
% Q = fixed pivot at (r,0), QA = crank,
% OB/OD = long links, ABPD = rhombus, P = end effector.
% P moves on the straight line x = k/(2r), with k = L^2 - s^2.

clear; clc; close all;

%% Link lengths
r = 60;     % crank |QA| = |OQ|
L = 150;    % long links |OB| = |OD|
s = 60;    % rhombus sides |AB| = |BP| = |PD| = |DA|

k         = L^2 - s^2;                    % |OA|*|OP| = k
x_line    = k / (2*r);                    % P always has this x
theta_max = 2*acosd((L - s) / (2*r));     % crank angle limit [deg]
y_max     = x_line * tand(theta_max/2);   % P's travel limit along the line



%% Print Link Lengths

fprintf('  r (|QA| = |OQ|)= %f\n', r);
fprintf('  L (|OB| = |OD|)= %f\n', L);
fprintf('  s (|AB| = |BP| = |PD| = |DA|) = %f\n', s);


%% Ask for forward or reverse kinematics
mech = upper(strtrim(input('Forward or Reverse Kinematics(F/R): ', 's')));

if mech == 'F'
    theta = input(sprintf('Enter theta in degrees (-%.2f to %.2f): ', theta_max, theta_max));
    P = [x_line, x_line*tand(theta/2)];
elseif mech == 'R'
    fprintf('P moves on the line x = %.2f\n', x_line);
    y = input(sprintf('Enter output position y_P (-%.2f to %.2f): ', y_max, y_max));
    P = [x_line, y];
    theta = 2*atand(2*r*y / k);
else
    error('Enter F or R.');
end

if abs(theta) > theta_max
    error('Out of reach: theta = %.2f deg, limit is +/- %.2f deg.', theta, theta_max);
end


%% Joint positions
O = [0 0];
Q = [r 0];
A = Q + r*[cosd(theta), sind(theta)];
a = norm(A);                         % |OA|
u = A / a;                           % unit vector along OA
n = [-u(2), u(1)];                   % perpendicular to OA
OM = (a^2 + L^2 - s^2) / (2*a);      % O to midpoint of BD
h  = sqrt(max(L^2 - OM^2, 0));       % half of BD
B = OM*u + h*n;
D = OM*u - h*n;

if(sqrt((P(1)-Q(1))^2+P(2)^2)<r)
    mode = 'inverted';
else
    mode = 'regular';
end

fprintf('\nAssembly mode of the mechanism is %s\n', mode);

%% Print results
fprintf('\ntheta = %.4f deg\n', theta);
fprintf('P     = (%.4f, %.4f)\n\n', P);
fprintf('Joint positions:\n');
fprintf('  O = (%9.4f, %9.4f)\n', O);
fprintf('  Q = (%9.4f, %9.4f)\n', Q);
fprintf('  A = (%9.4f, %9.4f)\n', A);
fprintf('  B = (%9.4f, %9.4f)\n', B);
fprintf('  D = (%9.4f, %9.4f)\n', D);
fprintf('  P = (%9.4f, %9.4f)\n', P);

%% Two-solution check (long links OB and OD driven as the inputs)
% With the crank as input, theta is unique. If the two long links are
% driven instead (alpha1, alpha2 measured from the x-axis), the circle of
% radius L about O and the circle of radius s about A cross at two points,
% so the same P is reached by two different sets of link angles:
%   alpha = phi +/- beta,  phi = atan2(yP, xP),
%   cos(beta) = (|OA|^2 + L^2 - s^2) / (2 L |OA|)
% Two solutions when  L - s < |OP| < L + s,  one (singular) on the boundary.
dOP  = norm(P);
tol  = 1e-6 * L;
phi  = atan2d(P(2), P(1));
cosb = (a^2 + L^2 - s^2) / (2*L*a);

fprintf('\nInverse kinematics with the long links as inputs:\n');
fprintf('  |OP| = %.4f, allowed range %.4f to %.4f\n', dOP, L - s, L + s);

if dOP > (L - s) + tol && dOP < (L + s) - tol
    beta = acosd(cosb);
    alphas = [phi + beta, phi - beta;     % solution 1: [alpha1 alpha2]
              phi - beta, phi + beta];    % solution 2: links swapped
    fprintf('  TWO SOLUTIONS (phi = %.4f deg, beta = %.4f deg)\n', phi, beta);
    for i = 1:2
        J1 = L*[cosd(alphas(i,1)), sind(alphas(i,1))];   % end of long link 1
        J2 = L*[cosd(alphas(i,2)), sind(alphas(i,2))];   % end of long link 2
        Pchk = J1 + J2 - A;                              % rhombus closes at P
        fprintf('\n  Solution %d:\n', i);
        fprintf('    alpha1 (link 1) = %9.4f deg  ->  joint at (%9.4f, %9.4f)\n', alphas(i,1), J1);
        fprintf('    alpha2 (link 2) = %9.4f deg  ->  joint at (%9.4f, %9.4f)\n', alphas(i,2), J2);
        fprintf('    P check         = (%9.4f, %9.4f)\n', Pchk);
    end
elseif abs(dOP - (L - s)) <= tol || abs(dOP - (L + s)) <= tol
    fprintf('  ONE SOLUTION (singular): both long links lie along OP.\n');
    fprintf('    alpha1 = alpha2 = %.4f deg\n', phi);
else
    fprintf('  NO SOLUTION: P is outside the reachable ring.\n');
end

%% Plot the mechanism
figure('Color', 'w'); hold on; axis equal; grid on;
plot([Q(1) A(1)], [Q(2) A(2)], 'Color', 'B', 'LineWidth', 4);   % crank
plot([O(1) B(1)], [O(2) B(2)], 'Color', 'K', 'LineWidth', 4);   % long links
plot([O(1) D(1)], [O(2) D(2)], 'Color', 'K', 'LineWidth', 4);
plot([A(1) B(1) P(1) D(1) A(1)], [A(2) B(2) P(2) D(2) A(2)], ...
    'Color', 'G', 'LineWidth', 4);                             % rhombus
plot([x_line x_line], [-y_max y_max], 'r--');                                 % P's line

pts   = [O; Q; A; B; D; P];
names = {'O', 'Q', 'A', 'B', 'D', 'P'};
plot(pts(:,1), pts(:,2), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
plot(P(1), P(2), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
text(pts(:,1) + 6, pts(:,2) + 10, names, 'FontWeight', 'bold');

xlabel('x'); ylabel('y');
title(sprintf('\\theta = %.2f^\\circ,  P = (%.2f, %.2f)', theta, P(1), P(2)));


%% Joint Angle vs Output Position

flThetaMax = floor(theta_max);

AngleArr = [-theta_max, (flThetaMax*-1:1:flThetaMax), theta_max];   
CorrespondingQ = x_line*tand(AngleArr/2);

% Singularities
singTheta = theta_max*[-1 1];                 % Stretched
singName  = {'Edge', 'Edge'};
if k < 4*r^2                                  % Folded over itself
    thetaFold = 2*acosd(sqrt(k)/(2*r));
    singTheta = [singTheta, thetaFold*[-1 1]];
    singName  = [singName, {'Fold', 'Fold'}];
end
singY = x_line*tand(singTheta/2);

fprintf('\nSingularities:\n');
for i = 1:numel(singTheta)
    fprintf('  %-6s theta = %9.4f deg,  P_y = %9.4f\n', singName{i}, singTheta(i), singY(i));
end

figure();
hold on;
plot(AngleArr, CorrespondingQ, 'LineWidth', 2);
plot(theta, P(2), 'r*', 'MarkerSize', 10);
plot(singTheta, singY, 'g*', 'MarkerSize', 12, 'LineWidth', 1.5);
text(singTheta + 3, singY, singName, 'Color', [0 0.5 0]);
xlabel("Theta (degrees)")
ylabel("Corresponding output y (Distance Units)")
legend("Output Curve", "Inputted Point (theta, P_y)", "Singularities", 'Location', 'best')
hold off;




