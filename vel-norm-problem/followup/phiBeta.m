load('../scripts/variables.mat');
source('../scripts/functions.m');

% Defining yp as anonymous function
y_p = @(Beta_p) sqrt(Beta_p^4*(3*pi + 8)^2 * (sin(2*phi + pi/2)).^2 ...
    - (Beta_p^4*(3*pi + 8)^2 - 8*(3*phi + 2)*(3*pi + 8)*Beta_p^2 + 16*(3*phi + 2).^2 + 2*(4*(sin(phi + pi/4)).^2 + 6) ...
    .* ((3*pi + 8)*Beta_p^2 - 4*(3*phi + 2)) .* sin(2*phi + pi/2) + (16*(sin(phi+pi/4)).^4 ...
    + 48*(sin(phi+pi/4)).^2 + 36) .* (sin(2*phi+pi/2)).^2)) ./ (Beta_p^2 * (3*pi+8) * sin(2*phi+pi/2)) * sqrt(mu/r_0) .* cot(phi + pi/4);

% return real part of a vector
function y = realBreakpoint(x)
    y = 0;
    for i = x,
        if ~isreal(i)
            break;
        end
        y = y + 1;
    end
end

% find 2 smallest values (y(1) = 0)
[y_sort idx_vec] = sort(y_p(Beta)(1:realBreakpoint(y_p(Beta))));

% Minimum under some given beta
phi_min = phi(idx_vec(2));
y_min = y_sort(2);

% find when y_p = 0 --> phi_m and redefine delta and beta as vectors
delta = (1:1:n)./1000; % n = 1000
Beta_vec = sqrt(1./delta);
phi_m1 = zeros(1,n);
w_min = zeros(1,n);

% define w_p
w_p = @(Beta_p) 1./cos(2*phi).*(1-(4/(Beta_p^2*(3*pi+8)))*(3*phi+2))+2/(Beta_p^2*(3*pi+8))*(2*sin(phi+pi/4).^2+3);

% populate phi_m by iterating through beta_vec
for i = 1:n
    y_real = y_p(Beta_vec(i))(1:realBreakpoint(y_p(Beta_vec(i))));
    [null idx] = sort(y_real);
    if (length(idx) <= 1)
        phi_m(i) = phi_m1(max(1,i-1));
    else
        w_min(i) = w_p(Beta_vec(i))(idx(2));
        phi_m1(i) = phi(idx(2));
    endif
end

% plot delta over phi_m
figure(1);
hold on;
plot(delta,phi_m1,'b', 'linewidth', 3);
xlabel("{\\it \\delta}");
ylabel("{\\it \\phi_m}");
title("{\\it \\phi_m} vs. {\\it \\delta} s.t. y_p(\\phi) = 0");
% legend({"{\\it \\phi_m}", "{\\it w(\\phi_m)}"}, 'location', 'northwest', 'orientation', 'vertical');
% legend boxoff;
exportPlot('phi_delta', 1);

% minimize w_p for different values of delta
phi_m2 = zeros(1,n);
for i = 1:n
    [null idx] = sort(w_p(Beta_vec(i)));
    phi_m2(i) = phi(idx(1));
end

figure(2);
hold on;
plot(delta,phi_m2,'r', 'linewidth', 3);
xlabel("{\\it \\delta}");
ylabel("{\\it \\phi_m}");
title("{\\it \\phi_m} vs. {\\it \\delta} for min(w_p(\\phi))");
% legend({"{\\it \\phi_m}", "{\\it w(\\phi_m)}"}, 'location', 'northwest', 'orientation', 'vertical');
% legend boxoff;
exportPlot('phi_delta_min', 2);
close all;