load('../scripts/variables.mat');
source('../scripts/functions.m');

% Defining yp as anonymous function
g_p = @(Beta_p) Beta_p^4*(3*pi + 8)^2 * (sin(2*phi + pi/2)).^2 ...
    - (Beta_p^4*(3*pi + 8)^2 - 8*(3*phi + 2)*(3*pi + 8)*Beta_p^2 + 16*(3*phi + 2).^2 + 2*(4*(sin(phi + pi/4)).^2 + 6) ...
    .* ((3*pi + 8)*Beta_p^2 - 4*(3*phi + 2)) .* sin(2*phi + pi/2) + (16*(sin(phi+pi/4)).^4 ...
    + 48*(sin(phi+pi/4)).^2 + 36) .* (sin(2*phi+pi/2)).^2);
y_p = @(Beta_p) sqrt([0 g_p(Beta_p)(2:end)]) ./ (Beta_p^2 * (3*pi+8) * sin(2*phi+pi/2)) * sqrt(mu/r_0) .* cot(phi + pi/4);

% return real part of a vector
function y = realBreakpoint(x)
    y = 0;
    for i = x,
        if iscomplex(i)
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
w_1 = zeros(1,n);
y_0 = zeros(1,n);

% define w_p
w_p = @(Beta_p) 1./cos(2*phi).*(1-(4/(Beta_p^2*(3*pi+8)))*(3*phi+2))+2/(Beta_p^2*(3*pi+8))*(2*sin(phi+pi/4).^2+3);

% populate phi_m by iterating through beta_vec
for i = 1:n
    y_real = y_p(Beta_vec(i))(1:realBreakpoint(y_p(Beta_vec(i))));
    [yvals idx] = sort(y_real);
    if (length(idx) > 1)
        w_1(i) = w_p(Beta_vec(i))(idx(2));
        phi_m1(i) = phi(idx(2));
        y_0(i) = yvals(2);
    endif
end

% plot delta over phi_m1
if (1)
figure(1);
hold on;
plot(delta,phi_m1,'b', 'linewidth', 3);
plot(delta,w_1,'k', 'linewidth', 1);
plot(delta,y_0,'r', 'linewidth', 1);
ylim([0 1.5]);
xlabel("{\\it \\delta}");
ylabel("{\\it \\phi_m}");
title("{\\it \\phi_m} vs. {\\it \\delta} s.t. {\\it y_p(\\phi_m) = 0}");
legend({"{\\it \\phi_m}", "{\\it w_p(\\phi_m)}", "{\\it y_p(\\phi_m)}"}, 'location', 'northwest', 'orientation', 'vertical');
legend boxoff;
exportPlot('phi_delta_y', 1);

% minimize w_p for different values of delta
phi_m2 = zeros(1,n);
w_min = zeros(1,n);
for i = 1:n
    [w_vals idx] = sort(w_p(Beta_vec(i)));
    phi_m2(i) = phi(idx(1));
    w_min(i) = w_vals(1);
end

figure(2);
hold on;
plot(delta,phi_m2,'b', 'linewidth', 3);
plot(delta,w_min,'r', 'linewidth', 1);
ylim([0 1.5]);
xlabel("{\\it \\delta}");
ylabel("{\\it \\phi_m}");
title("{\\it \\phi_m} vs. {\\it \\delta} for min({\\itw_p(\\phi)})");
legend({"{\\it \\phi_m}", "{\\it w_p(\\phi_m)}"}, 'location', 'northwest', 'orientation', 'vertical');
legend boxoff;
exportPlot('phi_delta_min', 2);

close all;
endif