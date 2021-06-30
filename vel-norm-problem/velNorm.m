% -------------------- DEFINITIONS --------------------
% plot variable
n = 1000;
inc = (pi/4)/n;
phi = (0:1:n-1) * inc;

% constants
r_0 = 1;
a_T = 0.2;
mu = 1;

width = 3;

B = sqrt(4 * mu/((3*pi + 8) * a_T * r_0^2));

% definition of r
r = 2 * r_0 .* sin(phi + pi/4) .^ 2;

% plot export function
function exportPlot(fileName, fig)
    print(fig, sprintf('%s.pdf', fileName), '-dpdfcairo');
    system(sprintf('mv %s.pdf plots/%s.pdf', fileName, fileName));
end

% toggle parts
A1 = 1;
B1 = 0;
C1 = 0;

% -------------------- PART A --------------------

w_r0 = r_0./sqrt(r.*(2*r_0 - r)) .* (1 - (a_T * r_0^2/mu) * (3 .* asin(sqrt(r./(2*r_0))) - 3*pi/4 + 2)) + (a_T*r_0/(2*mu)).*(r + 3*r_0);
w_1r = r_0./sqrt(r.*(2*r_0 - r));
w_2r = 1 - (a_T * r_0^2/mu) * (3 .* asin(sqrt(r./(2*r_0))) - 3*pi/4 + 2);
w_3r = (a_T*r_0/(2*mu)).*(r + 3*r_0);

w_p0 = (1 - 4*(3*phi +2) ./ ((3*pi + 8) * B^2)) ./ sin(2*phi + pi/2) + (4*(sin(phi + pi/4)).^2 + 6) ./ ((3*pi + 8) * B^2);
w_1p = 1./sin(2*phi + pi/2);
w_2p = 1 - 4*(3*phi + 2)/((3*pi + 8) * B^2); 
w_3p = 2/((3*pi + 8) * B^2) * (2*sin(phi + pi/4).^2 + 3);

if (A1)

    % plots
    figure(1);
    hold on;

    plot(r,w_r0,'m', 'linewidth', width);
    plot(r,w_1r,'r', 'linewidth', width/2, 'linestyle', '--');
    plot(r,w_2r,'b', 'linewidth', width/2, 'linestyle', '--');
    plot(r,w_3r,'k', 'linewidth', width/2, 'linestyle', '--');

    legend({"{\\it w_{r,0}(r)}", "{\\it w_1(r)}", "{\\it w_2(r)}", "{\\it w_3(r)}"}, 'location', 'northwest', 'orientation', 'vertical');
    legend boxoff;

    grid on;
    ylim([0 3]);
    xlim([r_0 2*r_0]);

    xlabel("{\\it r}");
    ylabel("{\\it w(r)}");
    title(sprintf('{\\it w(r)} vs. {\\it r} with {\\it r_0} = %d, {\\it a_T} = %d, {\\it \\mu} = %d', r_0, a_T, mu));

    figure(2);
    hold on;

    plot(r,w_p0,'m', 'linewidth', width);
    plot(r,w_1p,'r', 'linewidth', width/2, 'linestyle', '--');
    plot(r,w_2p,'b', 'linewidth', width/2, 'linestyle', '--');
    plot(r,w_3p,'k', 'linewidth', width/2, 'linestyle', '--');

    legend({"{\\it w_{p,0}(\\phi)}", "{\\it w_1(\\phi)}", "{\\it w_2(\\phi)}", "{\\it w_3(\\phi)}"}, 'location', 'northwest', 'orientation', 'vertical');
    legend boxoff;

    grid on;
    ylim([0 3]);
    xlim([r_0 2*r_0]);

    xlabel("{\\it r}");
    ylabel("{\\it w(\\phi)}");
    title(sprintf('{\\it w_p(\\phi)} vs. {\\it r} with {\\it r_0} = %d, {\\it a_T} = %d, {\\it \\mu} = %d', r_0, a_T, mu));

    exportPlot('partA_r', 1);
    exportPlot('partA_phi', 2);

    close all;

end

% -------------------- PART B --------------------

y_r = sqrt(mu * (2*r_0 - r)./(r_0 * r)) .* sqrt(1 - w_r0 .^ 2);

g_p = B^4*(3*pi + 8)^2 * (sin(2*phi + pi/2)).^2 - (B^4*(3*pi + 8)^2 - 8*(3*phi + 2)*(3*pi + 8)*B^2 + 16*(3*phi + 2).^2 + 2*(4*(sin(phi + pi/4)).^2 + 6) .* ((3*pi + 8)*B^2 - 4*(3*phi + 2)) .* sin(2*phi + pi/2) + (16*(sin(phi+pi/4)).^4 + 48*(sin(phi+pi/4)).^2 + 36) .* (sin(2*phi+pi/2)).^2);

y_p = sqrt(g_p) ./ (B^2 * (3*pi+8) * sin(2*phi+pi/2)) * sqrt(mu/r_0) .* cot(phi + pi/4);


% y_1r = sqrt(mu * (2*r_0 - r)/(r_0 * r));
% y_1p = sqrt(mu/r_0) .* cot(phi + pi/4);

% y_2r = sqrt(1 - w_r0 .^ 2);
% y_2p = sqrt(mu/r_0) .* cot(phi + pi/4);

if (B1)

    figure(1);
    hold on;

    subplot(1,2,1);
    
    plot(phi, y_r, 'r', 'linewidth', width);

    legend({"{\\it y(r)}"}, 'location', 'northwest', 'orientation', 'vertical');
    legend boxoff;

    grid on;
    ylim([0 0.25]);
    xlim([0 pi/4]);

    xlabel("{\\it \\phi}");
    ylabel("{\\it y(r)}");
    title(sprintf('{\\it y(r)} vs. {\\it \\phi} with {\\it r_0} = %d, {\\it a_T} = %d, {\\it \\mu} = %d', r_0, a_T, mu));

    subplot(1,2,2);

    plot(phi, y_p, 'b', 'linewidth', width);

    legend({"{\\it y_p(\\phi)}"}, 'location', 'northwest', 'orientation', 'vertical');
    legend boxoff;

    grid on;
    ylim([0 0.25]);
    xlim([0 pi/4]);

    xlabel("{\\it \\phi}");
    ylabel("{\\it y_p(\\phi)}");
    title(sprintf('{\\it y_p(\\phi)} vs. {\\it \\phi} with {\\it r_0} = %d, {\\it a_T} = %d, {\\it \\mu} = %d', r_0, a_T, mu));

    exportPlot('partB', 1);

    close all;
end

% -------------------- PART C --------------------

dw_r0 = diff(w_r0)./diff(r);

dw_r1 = -r_0*(r_0 - r)./(2*r_0*r - r.^2).^(3/2) ...
        .* (1 - (a_T * r_0^2/mu) * (3 .* asin(sqrt(r./(2*r_0))) - 3*pi/4 + 2)) ...
        + -3*a_T*r_0./(2*mu*(2*r_0*r - r.^2)) ...
        + a_T*r_0/(2*mu);

dw_p = (a_T*r_0^2*sin(2*phi+pi/2).^3 ...
        - 3*a_T*r_0^2*sin(2*phi+pi/2)-2*mu*r_0*cos(2*phi+pi/2).*(1-a_T*r_0^2/mu * (3*phi+2))) ...
        ./ (2*mu*r_0*sin(2*phi+pi/2).^3);

if (C1)

    figure(1);
    hold on;

    plot(phi(1:end-1),w_r0(1:end-1), 'r', 'linewidth', width);
    plot(phi(1:end-1), dw_p(1:end-1), 'b', 'linewidth', width);
    plot([0 pi/4], [0 0], 'linewidth', width/3, 'linestyle', '--');
    plot(0.2196, 0, 'k','linewidth', width/2);
    plot(0.2196, 0.96, 'k','linewidth', width/2);
    
    legend({"{\\it w_p\\'(\\phi)}","{\\it w_p(\\phi)}"}, 'location', 'northwest', 'orientation', 'vertical');
    legend boxoff;

    grid on;
    
    text(0.2168, 0.25, "{\\it w\\'(0.2168) = 0}");
    text(0.2168, 1.25, "{\\it w(0.2168) = 0.96}");

    ylim([-1 5]);
    xlim([0 pi/4]);

    xlabel("{\\it \\phi}");
    ylabel("{\\it f(\\phi)}");
    title(sprintf("{\\it w_p\\'(\\phi)} vs. {\\it \\phi} with {\\it r_0} = %d, {\\it a_T} = %d, {\\it \\mu} = %d", r_0, a_T, mu));

    exportPlot('partC_r', 1);

    close all;

end