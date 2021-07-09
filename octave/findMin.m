source('variables.m');
load('variables.mat');
source('functions.m');

% definitions
k = phi + pi/4;
r = 2*r_0*sin(k).^2;

d1 = tan(2*k)/(2*mu);
d2 = (1-a_T*r_0^2/mu * (3*phi+2))./(a_T*r_0^2*(sin(2*k).^2-3));

% minimization
intersect = find(abs(d1 - d2) <= min(abs(d1 - d2)));

ix = phi(intersect)
iy = mean([d1(intersect) d2(intersect)])
r_min = r(intersect)

% plot

if (1)

    figure(1);
    hold on;

    plot(phi,d1,'r', 'linewidth', width);
    plot(phi,d2,'b', 'linewidth', width);
    plot(ix, iy, 'k', 'linewidth', width);
    text(ix - 2*inc, iy - 0.05, sprintf("{\\it (%d, %d)}", ix, iy));

    ylim([-2 0]);
    xlim([0 pi/4]);

    xlabel("{\\it \\phi}");
    ylabel("{\\it d(\\phi)}");
    title("Numerical minimization of {\\it w_p(\\phi)} using {\\it d_2(\\phi)} and {\\it d_2(\\phi)}");

    legend({"{\\it d_1(\\phi)}", "{\\it d_2(\\phi)}"}, 'location', 'northwest', 'orientation', 'vertical');
    legend boxoff;

    exportPlot('partC_min', 1);

    close all;

end