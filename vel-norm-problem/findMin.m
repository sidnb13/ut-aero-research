k = phi + pi/4;

d_1 = tan(2*k)/(2*mu);
d_2 = (1-a_T*r_0^2/mu * (3*phi+2))./(a_T*r_0^2*(sin(2*k).^3-3));

figure(1);
hold on;

plot(phi,d_1,'r', 'linewidth', width);
plot(phi,d_2,'b', 'linewidth', width);

ylim([-2 0]);
xlim([0 pi/4]);

xlabel("{\\it \\phi}");
ylabel("{\\it d(\\phi)}");
title("Numerical minimization of {\\it w_p(\\phi)} using {\\it d_2(\\phi)} and {\\it d_2(\\phi)}");

legend({"{\\it d_1(\\phi)}", "{\\it d_2(\\phi)}"}, 'location', 'northwest', 'orientation', 'vertical');
legend boxoff;

exportPlot('partC_min', 1);