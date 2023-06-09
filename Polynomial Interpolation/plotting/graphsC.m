y1 = @(x) 1*x.^0  + 2.22045e-16*x.^1  + -3.54298*x.^2  + 2.7465*x.^4 
x = -1 : 0.01 : 1;
plot(x, y1(x))

hold on

y2 = @(x) 0.730822*x.^0  + -2.63684e-15*x.^1  + -4.81162*x.^2  ...
    + 3.06605e-14*x.^3  + 12.6193*x.^4  + -5.5704e-14*x.^5  + -14.0024*x.^6 ...
    + 9.10478e-14*x.^7  + 5.51277*x.^8  + -3.64196e-14*x.^9
plot(x, y2(x))

y3 = @(x) 1*x.^0  + 1.4988e-15*x.^1  + -17.3641*x.^2  + -7.10543e-15*x.^3  ...
    + 149.027*x.^4  + 3.69482e-13*x.^5  + -646.864*x.^6  + 1.76215e-12*x.^7  ...
    + 1510.61*x.^8  + 2.84217e-12*x.^9  + -1927.18*x.^10  + -9.09495e-13*x.^11  ...
    + 1264.42*x.^12  + -5.68434e-14*x.^13  + -333.619*x.^14
plot(x, y3(x))

y4 = @(x) 0.96241*x.^0  + 1.33568e-15*x.^1  + -16.5422*x.^2  + 6.65006e-14*x.^3  ...
    + 165.458*x.^4  + -1.95768e-12*x.^5  + -960.825*x.^6  + -1.72285e-11*x.^7  ...
    + 3379.02*x.^8  + -7.60981e-11*x.^9  + -7413.45*x.^10  + -6.12795e-11*x.^11  ...
    + 10195.5*x.^12  + -1.17966e-10*x.^13  + -8534.89*x.^14  + 3.93281e-12*x.^15  ...
    + 3973.16*x.^16  + -1.66573e-11*x.^17  + -788.326*x.^18  + 3.19307e-12*x.^19
plot(x, y4(x))

hold off