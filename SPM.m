% Implementation of elfouhaily spectra on small pertubations model.
U10 = [5 10 15];
U19 = 1.026*U10;
epsr = 20;
g = 9.82;
theta = linspace(0,80);
sigma = 0.4045
kl = [1 2 3];
ksigma = [0.1 0.2 0.3];
k = [ksigma./sigma];

N = (epsr*cosd(theta)+sqrt(epsr-sind(theta).^2)).^2;
T = sind(theta).^2-epsr.*(1+sind(theta).^2);
avv = (epsr-1).*(T./N);
W1 = ElfouhailyFun(2.*k(1).*sind(theta)); 
W2 = ElfouhailyFun(2.*k(2).*sind(theta));
W3 = ElfouhailyFun(2.*k(3).*sind(theta));

rcs1 = 10*log10((8.*k(1).^4.*sigma(1)^2).*cosd(theta).^4.*(abs(avv).^2).*W1)
rcs2 = 10*log10((8.*k(1).^4.*sigma(1)^2).*cosd(theta).^4.*(abs(avv).^2).*W2)
rcs3 = 10*log10((8.*k(1).^4.*sigma(1)^2).*cosd(theta).^4.*(abs(avv).^2).*W3)

figure(1)
hold on
subplot(2,2,1)
plot(theta,rcs1(1,:),theta,rcs1(2,:),theta,rcs1(3,:))
ylim([-75 0])
grid on
subplot(2,2,2)
plot(theta,rcs2(1,:),theta,rcs2(2,:),theta,rcs2(3,:))
ylim([-75 0])
grid on
subplot(2,2,3)
plot(theta,rcs3(1,:),theta,rcs3(2,:),theta,rcs3(3,:))
title('Small pertubations model, Elfouhaily surface spectra');
ylabel('NRCS (dB)');
xlabel('Angle of incidence (degree)');
legend('U_1_0 = 5m/s','U_1_0 = 10 m/s','U_1_5 = 15 m/s');
ylim([-75 0])
grid on

%subplot(2,2,3)
%title('Pierson-Moskowits roughness spectra')
%xlabel('Spatial frequency, k (rad/m)')
%loglog(k,Wpmplot);
%grid on
