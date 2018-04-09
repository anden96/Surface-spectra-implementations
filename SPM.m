
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
W = ElfouhailyFun(k);

rcs = 10*log10((8.*k.^4.*sigma^2).*cosd(theta).^4.*(abs(avv).^2).*W)

for i = 1:1:length(U19)
nrcs(:,i) = rcs(kl(i),l(i),theta,avv(theta));
end

nrcsPM = rcsPM(sigma(1),K,L(1),avv(theta)); 

figure(1)
hold on
subplot(2,2,1)
plot(theta,nrcs(:,1),theta,nrcs(:,2),theta,nrcs(:,3))
title('Small pertubations model, Gaussian surface spectra');
ylabel('NRCS (dB)');
xlabel('Angle of incidence (degree)');
legend('k \sigma = 0.1, k_l = 1','k \sigma = 0.2, k_l = 2','k \sigma = 0.3, k_l = 3');
grid on

subplot(2,2,2)
plot(theta,nrcsPM)
title('Small pertubations model, Pierson-Moskowits surface spectra');
ylabel('NRCS (dB)');
xlabel('Angle of incidence (degree)');
grid on

%subplot(2,2,3)
%title('Pierson-Moskowits roughness spectra')
%xlabel('Spatial frequency, k (rad/m)')
%loglog(k,Wpmplot);
%grid on
