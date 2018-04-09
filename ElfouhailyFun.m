%% Elfoulhaily (Wind speed variation)
function [ElfouOut] = ElfouhailyFun(in)

a = 0.0081;
b = 1.25;
g = 9.82;
U = [5 10 15];
w = linspace(-0.1,50*pi,5000);
k = in;
Om = [0.84 1 2 3 4 5]; % 0.84: Fully developed sea, 1: "mature" sea, 2-5: "young" sea
Cd = 0.00144; %Drag coefficient
u = [sqrt(Cd)*U(1) sqrt(Cd)*U(2) sqrt(Cd)*U(3)];
a_0 = 0.1733;
a_p = 4;
k_m = 370;
c_m = 0.23;


%prompt = 'Choose value for omega from 1 to 6 where 1 is fully developed ';
%i = input(prompt) %Chose value for Om from 1 to 6
i = 1;
gamma = [1.7 17+6*log10(Om(i))] ; %if Om <= 1 gamma = 1.7, else

if i <= 2
gamma = gamma(1);
else
    gamma = gamma(2);
end

sigma = 0.08.*[1+4*Om(i)^-3];
ap = 0.006.*[Om(i)].^(0.55);

S = zeros(length(U),length(k));

%Choose here windspeed by setting  j = 1,2,3
for j = 1:1:3
  
am1 = (0.01).*(1+3*log(u(j)./c_m));
a_m1 = (0.13/c_m).*u(j)
k_o1 = g.*[(1./U(j)).^2];

k_p = (k_o1.*[Om(i).^2]);
c_p = sqrt(g./[k_o1*Om(i).^2]);
c = sqrt((g./k).*(1+(k./k_m)).^2);

LPM_1 = @(k) exp(-1.25.*((k_o1.*[Om(i).^2])./k).^2);
Gamma_1 = @(k)exp(-(1./(2*sigma.^2))*(sqrt(k./(k_o1.*[Om(i).^2]))-1).^2);
J_p1 = @(Gamma_1) gamma.^Gamma_1(k); 
F_p1 = @(k) LPM_1(k).*J_p1(Gamma_1).*exp(-0.3162.*Om(i).*(sqrt(k./k_p)-1));
F_m1 = @(k) LPM_1(k).*J_p1(Gamma_1).*exp(-(1/4).*((k./k_m)-1).^2);
Bl1 = 0.5.*ap.*(c_p./c).*F_p1(k);
Bh1 = 0.5.*am1.*(c_m./c).*F_m1(k);

    S(j,:) = (Bl1 + Bh1).*(1./(k.^3));
end

ElfouOut = S(:,:);

end

%%
% subplot(2,2,1)
% loglog(k,S(1,:),'b',k,S(2,:),'r',k,S(3,:),'g')
% ylim([10^(-15) 100])
% xlabel('Spatial frequency k [rad/m]')
% ylabel('Variance spectrum S(k) [m^2/(rad/m)]');
% legend('U_1_0 = 5.0','10','15');
% title('Elfouhaily et al., constant Omega_c')
% 
% subplot(2,2,2)
% loglog(k,k.^2.*S(1,:),'b',k,k.^2.*S(2,:),'r',k,k.^2.*S(3,:),'g')
% ylim([10^(-6) 0.1])
% xlabel('Spatial frequency k [rad/m]')
% ylabel('Slope spectrum k^2S(k) [m^2/(rad/m)]');
% title('Elfouhaily et al. constant Omega_c')
% 
% %% Elfouhaliy (Omega variation)
% clc
% clear all
% 
% a = 0.0081;
% b = 1.25;
% g = 9.82;
% U = [5 10 15];
% w = linspace(-0.1,50*pi,5000);
% k = (w.^2)/g;
% Om = [0.84 1 2 3 4 5]; % 0.84: Fully developed sea, 1: "mature" sea, 2-5: "young" sea
% Cd = 0.00144; %Drag coefficient
% u = [sqrt(Cd)*U(1) sqrt(Cd)*U(2) sqrt(Cd)*U(3)];
% a_0 = 0.1733;
% a_p = 4;
% k_m = 370;
% c_m = 0.23;
% j = 2; %Choose here a constant windspeed 
% am1 = (0.01).*(1+3*log(u(j)./c_m));
% a_m1 = (0.13/c_m).*u(j);
% k_o1 = g.*[(1./U(j)).^2];
% 
% %prompt = 'Choose value for omega from 1 to 6 where 1 is fully developed ';
% %i = input(prompt) %Chose value for Om from 1 to 6
% i=1;
% 
%     k_p = (k_o1.*[Om(i).^2]);
%     c_p = sqrt(g./[k_o1*Om(i).^2]);
%     c = sqrt((g./k).*(1+(k./k_m)).^2);
% 
%     gamma = [1.7 17+6*log10(Om(i))] ; %if Om <= 1 gamma = 1.7, else
% 
% 
% S = zeros(length(Om),length(k));
% for i =1:1:6
%     
% if i <= 2
% gamma = gamma(1);
% else
%     k_p = (k_o1.*[Om(i).^2]);
%     c_p = sqrt(g./[k_o1*Om(i).^2]);
%     c = sqrt((g./k).*(1+(k./k_m)).^2);
% 
%     gamma = [1.7 17+6*log10(Om(i))] ; %if Om <= 1 gamma = 1.7, else
%     gamma = gamma(2);
% end
% sigma = 0.08.*[1+4*Om(i)^-3];
% ap = 0.006.*[Om(i)].^(0.55);
% 
% LPM_1 = @(k) exp(-1.25.*((k_o1.*[Om(i).^2])./k).^2);
% Gamma_1 = @(k)exp(-(1./(2*sigma.^2))*(sqrt(k./(k_o1.*[Om(i).^2]))-1).^2);
% J_p1 = @(Gamma_1) gamma.^Gamma_1(k); 
% F_p1 = @(k) LPM_1(k).*J_p1(Gamma_1).*exp(-0.3162.*Om(i).*(sqrt(k./k_p)-1));
% F_m1 = @(k) LPM_1(k).*J_p1(Gamma_1).*exp(-(1/4).*((k./k_m)-1).^2);
% Bl1 = 0.5.*ap.*(c_p./c).*F_p1(k);
% Bh1 = 0.5.*am1.*(c_m./c).*F_m1(k);
% 
% S(i,:) = (Bl1 + Bh1).*(1./(k.^3));
% end
% 
% subplot(2,2,3)
% loglog(k,S(1,:),'b',k,S(2,:),'r',k,S(3,:),'g',k,S(4,:),k,S(5,:),k,S(6,:))
% ylim([10^(-15) 100])
% xlabel('Spatial frequency k [rad/m]')
% ylabel('Variance spectrum S(k) [m^2/(rad/m)]');
% legend('Omega_c = 0.84','1','2','3','4','5');
% title('Elfouhaily et al., constant U_1_0')
% 
% subplot(2,2,4)
% loglog(k,k.^2.*S(1,:),'b',k,k.^2.*S(2,:),'r',k,k.^2.*S(3,:),'g',k,k.^2.*S(4,:),k,k.^2.*S(5,:),k,k.^2.*S(6,:))
% ylim([10^(-6) 0.1])
% xlabel('Spatial frequency k [rad/m]')
% ylabel('Slope spectrum k^2S(k) [m^2/(rad/m)]');
% title('Elfouhaily et al. constant U_1_0')
% 
% end
