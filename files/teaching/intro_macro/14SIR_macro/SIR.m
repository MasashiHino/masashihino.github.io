%-----------------
% code for SIRD model
% initially provided by Tai Nakata
% modified by Masashi Hino
%-----------------

clc;clear all;
close all;

% Parameter setting
T=150; % end of period
beta = 0.5852; % infection rate
gamma_D = 7/18*0.005; %death probability 
gamma_R   =7/18*(1-0.005) ; % recovery probability

S=zeros(1,T);%population of susceptible at all t 
I=zeros(1,T);%population of infectious at all t
R=zeros(1,T);% population of recovered at all t
D=zeros(1,T);% number of death at all t

X0    = beta/gamma_R;
Sstar = 1/X0;
Rstar = 1-Sstar;

% Initial Condition
I(1) = 10^(-6);
S(1) = 1 - I(1);
R(1) = 0;
D(1) =0;

for t=2:T
	S(t) = S(t-1) - beta*S(t-1)*I(t-1);	
	I(t) = I(t-1) + beta*S(t-1)*I(t-1) - gamma_R*I(t-1)-gamma_D*I(t-1);	
	R(t) = R(t-1) + gamma_R*I(t-1);	
    D(t) = D(t-1) + gamma_D*I(t-1);
end
    Trans=beta*S.*I;

%%
figure    
plot([1:T],Trans*100,'k','LineWidth',1.5)
title('êVãKä¥êıé“êî(T)','FontSize',20)
xlabel('éûä‘(èT)','FontSize',14,'interpreter','latex')

%%
figure
subplot(2,2,1)
plot([1:T],I*100,'k','LineWidth',1.5);
title('ä¥êıé“(I)','FontSize',20,'interpreter','latex')
xlabel('éûä‘(èT)','FontSize',14,'interpreter','latex')

subplot(2,2,2)
plot([1:T],S*100,'k','LineWidth',1.5);
hold on
%plot([1:T],ones(1,T)*Sstar,'k--','LineWidth',1.5);
title('åíçNÇ»êl(S)','FontSize',20,'interpreter','latex')




subplot(2,2,3)
plot([1:T],R*100,'k','LineWidth',1.5);
hold on
%plot([1:T],ones(1,T)*Rstar,'k--','LineWidth',1.5);
title('âÒïúÇµÇΩêl(R)','FontSize',20,'interpreter','latex')

subplot(2,2,4)
plot([1:T],D*100,'k','LineWidth',1.5);
hold on
title('éÄñSé“êî(D)','FontSize',20,'interpreter','latex')


%%

% subplot(2,3,5)
% plot([1:T],X0*ones(1,t),'k','LineWidth',1.5);
% title('Basic reproductive number','FontSize',20)
% xlabel('Time')
% subplot(2,3,6)
% plot([1:T],X0*S,'k','LineWidth',1.5);
% xline(54)
% title('Effective reproductive number','FontSize',20)


%print(figure(1),'-dpdf','sir.pdf')
print(figure(1),'-dpng','sir.png')

