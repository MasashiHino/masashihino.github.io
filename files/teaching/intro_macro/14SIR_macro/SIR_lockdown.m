%-----------------
% code for SIR+Death
% initially provided by Taisuke Nakata
% modified by Masashi Hino
%-----------------

clc;clear all;
close all;
T=150;
beta = 0.5852;

%% setting parameters
gamma_D = 7/18*0.005;% probability of death within a week
gamma_R   =7/18*(1-0.005) ;% probability of recovery within a week
gamma_S = 7/18*0.0;% probability of getting covid "again".

S=zeros(1,T);%population of susceptible at all t 
I=zeros(1,T);%population of infectious at all t
R=zeros(1,T);% population of recovered at all t
D=zeros(1,T);% number of death at all t

X0    = beta/gamma_R; % effective reproduction number (not used)
Sstar = 1/X0;
Rstar = 1-Sstar;

% Initial Condition
I(1) = 10^(-6);
S(1) = 1 - I(1);
R(1) = 0;
D(1) =0;

T1=61;% timing when lockdown begins
T2=65;% timing when lockdown ends 

t_L=zeros(1,T);
t_L(T1:T2)=1;
% lockdown
for i=1:2
    if i==1
        L=zeros(1,T);
    else
        L(T1:T2)=0.2;
    end

for t=2:T
	S(t) = S(t-1) - beta*(1-L(t-1))*S(t-1)*(1-L(t-1))*I(t-1)+gamma_S*R(t-1);	
	I(t) = I(t-1) + beta*(1-L(t-1))*S(t-1)*(1-L(t-1))*I(t-1)- gamma_R*I(t-1)-gamma_D*I(t-1);	
	R(t) = R(t-1) + gamma_R*I(t-1) - gamma_S*R(t-1);	
    D(t) = D(t-1) + gamma_D*I(t-1);
end
    Trans=beta*S.*I;

    %%
figure(1)    
hold on
if i==1
    area(t_L*100)
end
plot([1:T],Trans*100,'k','LineWidth',1.5)
ylim([min(Trans*100),max(Trans*100)+0.5])
title('�V�K�����Ґ�(T)','FontSize',20)
xlabel('����(�T)','FontSize',14,'interpreter','latex')
ylabel('(%)','FontSize',12)
%%

figure(2)
subplot(2,2,1)
hold on
if i==1
    area(t_L*100)
end
plot([1:T],I*100,'k','LineWidth',1.5);
ylim([min(I*100),max(I*100)+2])
title('������(I)','FontSize',20,'interpreter','latex')
xlabel('����(�T)','FontSize',14,'interpreter','latex')
ylabel('(%)','FontSize',12)

subplot(2,2,2)
hold on
if i==1
    area(t_L*100)
end
plot([1:T],S*100,'k','LineWidth',1.5);
ylim([min(S*100),100])
%plot([1:T],ones(1,T)*Sstar,'k--','LineWidth',1.5);
title('���N�Ȑl(S)','FontSize',20,'interpreter','latex')


subplot(2,2,3)
hold on
if i==1
    area(t_L*100)
end
plot([1:T],R*100,'k','LineWidth',1.5);
ylim([min(R*100),max(R*100)+5])
%plot([1:T],ones(1,T)*Rstar,'k--','LineWidth',1.5);
title('�񕜂����l(R)','FontSize',20,'interpreter','latex')

subplot(2,2,4)
hold on
if i==1
    area(t_L*100)
end
plot([1:T],D*100,'k','LineWidth',1.5);
ylim([min(D*100),max(D*100)+0.01])
title('���S�Ґ�(D)','FontSize',20,'interpreter','latex')

end
%%


