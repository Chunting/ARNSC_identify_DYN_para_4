t = timer.Data;
mT(1)= 0;
bT(1) =0;
IT(1) = 0;
figure(1);
subplot(3,1,1)
semilogx(t,mT,'LineWidth',2)
axis([0 20 0 300]);
ylabel('Mass (kg)','fontsize',12);
xlabel('Time (s)','fontsize',12);
subplot(3,1,2)
semilogx(t,bT,'LineWidth',2);
ylabel('Mass center (m)','fontsize',12);
axis([0 20 -0.5 0.5]);

%figure(6);
subplot(3,1,3);
IT(1)=0;

% plot(t(1:4001),IT(1:4001),'LineWidth',2);
semilogx(t,IT,'LineWidth',2);
% grid on
axis([0 20 -100 300]);
%set(gca,'xticklabel',0:2:20);
%set(gca,'YTick',100:100:400);
xlabel('Time (s)','fontsize',12);
ylabel('Inertia (kg.m^2)','fontsize',12);
%title(' Inertia of the target');
%{
figure(1);

plot(t,mT,'LineWidth',2)
%plot(t,mT,'LineWidth',2)
% grid on
axis([0 20 -50 250]);
set(gca,'xticklabel',0:2:20);
%set(gca,'YTick',210:10:230);
ylabel('Mass (kg)','fontsize',12);
xlabel('Time (s)','fontsize',12);
%figure(5);
figure(2);
bT(1) = 0;

plot(t,bT,'LineWidth',2);

% grid on
ylabel('Mass center (m)','fontsize',12);
xlabel('Time (s)','fontsize',12);
axis([0 20 -0.5 0.5]);
set(gca,'xticklabel',0:2:20);
%set(gca,'YTick',0.25:0.05:0.35);
%figure(6);
figure(3);
IT(1)=0;

plot(t,IT,'LineWidth',2);
%plot(t,IT,'LineWidth',2);
% grid on
axis([0 20 -100 500]);
set(gca,'xticklabel',0:2:20);

xlabel('Time (s)','fontsize',12);
ylabel('Inertia (kg.m^2)','fontsize',12);

figure(4);
plot(t,mT,'LineWidth',2)
axis([0 0.1 -100 300]);
set(gca,'xticklabel',0:0.01:0.1);

figure(5);
bT(1) = 0;
plot(t,bT,'LineWidth',2);
axis([0 0.1 -0.5 0.5]);
set(gca,'xticklabel',0:0.01:0.1);

figure(6);
IT(1)=0;
plot(t,IT,'LineWidth',2);
axis([0 0.1 -100 300]);
set(gca,'xticklabel',0:0.01:0.1);
%}

