figure(20)
clf
cla

AA=load('save_ex01_M24_history068.txt');

t=AA(:,2);
rb1=AA(:,3);
rb2=AA(:,4);
rb3=AA(:,5);
Tmax=AA(:,12);
equiv=AA(:,15);
rmax=AA(:,16);
total_area=AA(:,17);
amass_total=AA(:,18);
nn=length(t);
rb1_mean=mean(rb1(ceil(0.5*nn):nn));
rb2_mean=mean(rb2(ceil(0.5*nn):nn));
rb3_mean=mean(rb3(ceil(0.5*nn):nn));
[rb1_mean rb2_mean rb3_mean]'

plot(t,rb1,'k','LineWidth',[2])
hold on

BB=load('save_ex01_M24_history020.txt');

t3=BB(:,2);
rb3save=BB(:,3);
rb4save=BB(:,4);
rb5save=BB(:,5);
T5max=BB(:,12);
equiv3=BB(:,15);
rmax3=BB(:,16);
total_area3=BB(:,17);
amass_total3=BB(:,18);

plot(t3,rb3save,'b','LineWidth',[2])
nn=length(t3);
rb1_mean=mean(rb3save(ceil(0.5*nn):nn));
rb2_mean=mean(rb4save(ceil(0.5*nn):nn));
rb3_mean=mean(rb5save(ceil(0.5*nn):nn));
[rb3_mean rb2_mean rb3_mean]'


hold off
set(gca,'FontSize',[18],'FontWeight','bold')
xlabel('t (s)')
ylabel('r_b (cm/s)')
title('M24')
h3=legend('68 atm','30 atm','Location','NE');
set(h3,'FontSize',[14],'FontWeight','bold')
grid on

print ex01_M24_fig01 -dpdf

figure(21)
clf
cla
plot(t,Tmax,'k','LineWidth',[2])
hold on
plot(t3,T5max,'b','LineWidth',[2])
hold off
set(gca,'FontSize',[18],'FontWeight','bold')
xlabel('t (s)')
ylabel('T_{surf,max} (K)')
grid on
print ex01_M24_fig02 -dpdf

figure(22)
plot(t,rmax,'k','LineWidth',[2])
hold on
plot(t3,rmax3,'b','LineWidth',[2])
hold off
set(gca,'FontSize',[18],'FontWeight','bold')
xlabel('t (s)')
ylabel('Max Total Heat Release')
grid on
print ex01_M24_fig03 -dpdf

figure(23)
plot(t,total_area,'k','LineWidth',[2])
hold on
plot(t3,total_area3,'b','LineWidth',[2])
hold off
set(gca,'FontSize',[18],'FontWeight','bold')
xlabel('t (s)')
ylabel('Surface Area')
grid on
print ex01_M24_fig04 -dpdf

figure(24)
plot(t,amass_total,'k','LineWidth',[2])
hold on
plot(t3,amass_total3,'b','LineWidth',[2])
hold off
set(gca,'FontSize',[18],'FontWeight','bold')
xlabel('t (s)')
ylabel('Mass Flux (gm/cm^2-s)')
grid on
print ex01_M24_fig05 -dpdf

figure(25)
plot(t,equiv,'k','LineWidth',[2])
hold on
plot(t3,equiv3,'b','LineWidth',[2])
hold off
set(gca,'FontSize',[18],'FontWeight','bold')
xlabel('t (s)')
ylabel('Equiv Ratio')
grid on
print ex01_M24_fig06 -dpdf

