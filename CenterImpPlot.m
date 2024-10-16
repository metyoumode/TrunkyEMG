close all
clear
clc

data=load("CenterFixImComplete.mat");
%%
clc
posM1=data.data{17}.Values.PDO_M1.pos_actual_1;
posM2=data.data{17}.Values.PDO_M2.pos_actual_2;
posM3=data.data{17}.Values.PDO_M3.pos_actual_3;
posM4=data.data{17}.Values.PDO_M4.pos_actual_4;
posM5=data.data{17}.Values.PDO_M5.pos_actual_5;

% target=data.data{3}.Values.Data;
% time=data.data{3}.Values.Time;
% 
% targetM1=data.data{3}.Values.Data(:,1);
% targetM2=data.data{3}.Values.Data(:,2);
% targetM3=data.data{3}.Values.Data(:,3);
% targetM4=data.data{3}.Values.Data(:,4);
% targetM5=data.data{3}.Values.Data(:,5);
% 


figure(1)
plot(posM1)
% hold on
% grid on
% plot(time,targetM1)

figure(2)
plot(posM2)
% hold on
% grid on
% plot(time,targetM2)

figure(3)
plot(posM3)
% hold on
% grid on
% plot(time,targetM3)
figure(4)
plot(posM4)
% hold on
% grid on
% plot(time,targetM4)

figure(5)
plot(posM5)
% hold on
% grid on
% plot(time,targetM5)







