close all
clear
clc

data=load("ImpFixJ4.mat");
%%
clc
posM1=data.data{17}.Values.PDO_M1.pos_actual_1;
posM2=data.data{17}.Values.PDO_M2.pos_actual_2;
posM3=data.data{17}.Values.PDO_M3.pos_actual_3;
posM4=data.data{17}.Values.PDO_M4.pos_actual_4;
posM5=data.data{17}.Values.PDO_M5.pos_actual_5;

target=data.data{3}.Values.Data;
time=data.data{3}.Values.Time;

targetM1=data.data{3}.Values.Data(:,1)/(2276*2.89);
targetM2=data.data{3}.Values.Data(:,2)/(2276*2.89);
targetM3=data.data{3}.Values.Data(:,3)/(2276*2.89);
targetM4=data.data{3}.Values.Data(:,4)/5461;
targetM5=data.data{3}.Values.Data(:,5)/5461;


%
% targetM1=target(:,1);
% targetM2=target(:,2);
% targetM3=target(:,3);
% targetM4=target(:,4);
% targetM5=target(:,5);



%
% State=data.data{11}.Values;
%
% StartFrontFlex=[0 0 0]; % 1 è High imp, 2 medium, 3 low
% EndFrontFlex=[ 0 0 0];
% StartLatFlex=[ 0 0 0];
% EndLatFlex=[ 0 0 0];
%
% repetition=0;
% j=1;
% i=1;
% startExercise=1;
% entrato=0;
%
% for t=1:length(State.Data)
%     %se ho iniziato una flessione e sono sotto il numero di ripetizioni
%     if (State.Data(t) == 18 && startExercise==1 )
%        StartFrontFlex(j)=State.Time(t);
%         EndFrontFlex(j)=StartFrontFlex(j)+40;
%         j=j+1;
%         startExercise=0;
%         repetition=repetition+1;
%     end
%     %se ho iniziato una flessione laterale e sono sotto il numero di ripetizioni
%     if (State.Data(t) == 19 && startExercise==1 )
%         StartLatFlex(i)=State.Time(t);
%         EndLatFlex(i)=StartLatFlex(i)+40;
%         i=i+1;
%         startExercise=0;
%         repetition=repetition+1;
%     end
%      % if the exercise is finished we enable the starting of the ext one
%     if State.Data(t)== 22
%         startExercise=1;
%         entrato=1;
%     end
%
% end
%

%% inc to deg conversion
pM1=(posM1.Data)/(2276*2.89);
pM2=(posM2.Data)/(2276*2.89);
pM3=(posM3.Data)/(2276*2.89);
pM4=(posM4.Data)/5461;
pM5=(posM5.Data)/5461;
%% plot position
close all
figure(1)
plot(posM1.Time,pM1)
grid on
set(gca, 'FontSize', 16);  % Set font size for tick labels


figure(2)
plot(posM2.Time,pM2)
grid on

figure(3)
plot(posM3.Time,pM3)
grid on

figure(4)
plot(posM4.Time,pM4,'LineWidth',2)
grid on
set(gca, 'FontSize', 16);  % Set font size for tick labels
ylabel('Position [°]',FontSize=14)
xlabel('Time [s]',FontSize=14)
set(gca, 'FontSize', 14);  % Set font size for tick labels

figure(5)
plot(posM5.Time,pM5,'LineWidth',2)
grid on
set(gca, 'FontSize', 16);  % Set font size for tick labels
ylabel('Position [°]',FontSize=14)
xlabel('Time [s]',FontSize=14)
set(gca, 'FontSize', 14);  % Set font size for tick labels


%% plot position and target
close all
figure(1)
plot(posM1.time,pM1,'LineWidth',4.5)
hold on
grid on
plot(time,targetM1,'LineWidth',1.5,Color="#FF00FF")

figure(2)
plot(posM2.time,pM2,'LineWidth',4.5)
hold on
grid on
plot(time,targetM2,'LineWidth',1.5,Color="#FF00FF")

figure(3)
plot(posM3.time,pM3,'LineWidth',4.5)
hold on
grid on
plot(time,targetM3,'LineWidth',1.5,Color="#FF00FF")

figure(4)
plot(posM4.time,pM4,'LineWidth',4.5)
hold on
grid on
plot(time,targetM4,'LineWidth',1.5,Color="#FF00FF")
ylabel('Position [°]',FontSize=14)
xlabel('Time [s]',FontSize=14)
set(gca, 'FontSize', 14);  % Set font size for tick labels

figure(5)
plot(posM5.time,pM5,'LineWidth',4.5)
hold on
grid on
plot(time,targetM5,LineWidth=1.5,Color="#FF00FF")
ylabel('Position [°]',FontSize=14)
xlabel('Time [s]',FontSize=14)
set(gca, 'FontSize', 14);  % Set font size for tick labels








