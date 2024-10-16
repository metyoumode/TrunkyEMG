%% TRUNKY KINEMATIC STUDY AND TRAJECTORIES ANALYSIS
clear
close all
clc


%% load mat file
data=load("ModeLog.mat");
%% processing mat data
State=data.logsout{11}.Values;
posM1=data.logsout{17}.Values.PDO_M1.pos_actual_1;
posM2=data.logsout{17}.Values.PDO_M2.pos_actual_2;
posM3=data.logsout{17}.Values.PDO_M3.pos_actual_3;
posM4=data.logsout{17}.Values.PDO_M4.pos_actual_4;
posM5=data.logsout{17}.Values.PDO_M5.pos_actual_5;

StartFrontFlex=[0 0 0]; % 1 è High imp, 2 medium, 3 low
EndFrontFlex=[ 0 0 0];
StartLatFlex=[ 0 0 0];
EndLatFlex=[ 0 0 0];

j=1;
i=1;
startExercise=1;
entrato=0;

for t=1:length(State.Data)
    %I have started a frontal flexion and I am at the beginning off the
    %movement
    if (State.Data(t) == 18 && startExercise==1 )
        StartFrontFlex(j)=t;
        EndFrontFlex(j)=StartFrontFlex(j)+12/0.002;
        j=j+1;
        startExercise=0;
    end
    %I have started a lateral flexion and I am at the beginning off the
    %movement
    if (State.Data(t) == 19 && startExercise==1 )
        StartLatFlex(i)=t;
        EndLatFlex(i)=StartLatFlex(i)+12/0.002;
        i=i+1;
        startExercise=0;
    end
    % if the exercise is finished we enable the starting of the ext one
    if State.Data(t)== 22
        startExercise=1;
        entrato=1;
    end

end

% figure(1)
% plot(State)
%
% for i=1:2
% figure()
% plot(posM1.Time(StartFrontFlex:EndFrontFlex),posM1.Data(StartFrontFlex:EndFrontFlex))
%
% end

%% SELECT THE EXERCISES TO BE ANALYZED
close all
%EXERCISE=1; % uncomment for FRONTAL FLEXION
EXERCISE=2; % uncomment for LATERAL FLEXION


% DIRECT KINEMATIC TRUNKY
% syms qref1 qref2 qref3 qref4 qref5 real; % Uncomment to obtain symbolic
% DK equations

% Homing position RBT, comment to obtain symbolic equations
qref1=deg2rad(0);
qref2=deg2rad(0);
qref3=deg2rad(0);
qref4=deg2rad(0);
qref5=deg2rad(0);

% Rotations around z axes (for the homing config they will be null)
%trotz compute the rotation matrix around the given angle
Tq1=trotz(qref1);
Tq2=trotz(qref2);
Tq3=trotz(qref3);
Tq4=trotz(qref4);
Tq5=trotz(qref5);

% Homogeneous transformation matrices from CREO model
T01= [1	  0	  0	   0
    0	 1	0	0
    0	0	1   0
    0    0   0   1];

T12=Tq1*[-0.853551	0.521010	0.0000000000	240.000e-3
    -0.521010	-0.853551	0.0000000000	0.0000000000e-3
    0.0000000000	0.0000000000	1.00000	75.0000e-3
    0     0     0   1];

T23=Tq2*[0.270600	-0.962692	0.0000000000	240.000e-3
    0.962692	0.270600	0.0000000000	0.0000000000e-3
    0.0000000000	0.0000000000	1.00000	33.0000e-3
    0     0     0   1];

T34=Tq3*[0.0000000000	0.0000000000	1.00000	532.953e-3
    0.0000000000	-1.00000	0.0000000000	0.0000000000e-3
    1.00000	0.0000000000	0.0000000000	82.0000e-3
    0     0     0   1];

T45=Tq4*[-1.00000	0.0000000000	0.0000000000	0.0000000000e-3
    0.0000000000	0.0000000000	-1.00000	0.0000000000e-3
    0.0000000000	-1.00000	0.0000000000	100.500e-3
    0     0     0   1];

T5ee=Tq5*[1.00000	0.0000000000	0.0000000000	100.000e-3
    0.0000000000	1.00000	0.0000000000	0.0000000000e-3
    0.0000000000	0.0000000000	1.00000	0.0000000000e-3
    0     0     0   1];

% Compute the final transformation -> get x y z equations
T0EE=T01*T12*T23*T34*T45*T5ee;
% x = T0EE(1,4);
% y = T0EE(2,4);
% z = T0EE(3,4);


%inizialize posizion matrix for each joints
pM1=zeros(6001,3);
pM2=zeros(6001,3);
pM3=zeros(6001,3);
pM4=zeros(6001,3);
pM5=zeros(6001,3);

% Trajectory analysis, load the joint data obtained from CAD simulations
switch EXERCISE

    case 1
        load q1FFlex.mat; load q2FFlex.mat; load q3FFlex.mat; load q4FFlex.mat; load q5FFlex.mat
        q1Tg=q1FFlex; q2Tg=q2FFlex; q3Tg=q3FFlex; q4Tg=q4FFlex; q5Tg=q5FFlex;
        for j=1:3
            pM1(:,j)=(posM1.Data(StartFrontFlex(j):EndFrontFlex(j))/(2276*2.89)+65);
            pM2(:,j)=(posM2.Data(StartFrontFlex(j):EndFrontFlex(j))/(2276*2.89)-130);
            pM3(:,j)=-(posM3.Data(StartFrontFlex(j):EndFrontFlex(j))/(2276*2.89)+65);
            pM4(:,j)=-(posM4.Data(StartFrontFlex(j):EndFrontFlex(j))/5461);
            pM5(:,j)=-(posM5.Data(StartFrontFlex(j):EndFrontFlex(j))/5461);

        end
    case 2
        load q1LFlex.mat; load q2LFlex.mat; load q3LFlex.mat; load q4LFlex.mat; load q5LFlex.mat
        q1Tg=q1LFlex; q2Tg=q2LFlex; q3Tg=q3LFlex; q4Tg=q4LFlex; q5Tg=q5LFlex;
        for j=1:3
            pM1(:,j)=(posM1.Data(StartLatFlex(j):EndLatFlex(j))/(2276*2.89)+65);
            pM2(:,j)=(posM2.Data(StartLatFlex(j):EndLatFlex(j))/(2276*2.89)-130);
            pM3(:,j)=-(posM3.Data(StartLatFlex(j):EndLatFlex(j))/(2276*2.89)+65);
            pM4(:,j)=-(posM4.Data(StartLatFlex(j):EndLatFlex(j))/5461);
            pM5(:,j)=-posM5.Data(StartLatFlex(j):EndLatFlex(j))/5461;

        end
end


qref1H=deg2rad(65);
qref2H=deg2rad(-130);
qref3H=deg2rad(65);
qref4H=deg2rad(0);
qref5H=deg2rad(90);

[x,y,z]=Joint2Cartesian(q1Tg,q2Tg,q3Tg,q4Tg,q5Tg,T0EE);
[xH,yH,zH]=Joint2Cartesian(pM1(:,1),pM2(:,1),pM3(:,1),pM4(:,1),pM5(:,1),T0EE);
[xM,yM,zM]=Joint2Cartesian(pM1(:,2),pM2(:,2),pM3(:,2),pM4(:,2),pM5(:,2),T0EE);
[xL,yL,zL]=Joint2Cartesian(pM1(:,3),pM2(:,3),pM3(:,3),pM4(:,3),pM5(:,3),T0EE);


figure(2)
plot3(x,y,z,LineWidth=2,Color="#00FF00");%green
hold on
grid on
plot3(xH,yH,zH,LineWidth=2,Color="#FF0000"); %red
plot3(xM,yM,zM,LineWidth=2,Color="#0000FF"); %blue
plot3(xL,yL,zL,LineWidth=2,Color="#FF00FF"); %magenta
% Add a legend
legend('Ideal Trajectory', 'High Impedance Trajectory', 'High Impedance Trajectory', 'High Impedance Trajectory')
xlabel('x [m]','fontname','times','FontSize',12);ylabel('y [m]','fontname','times','FontSize',12);zlabel('z [m]','fontname','times','FontSize',12)
xlim([0.1 0.4]); ylim([-1.3 -0.6]); zlim([0.1 0.3]);
