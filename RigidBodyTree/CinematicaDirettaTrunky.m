%% TRUNKY KINEMATIC STUDY AND TRAJECTORIES ANALYSIS
clear 
close all
clc

% SELECT THE EXERCISES TO BE ANALYZED

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
x = T0EE(1,4);
y = T0EE(2,4);
z = T0EE(3,4);

% Build the kinematic chain model with RBT
Trunky = rigidBodyTree;
Trunky.Gravity=[0 0 -9.81];
bodybase = rigidBody('base');
basename = Trunky.BaseName;

% Joint1: planar complex
body1 = rigidBody('body1');
jnt1 = rigidBodyJoint('jnt1', 'revolute');
setFixedTransform(jnt1, T01);
body1.Joint = jnt1;
addBody(Trunky, body1, basename);

% Joint2: planar complex
body2 = rigidBody('body2');
jnt2 = rigidBodyJoint('jnt2', 'revolute');
setFixedTransform(jnt2, T12);
body2.Joint = jnt2;
addBody(Trunky, body2, 'body1');

% Joint3: planar complex
body3 = rigidBody('body3');
jnt3 = rigidBodyJoint('jnt3', 'revolute');
setFixedTransform(jnt3, T23);
body3.Joint = jnt3;
addBody(Trunky, body3, 'body2');

% Joint4: trunk lateral flexion
body4 = rigidBody('body4');
jnt4 = rigidBodyJoint('jnt4', 'revolute');
setFixedTransform(jnt4, T34);
body4.Joint = jnt4;
addBody(Trunky, body4, 'body3');

% Joint5: trunk frontal flexion
body5 = rigidBody('body5');
jnt5 = rigidBodyJoint('jnt5', 'revolute');
setFixedTransform(jnt5, T45);
body5.Joint = jnt5;
addBody(Trunky, body5, 'body4');

% End-effector
bodyEndEffector = rigidBody('endeffector');
setFixedTransform(bodyEndEffector.Joint, T5ee);
addBody(Trunky, bodyEndEffector, 'body5');

show(Trunky)
xlim([-0.5 0.5]); ylim([-1.5 0.5]); zlim([-0.5 0.5]);
xlabel('x [m]','fontname','times','FontSize',12);ylabel('y [m]','fontname','times','FontSize',12);zlabel('z [m]','fontname','times','FontSize',12);
title('Trunky RBT');

% Trajectory analysis, load the joint data obtained from CAD simulations
switch EXERCISE
    case 1
        load q1FFlex.mat; load q2FFlex.mat; load q3FFlex.mat; load q4FFlex.mat; load q5FFlex.mat
        q1Tg=q1FFlex; q2Tg=q2FFlex; q3Tg=q3FFlex; q4Tg=q4FFlex; q5Tg=q5FFlex;
    case 2
        load q1LFlex.mat; load q2LFlex.mat; load q3LFlex.mat; load q4LFlex.mat; load q5LFlex.mat
        q1Tg=q1LFlex; q2Tg=q2LFlex; q3Tg=q3LFlex; q4Tg=q4LFlex; q5Tg=q5LFlex;
end


% Corresponig Homing CREO values
qref1H=deg2rad(65);
qref2H=deg2rad(-130);
qref3H=deg2rad(65);
qref4H=deg2rad(0);
qref5H=deg2rad(90);
i=1;

% Compute the Cartesian Trajectory
while(i<=length(q1Tg))

    qref1=deg2rad(q1Tg(i))-qref1H;
    qref2=deg2rad(q2Tg(i))-qref2H;
    qref3=deg2rad(-q3Tg(i))-qref3H;
    qref4=deg2rad(q4Tg(i))-qref4H;
    qref5=deg2rad(q5Tg(i))-qref5H;

    x(i)=(6*cos(qref1))/25 + (sin(qref5)*(cos(qref3)*(((1353*cos(qref2))/5000 - (8671158664945115*sin(qref2))/9007199254740992)*((1922025982770857*cos(qref1))/2251799813685248 - (52101*sin(qref1))/100000) - ((8671158664945115*cos(qref2))/9007199254740992 + (1353*sin(qref2))/5000)*((52101*cos(qref1))/100000 + (1922025982770857*sin(qref1))/2251799813685248)) - sin(qref3)*(((1353*cos(qref2))/5000 - (8671158664945115*sin(qref2))/9007199254740992)*((52101*cos(qref1))/100000 + (1922025982770857*sin(qref1))/2251799813685248) + ((8671158664945115*cos(qref2))/9007199254740992 + (1353*sin(qref2))/5000)*((1922025982770857*cos(qref1))/2251799813685248 - (52101*sin(qref1))/100000))))/10 - (6*cos(qref2)*((1922025982770857*cos(qref1))/2251799813685248 - (52101*sin(qref1))/100000))/25 + (6*sin(qref2)*((52101*cos(qref1))/100000 + (1922025982770857*sin(qref1))/2251799813685248))/25 - (89150584211147589*cos(qref3)*(((1353*cos(qref2))/5000 - (8671158664945115*sin(qref2))/9007199254740992)*((1922025982770857*cos(qref1))/2251799813685248 - (52101*sin(qref1))/100000) - ((8671158664945115*cos(qref2))/9007199254740992 + (1353*sin(qref2))/5000)*((52101*cos(qref1))/100000 + (1922025982770857*sin(qref1))/2251799813685248)))/140737488355328000 + (89150584211147589*sin(qref3)*(((1353*cos(qref2))/5000 - (8671158664945115*sin(qref2))/9007199254740992)*((52101*cos(qref1))/100000 + (1922025982770857*sin(qref1))/2251799813685248) + ((8671158664945115*cos(qref2))/9007199254740992 + (1353*sin(qref2))/5000)*((1922025982770857*cos(qref1))/2251799813685248 - (52101*sin(qref1))/100000)))/140737488355328000 + (cos(qref5)*sin(qref4)*(cos(qref3)*(((1353*cos(qref2))/5000 - (8671158664945115*sin(qref2))/9007199254740992)*((52101*cos(qref1))/100000 + (1922025982770857*sin(qref1))/2251799813685248) + ((8671158664945115*cos(qref2))/9007199254740992 + (1353*sin(qref2))/5000)*((1922025982770857*cos(qref1))/2251799813685248 - (52101*sin(qref1))/100000)) + sin(qref3)*(((1353*cos(qref2))/5000 - (8671158664945115*sin(qref2))/9007199254740992)*((1922025982770857*cos(qref1))/2251799813685248 - (52101*sin(qref1))/100000) - ((8671158664945115*cos(qref2))/9007199254740992 + (1353*sin(qref2))/5000)*((52101*cos(qref1))/100000 + (1922025982770857*sin(qref1))/2251799813685248))))/10;
    y(i)=(6*sin(qref1))/25 + (sin(qref5)*(cos(qref3)*(((1353*cos(qref2))/5000 - (8671158664945115*sin(qref2))/9007199254740992)*((52101*cos(qref1))/100000 + (1922025982770857*sin(qref1))/2251799813685248) + ((8671158664945115*cos(qref2))/9007199254740992 + (1353*sin(qref2))/5000)*((1922025982770857*cos(qref1))/2251799813685248 - (52101*sin(qref1))/100000)) + sin(qref3)*(((1353*cos(qref2))/5000 - (8671158664945115*sin(qref2))/9007199254740992)*((1922025982770857*cos(qref1))/2251799813685248 - (52101*sin(qref1))/100000) - ((8671158664945115*cos(qref2))/9007199254740992 + (1353*sin(qref2))/5000)*((52101*cos(qref1))/100000 + (1922025982770857*sin(qref1))/2251799813685248))))/10 - (6*cos(qref2)*((52101*cos(qref1))/100000 + (1922025982770857*sin(qref1))/2251799813685248))/25 - (6*sin(qref2)*((1922025982770857*cos(qref1))/2251799813685248 - (52101*sin(qref1))/100000))/25 - (89150584211147589*cos(qref3)*(((1353*cos(qref2))/5000 - (8671158664945115*sin(qref2))/9007199254740992)*((52101*cos(qref1))/100000 + (1922025982770857*sin(qref1))/2251799813685248) + ((8671158664945115*cos(qref2))/9007199254740992 + (1353*sin(qref2))/5000)*((1922025982770857*cos(qref1))/2251799813685248 - (52101*sin(qref1))/100000)))/140737488355328000 - (89150584211147589*sin(qref3)*(((1353*cos(qref2))/5000 - (8671158664945115*sin(qref2))/9007199254740992)*((1922025982770857*cos(qref1))/2251799813685248 - (52101*sin(qref1))/100000) - ((8671158664945115*cos(qref2))/9007199254740992 + (1353*sin(qref2))/5000)*((52101*cos(qref1))/100000 + (1922025982770857*sin(qref1))/2251799813685248)))/140737488355328000 - (cos(qref5)*sin(qref4)*(cos(qref3)*(((1353*cos(qref2))/5000 - (8671158664945115*sin(qref2))/9007199254740992)*((1922025982770857*cos(qref1))/2251799813685248 - (52101*sin(qref1))/100000) - ((8671158664945115*cos(qref2))/9007199254740992 + (1353*sin(qref2))/5000)*((52101*cos(qref1))/100000 + (1922025982770857*sin(qref1))/2251799813685248)) - sin(qref3)*(((1353*cos(qref2))/5000 - (8671158664945115*sin(qref2))/9007199254740992)*((52101*cos(qref1))/100000 + (1922025982770857*sin(qref1))/2251799813685248) + ((8671158664945115*cos(qref2))/9007199254740992 + (1353*sin(qref2))/5000)*((1922025982770857*cos(qref1))/2251799813685248 - (52101*sin(qref1))/100000))))/10;
    z(i)=19/100 - (cos(qref4)*cos(qref5))/10;
    i=i+1;
end

% Compute the exercise starting configuration
qref1StFF=deg2rad(q1Tg(1))-qref1H;
qref2StFF=deg2rad(q2Tg(1))-qref2H;
qref3StFF=deg2rad(-q3Tg(1))-qref3H;
qref4StFF=deg2rad(q4Tg(1))-qref4H;
qref5StFF=deg2rad(q5Tg(1))-qref5H;
TgStartPosConfiguration=[qref1StFF qref2StFF qref3StFF qref4StFF qref5StFF]';

% Visualization of the performed motion
figure();
Trunky.DataFormat='column';
show(Trunky,TgStartPosConfiguration,Visuals="off");

hold on
plot3(x,y,z,LineWidth=2);
switch EXERCISE
    case 1
        title('Frontal Flexion');
        framesPerSecond = 1000;
    case 2
        title('Lateral Flexion');
        framesPerSecond = 10;
end

xlabel('x [m]','fontname','times','FontSize',12);ylabel('y [m]','fontname','times','FontSize',12);zlabel('z [m]','fontname','times','FontSize',12)
xlim([-0.5 0.5]); ylim([-1.5 0.5]); zlim([-0.5 0.5]);

rate = rateControl(framesPerSecond);
Trunky.DataFormat='row';
Trunky.Gravity=[0 0 -9.81];

count=length(q1Tg);

for i = 1:count
    qref1=(deg2rad(q1Tg(i))-qref1H);
    qref2=(deg2rad(q2Tg(i))-qref2H);
    qref3=deg2rad(-q3Tg(i))-qref3H;
    qref4=deg2rad(q4Tg(i))-qref4H;
    qref5=deg2rad(q5Tg(i))-qref5H;
    show(Trunky,[qref1 qref2 qref3 qref4 qref5],'PreservePlot',false, Visuals="off");
    drawnow
    waitfor(rate);
end