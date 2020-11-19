%以下程序是机械臂正动力学程序
%开始运行之前需要得到关节力矩tau的数据才能执行
thetalist = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1];
dthetalist = [0.1; 0.2; 0.3; 0.4; 0.5; 0.6];
ddthetalist = [0 ; 0 ; 0 ; 0 ; 0 ; 0];
taumat = tau.Data;

g = [0; 0; -9.8];
Ftipmat = zeros(size(taumat, 1), 6);

dt = 0.1;
intRes = 8;    %积分分辨率是在每个时间步之间进行积分（Euler）的次数

[m,n] = size(taumat);
disp('程序执行时间约需要 ')
time = m * intRes / 50 

M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
M45 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
M56 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
M67 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
Mlist = cat(3, M01, M12, M23, M34, M45, M56, M67); 

%连杆i的空间惯性矩阵
G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7]);
G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393]);
G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
G4 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
G5 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
G6 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
Glist = cat(3, G1, G2, G3, G4, G5, G6);

w1 = [0;0;1]; %沿关节轴的单位向量 相对于基坐标系
w2 = [0;-1;0];
w3 = [0;-1;0];
w4 = [1;0;0];
w5 = [0;1;0];
w6 = [0;0;-1];

q1 = [0; 0; 0]; %qn是关节轴上任一点 坐标值在基坐标系中进行度量 
q2 = [169.9876; 0; -0.0121] ;
q3 = [170.0136; 0; 613.9879] ;
q4 = [0; -0.0112; 814.0441] ;
q5 = [810.0109; 0; 814.0039] ;
q6 = [810.0018; 0; 0];

v1 = -cross(w1,q1); %线速度
v2 = -cross(w2,q2);
v3 = -cross(w3,q3);
v4 = -cross(w4,q4);
v5 = -cross(w5,q5);
v6 = -cross(w6,q6);

%运动旋量S1-S6
Slist = [[w1;v1]  [w2;v2]  [w3;v3]  [w4;v4]  [w5;v5]  [w6;v6]];

taumat = taumat';
Ftipmat = Ftipmat';
thetamat = taumat;
thetamat(:, 1) = thetalist;
dthetamat = taumat;
dthetamat(:, 1) = dthetalist;
ddthetamat(:, 1) = ddthetalist;
for i = 1: size(taumat, 2) - 1
    for j = 1: intRes
       ddthetalist = ForwardDynamics(thetalist, dthetalist, taumat(:,i), g,  Ftipmat(:, i), Mlist, Glist, Slist);     
       [thetalist, dthetalist] = EulerStep(thetalist, dthetalist, ddthetalist, dt / intRes);
    end
    thetamat(:, i + 1) = thetalist;
    dthetamat(:, i + 1) = dthetalist;
    ddthetamat(:, i + 1) = ddthetalist;
end
thetamat = thetamat'; %关节角度
dthetamat = dthetamat'; %关节角速度
ddthetamat = ddthetamat';%关节角加速度

subplot(231)
plot(thetamat(:,1))
title('第一关节角度变化')
subplot(232)
plot(thetamat(:,2))
title('第二关节角度变化')
subplot(233)
plot(thetamat(:,3))
title('第三关节角度变化')
subplot(234)
plot(thetamat(:,4))
title('第四关节角度变化')
subplot(235)
plot(thetamat(:,5))
title('第五关节角度变化')
subplot(236)
plot(thetamat(:,6))
title('第六关节角度变化')


figure
subplot(231)
plot(dthetamat(:,1))
title('第一关节角速度变化')
subplot(232)
plot(dthetamat(:,2))
title('第二关节角速度变化')
subplot(233)
plot(dthetamat(:,3))
title('第三关节角速度变化')
subplot(234)
plot(dthetamat(:,4))
title('第四关节角速度变化')
subplot(235)
plot(dthetamat(:,5))
title('第五关节角速度变化')
subplot(236)
plot(dthetamat(:,6))
title('第六关节角速度变化')

