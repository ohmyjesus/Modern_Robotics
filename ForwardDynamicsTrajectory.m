function [thetamat, dthetamat]  = ForwardDynamicsTrajectory(thetalist, dthetalist, taumat, g, Ftipmat, Mlist, Glist, Slist, dt, intRes)
% *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
% Takes thetalist: n维初始关节角,
%       dthetalist:n维初始关节角速度,
%       taumat: 力/力矩的N×N矩阵，其中每一行是任何时间步的合力
%       g: Gravity vector g,
%       Ftipmat: 末端执行器施加的N x 6空间力矩阵（如果没有尖端力，用户应输入零，将使用零矩阵）
%       Mlist: List of link frames {i} relative to {i-1} at the home
%              position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns,
%       dt: 连续合力/力矩之间的时间步长
%       intRes: 积分分辨率是在每个时间步之间进行积分（Euler）的次数。必须是大于或等于1的整数值.
% Returns thetamat: 由指定关节力/力矩产生的机器人关节角的n x n矩阵,
%         dthetamat: 由指定关节力/力矩产生的机器人关节角速度的n x n矩阵,
% 该函数模拟给定关节力/转矩开环历史的串行链的运动。它调用使用前向动力学的数值积分过程。
% Example Inputs (3 Link Robot):
% 
% clc; clear;
% thetalist = [0.1; 0.1; 0.1];
% dthetalist = [0.1; 0.2; 0.3];
% taumat = [[3.63, -6.58, -5.57]; [3.74, -5.55, -5.5]; ...
%         [4.31, -0.68, -5.19]; [5.18, 5.63, -4.31]; ...
%         [5.85, 8.17, -2.59]; [5.78, 2.79, -1.7]; ...
%         [4.99, -5.3, -1.19]; [4.08, -9.41, 0.07]; ...
%         [3.56, -10.1, 0.97]; [3.49, -9.41, 1.23]];
% %Initialise robot description (Example with 3 links)
% g = [0; 0; -9.8];
% Ftipmat = ones(size(taumat, 1), 6);
% M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
% M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
% M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
% M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
% G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7]);
% G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393]);
% G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
% Glist = cat(3, G1, G2, G3);
% Mlist = cat(3, M01, M12, M23, M34);  
% Slist = [[1; 0; 1;      0; 1;     0], ...
%        [0; 1; 0; -0.089; 0;     0], ...
%        [0; 1; 0; -0.089; 0; 0.425]];
% dt = 0.1;
% intRes = 8;
% [thetamat, dthetamat] ...
% = ForwardDynamicsTrajectory(thetalist, dthetalist, taumat, g, ...
%                           Ftipmat, Mlist, Glist, Slist, dt, intRes);
% %Output using matplotlib to plot the joint forces/torques
% Tf = size(taumat, 1);
% time=0: (Tf / size(thetamat, 1)): (Tf - (Tf / size(thetamat, 1)));
% plot(time,thetamat(:, 1),'b')
% hold on
% plot(time,thetamat(:, 2), 'g')
% plot(time,thetamat(:, 3), 'r')
% plot(time,dthetamat(:, 1), 'c')
% plot(time,dthetamat(:, 2), 'm')
% plot(time,dthetamat(:, 3), 'y')
% title('Plot of Joint Angles and Joint Velocities')
% xlabel('Time')
% ylabel('Joint Angles/Velocities')
% legend('Theta1', 'Theta2', 'Theta3', 'DTheta1', 'DTheta2', 'DTheta3')
%

taumat = taumat';
Ftipmat = Ftipmat';
thetamat = taumat;
thetamat(:, 1) = thetalist;
dthetamat = taumat;
dthetamat(:, 1) = dthetalist;

for i = 1: size(taumat, 2) - 1
    for j = 1: intRes
       ddthetalist = ForwardDynamics(thetalist, dthetalist, taumat(:,i), g,  Ftipmat(:, i), Mlist, Glist, Slist);     
       [thetalist, dthetalist] = EulerStep(thetalist, dthetalist, ddthetalist, dt / intRes);
    end
    thetamat(:, i + 1) = thetalist;
    dthetamat(:, i + 1) = dthetalist;
    
end
thetamat = thetamat';
dthetamat = dthetamat';

end