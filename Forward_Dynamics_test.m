%���³����ǻ�е��������ѧ����
%��ʼ����֮ǰ��Ҫ�õ��ؽ�����tau�����ݲ���ִ��
thetalist = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1];
dthetalist = [0.1; 0.2; 0.3; 0.4; 0.5; 0.6];
ddthetalist = [0 ; 0 ; 0 ; 0 ; 0 ; 0];
taumat = tau.Data;

g = [0; 0; -9.8];
Ftipmat = zeros(size(taumat, 1), 6);

dt = 0.1;
intRes = 8;    %���ֱַ�������ÿ��ʱ�䲽֮����л��֣�Euler���Ĵ���

[m,n] = size(taumat);
disp('����ִ��ʱ��Լ��Ҫ ')
time = m * intRes / 50 

M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
M45 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
M56 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
M67 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
Mlist = cat(3, M01, M12, M23, M34, M45, M56, M67); 

%����i�Ŀռ���Ծ���
G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7]);
G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393]);
G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
G4 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
G5 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
G6 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
Glist = cat(3, G1, G2, G3, G4, G5, G6);

w1 = [0;0;1]; %�عؽ���ĵ�λ���� ����ڻ�����ϵ
w2 = [0;-1;0];
w3 = [0;-1;0];
w4 = [1;0;0];
w5 = [0;1;0];
w6 = [0;0;-1];

q1 = [0; 0; 0]; %qn�ǹؽ�������һ�� ����ֵ�ڻ�����ϵ�н��ж��� 
q2 = [169.9876; 0; -0.0121] ;
q3 = [170.0136; 0; 613.9879] ;
q4 = [0; -0.0112; 814.0441] ;
q5 = [810.0109; 0; 814.0039] ;
q6 = [810.0018; 0; 0];

v1 = -cross(w1,q1); %���ٶ�
v2 = -cross(w2,q2);
v3 = -cross(w3,q3);
v4 = -cross(w4,q4);
v5 = -cross(w5,q5);
v6 = -cross(w6,q6);

%�˶�����S1-S6
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
thetamat = thetamat'; %�ؽڽǶ�
dthetamat = dthetamat'; %�ؽڽ��ٶ�
ddthetamat = ddthetamat';%�ؽڽǼ��ٶ�

subplot(231)
plot(thetamat(:,1))
title('��һ�ؽڽǶȱ仯')
subplot(232)
plot(thetamat(:,2))
title('�ڶ��ؽڽǶȱ仯')
subplot(233)
plot(thetamat(:,3))
title('�����ؽڽǶȱ仯')
subplot(234)
plot(thetamat(:,4))
title('���ĹؽڽǶȱ仯')
subplot(235)
plot(thetamat(:,5))
title('����ؽڽǶȱ仯')
subplot(236)
plot(thetamat(:,6))
title('�����ؽڽǶȱ仯')


figure
subplot(231)
plot(dthetamat(:,1))
title('��һ�ؽڽ��ٶȱ仯')
subplot(232)
plot(dthetamat(:,2))
title('�ڶ��ؽڽ��ٶȱ仯')
subplot(233)
plot(dthetamat(:,3))
title('�����ؽڽ��ٶȱ仯')
subplot(234)
plot(dthetamat(:,4))
title('���Ĺؽڽ��ٶȱ仯')
subplot(235)
plot(dthetamat(:,5))
title('����ؽڽ��ٶȱ仯')
subplot(236)
plot(dthetamat(:,6))
title('�����ؽڽ��ٶȱ仯')

