function [J] = Jacoby_POE(q)
%���³����������ſɱȾ�������6*6����
%�ο��ִ�������ѧ5.1��,�ſɱȾ����������˶�ѧ��ָ�������
%����ȡ�ռ��ſɱȾ�����ת��Ϊ�����ſɱ�
zer = zeros(3,3);
%��ʼλ��ʱ��ĩ��λ��M����ڻ�����ϵ
M = [1 0 0 -0.0261;
     0 1 0 -1.1863;  
     0 0 1 -0.6540; 
     0 0 0 1];
w1 = [0;0;1]; %�عؽ�������ĵ�λ���� ����ڻ�����ϵ
w2 = [0;-1;0];
w3 = [0;-1;0];
w4 = [-1;0;0];
w5 = [0;1;0];
w6 = [0;0;1];

q1 = [0;0;0]; %qn�ǹؽ�������һ�� ����ֵ�ڻ�����ϵ�н��ж���
q2 = [-0.1573;0.0653;0] ;
q3 = [-0.2037;-0.5483;-0.0165] ;
q4 = [-0.1518;-0.7154;-0.3140] ;
q5 = [-0.1081;-0.6714;-0.6685] ;
q6 = [-0.0261;-1.1863;-0.6540] ;

v1 = -cross(w1,q1); %���ٶ�
v2 = -cross(w2,q2);
v3 = -cross(w3,q3);
v4 = -cross(w4,q4);
v5 = -cross(w5,q5);
v6 = -cross(w6,q6);

R1 = [0 -w1(3) w1(2);w1(3) 0 -w1(1);-w1(2) w1(1) 0]; %ת��ķ��Գƾ���
R2 = [0 -w2(3) w2(2);w2(3) 0 -w2(1);-w2(2) w2(1) 0];
R3 = [0 -w3(3) w3(2);w3(3) 0 -w3(1);-w3(2) w3(1) 0];
R4 = [0 -w4(3) w4(2);w4(3) 0 -w4(1);-w4(2) w4(1) 0];
R5 = [0 -w5(3) w5(2);w5(3) 0 -w5(1);-w5(2) w5(1) 0];
R6 = [0 -w6(3) w6(2);w6(3) 0 -w6(1);-w6(2) w6(1) 0];

G1 = eye(3)*q(1) + (1-cos(q(1)))*R1 + (q(1)-sin(q(1)))*R1^2;  
G2 = eye(3)*q(2) + (1-cos(q(2)))*R2 + (q(2)-sin(q(2)))*R2^2;
G3 = eye(3)*q(3) + (1-cos(q(3)))*R3 + (q(3)-sin(q(3)))*R3^2;
G4 = eye(3)*q(4) + (1-cos(q(4)))*R4 + (q(4)-sin(q(4)))*R4^2;
G5 = eye(3)*q(5) + (1-cos(q(5)))*R5 + (q(5)-sin(q(5)))*R5^2;
G6 = eye(3)*q(6) + (1-cos(q(6)))*R6 + (q(6)-sin(q(6)))*R6^2;

Rot1 = eye(3) + sin(q(1))*R1 + (1-cos(q(1)))*R1^2;  %����ת���ľ���ָ��
Rot2 = eye(3) + sin(q(2))*R2 + (1-cos(q(2)))*R2^2; 
Rot3 = eye(3) + sin(q(3))*R3 + (1-cos(q(3)))*R3^2;
Rot4 = eye(3) + sin(q(4))*R4 + (1-cos(q(4)))*R4^2; 
Rot5 = eye(3) + sin(q(5))*R5 + (1-cos(q(5)))*R5^2;
Rot6 = eye(3) + sin(q(6))*R6 + (1-cos(q(6)))*R6^2; 

eS1 = [Rot1 G1*v1 ;0 0 0 1];%�����˶��ľ���ָ��
eS2 = [Rot2 G2*v2 ;0 0 0 1];
eS3 = [Rot3 G3*v3 ;0 0 0 1];
eS4 = [Rot4 G4*v4 ;0 0 0 1];
eS5 = [Rot5 G5*v5 ;0 0 0 1];
eS6 = [Rot6 G6*v6 ;0 0 0 1];

Tsb = eS1*eS2*eS3*eS4*eS5*eS6*M; %���ڻ�����ϵ�ı�ʾ
%���������˶�ѧ��ָ������ʽ

%�����ǿռ��ſ˱Ⱦ�������
e_S1 = eS1;
e_S2 = eS1*eS2;
e_S3 = eS1*eS2*eS3;
e_S4 = eS1*eS2*eS3*eS4;
e_S5 = eS1*eS2*eS3*eS4*eS5;

e_S1_R = e_S1(1:3,1:3);
e_S2_R = e_S2(1:3,1:3);
e_S3_R = e_S3(1:3,1:3);
e_S4_R = e_S4(1:3,1:3);
e_S5_R = e_S5(1:3,1:3);

e_S1_P = e_S1(1:3,4);
e_S2_P = e_S2(1:3,4);
e_S3_P = e_S3(1:3,4);
e_S4_P = e_S4(1:3,4);
e_S5_P = e_S5(1:3,4);

e_S1_P_frame = [0 -e_S1_P(3) e_S1_P(2);e_S1_P(3) 0 -e_S1_P(1);-e_S1_P(2) e_S1_P(1) 0];
e_S2_P_frame = [0 -e_S2_P(3) e_S2_P(2);e_S2_P(3) 0 -e_S2_P(1);-e_S2_P(2) e_S2_P(1) 0];
e_S3_P_frame = [0 -e_S3_P(3) e_S3_P(2);e_S3_P(3) 0 -e_S3_P(1);-e_S3_P(2) e_S3_P(1) 0];
e_S4_P_frame = [0 -e_S4_P(3) e_S4_P(2);e_S4_P(3) 0 -e_S4_P(1);-e_S4_P(2) e_S4_P(1) 0];
e_S5_P_frame = [0 -e_S5_P(3) e_S5_P(2);e_S5_P(3) 0 -e_S5_P(1);-e_S5_P(2) e_S5_P(1) 0];

S1 = [w1;v1];
S2 = [w2;v2];
S3 = [w3;v3];
S4 = [w4;v4];
S5 = [w5;v5];
S6 = [w6;v6];

%�������6*6
Ad_es1 = [e_S1_R  zer;e_S1_P_frame*e_S1_R e_S1_R];
Ad_es1es2 = [e_S2_R  zer;e_S2_P_frame*e_S2_R e_S2_R];
Ad_es1es2es3 = [e_S3_R  zer;e_S3_P_frame*e_S3_R e_S3_R];
Ad_es1es2es3es4 = [e_S4_R  zer;e_S4_P_frame*e_S4_R e_S4_R];
Ad_es1es2es3es4es5 = [e_S5_R  zer;e_S5_P_frame*e_S5_R e_S5_R];

Js1 = S1;
Js2 = Ad_es1*S2;
Js3 = Ad_es1es2*S3;
Js4 = Ad_es1es2es3*S4;
Js5 = Ad_es1es2es3es4*S5;
Js6 = Ad_es1es2es3es4es5*S6;

%�ռ��ſ˱ȵ����
Js = [Js1 Js2 Js3 Js4 Js5 Js6] ;%6*6
%�����ſ˱ȵ����
Tbs = inv(Tsb);
R_Tbs = Tbs(1:3,1:3);
P_Tbs = Tbs(1:3,4);
P_Tbs_frame = [0 -P_Tbs(3) P_Tbs(2);P_Tbs(3) 0 -P_Tbs(1);-P_Tbs(2) P_Tbs(1) 0];
Ad_Tbs = [R_Tbs zer;P_Tbs_frame*R_Tbs R_Tbs];
Jb = Ad_Tbs * Js;


J = Jb;

end
