function[q] = Modern_Robotics_IKsolver61(q0,Td)
%�ִ�������ѧ����6.1�����˶�ѧ���
theta = q0;
for i = 1:10
   disp(i)
    T_sb = Modern_Robotics_FKsolver61(theta);
    T_bd = inv(T_sb)*Td;  %��eΪ��
    Vb_frame = logm(T_bd);
    v = Vb_frame(1:3,4);
    w_frame = Vb_frame(1:3,1:3);
    %�ٶ�����Vb ά��6*1 
    Vb = [w_frame(6);w_frame(7);w_frame(2);v];
    xitew = Vb(1:3);
    xitev = Vb(4:6);
    errorw = norm(xitew)
    errorv = norm(xitev)
    if errorw < 0.001 || errorv <0.0001
       q = theta;
       break
    end 
    %��������趨�������Сֵ��ֹͣѭ��
    %�������������޸�theta��ֵ
    %���㵱ǰ�Ƕȵ������ſ˱Ⱦ���
    J = Modern_Robotics_Jacobe61(theta);
    pinv(J)*Vb;
    theta = theta + pinv(J)*Vb;
    for j = 1:2
       while theta(j) < -pi
           theta(j) = theta(j) + pi;
       end
       while theta(j) > pi
           theta(j) = theta(j) - pi;
       end
    end
   theta*180/pi
end

end

