function[q] = Inverse_kinematics_POE(q0,Td)
%�˳��������˶�ѧ��ֵ���,���õķ����ǻ��������ſɱȾ����ţ��-����ɭ��
%��ʼ������֪Ŀ�����Td����ʼ�Ƕȹ���ֵq0���趨i=0
%�����ʼ����ֵ����ʵֵû���㹻�ӽ�����������̿��ܲ�����
%������ĽǶ����ʼ�Ƕ��й�
theta = q0;
for i = 1:1000
   disp(i)
   T_sb = Forward_kinematics(theta);
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
   %��������趨�������Сֵ��ֹͣѭ��
   if errorw < 0.00001 || errorv <0.00001
       q = theta;
       break
   end
   %�������������޸�theta��ֵ
   %���㵱ǰ�Ƕȵ������ſ˱Ⱦ���
   J = Jacoby(theta);
   theta = theta + pinv(J)*Vb;
   for j = 1:6
       while theta(j) < -pi
           theta(j) = theta(j) + pi;
       end
       while theta(j) > pi
           theta(j) = theta(j) - pi;
       end
   end
end

end