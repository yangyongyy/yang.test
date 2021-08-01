function [proj,scale,translation,freq,phase]=select_best(signal_r,N,a_base,j_min,j_max,u_base,p_min,v_base,k_min,w_base,i_min,i_max); 

% this subroutine is to select in the dictionary the best atom suited the siganl or the residual of the signal 

% INPUT 
% the signal_r: the signal or the residual of the signal to be decomposed 
% the N: the longth of the signal or of the residual of the signal or the length of the atoms 

% ����
% the signal_r: ʣ��Ĵ��ֽ��ź�
% N���źţ�ʣ���źţ�ԭ�ӵĳ���

% parameters :the parameter to construct the dictionary , it have much influence on the speed of the decomposition 
% a_base=2 
% j_min=0; 
% j_max=log2(N); 
 
%OUTPUT 
% proj: the projection of the signal or the residual of the signal on the best atom 
% the scale: the scale of the best atom (s in the formula) 
% the translation : the translation of the best atom (u in the formula) 
% the freq: the frequency of the best atom (v in the formula) 
% phase: the phase of the best atom (w in the formula) 

% ���
% proj: �����ԭ�ӻ��ϵ�ʣ���ź�
% the scale: ���ԭ�ӻ������� 
% the translation : ���ԭ�ӻ���ƽ��
% the freq: ���ԭ�ӻ���Ƶ��
% phase: ���ԭ�ӻ�����λ

% proj_trans :to determine which projection is biggest 
% proj_trans ���ж��ĸ��������õ�

proj_trans=0; 
proj=0; 
  
for j=j_min:j_max 
   for p=p_min:N*2^(-j+1) 
      for k=k_min:2^(j+1) 
         for i=i_min:i_max 
            s=a_base^j; % ����������ֵ��2Ϊ���� 
            u=p*s*u_base; % ����ƽ�Ƶ�ֵ
            v=k*(1/s)*v_base; % ����Ƶ��ֵ
            w=i*w_base; % ������λֵ
            t=0:N-1;
            t=(t-u)/s; % ��������ƽ�Ʊ任
             
            g=(1/sqrt(s))*exp(-pi*t.*t).*cos(v*t+w); % �Ե���ԭ���źŽ�������ƽ�ƺ͵���
             
            g=g/sqrt(sum(g.*g)); % ����һ��
             
            proj_trans=sum(signal_r.*g); % �����ؽ��ź�
             
            if abs(proj_trans)>abs(proj) % �ж��ؽ��ź��Ƿ�����ѣ����������Ӧ����ֵ������Ӧ������
               proj=proj_trans; 
               scale=s; 
               translation=u; 
               freq=v; 
               phase=w; 
            end 
         end 
      end 
   end 
end 