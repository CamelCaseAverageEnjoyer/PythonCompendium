function [wn, qsatN]=rightSide(w,qsatP,u_cur, t_cur)

global Jreal w02 w0 mm h_wheel

%���������� ������� ������������ ���������
A=quat2dcm(qsatP');

%��������� ����
Bmag = Bmag_orb(t_cur,u_cur);
Bx = A*Bmag;

%����������� �� �������� ���������� �������
m0=10;
%����������� �� �������� ������������� �������
M0=1*m0*norm(Bx);

%������������� ��� ��-���������
Mm=cross(mm,Bx);

%����������� �������� ����������
if max(abs(Mm))>M0
    Mm=M0*Mm/max(abs(Mm));
end
%��� ������������� ������������� ������� �� ��������
% Mm=M0*Mm/max(abs(Mm));

%�������������� ����
r_sat_bound=A*[0;0;1];
Mgr=3*w02^2*cross(r_sat_bound,Jreal*r_sat_bound);
% Mgr=zeros(3,1);
    
%���������� ��������� ����������
%������� �����������
Mdist=normrnd(zeros(3,1),ones(3,1));
Mdist=5e-7*Mdist/norm(Mdist);
%�� ��������� �� ������ 1e-4, ��������� �� -3 �� 5 e-4
% Mdist=normrnd(zeros(3,1),ones(3,1));
% Mdist=4e-7*Mdist/norm(Mdist)+1e-7*ones(3,1)/3^.5;
%���������� ������������ ��������
% Mdist=5e-7*ones(3,1)/3^.5;
% Mdist=zeros(3,1);

%����������
qsatN(1,1)=-0.5*w'*qsatP(2:4);
qsatN(2,1)=0.5*(w(1)*qsatP(1)+w(3)*qsatP(3)-w(2)*qsatP(4));
qsatN(3,1)=0.5*(w(2)*qsatP(1)-w(3)*qsatP(2)+w(1)*qsatP(4));
qsatN(4,1)=0.5*(w(3)*qsatP(1)+w(2)*qsatP(2)-w(1)*qsatP(3));

%����������� �������
w0n=A*w0;
W=[0,w(3),-w(2);-w(3),0,w(1);w(2),-w(1),0;];
wn=-inv(Jreal)*(cross(w+w0n,h_wheel)-Mm-Mgr-Mdist+cross(w,Jreal*w0n)+cross(w,Jreal*w)+cross(w0n,Jreal*w0n)+cross(w0n,Jreal*w)+Jreal*W*w0n);