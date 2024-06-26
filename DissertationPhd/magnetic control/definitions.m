%������ �������
global J Jreal
%��������� ������
global w02 w0 incl F p
%��������� ������ �������
global t0scal
%��������� ������ ���������� ����
global B0 k0
%������� �������� � �������
global DEG_TO_RAD
%������������ ������ ��������
global h_wheel

DEG_TO_RAD=pi/180;

%�������
h_wheel=[0 ;0.0; 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������ �������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%"���������" ������ ������� �� �� - ������������ � ������ �������� � ��
J=zeros(3,3);
% J(1,1)=3.600;
% J(2,2)=5.836;
% J(3,3)=2.468;
J(1,1)=0.01;
J(2,2)=0.01;
J(3,3)=0.007;
%"��������" ������ ������� � ����������� 10 ��������� - ������������ � �������������
%������ ���������
% Jreal=zeros(3,3);
% Jreal(1,1)=5836;
% Jreal(2,2)=2468;
% Jreal(3,3)=3600;
%��������� �������, �� ���� �� �������� ������� ����������� ���������� ��
%����� ��� �� 10%
add=normrnd(zeros(3,1),ones(3,1));
add=.1*add/max(abs(add));
% Jreal=zeros(3,3);
% % Jreal(1,1)=(1+add(1))*J(1,1);
% % Jreal(2,2)=(1+add(2))*J(2,2);
% % Jreal(3,3)=(1+add(3))*J(3,3);
% Jreal(1,1)=(0.9)*J(1,1);
% % Jreal(2,2)=(0.9)*J(2,2);
% Jreal(3,3)=(0.9)*J(3,3);
% Jreal(2,2)=Jreal(3,3)+Jreal(1,1)-0.1;
%���������� �������������
% Jreal=normrnd(J,.1*J/3)
%������ �������� �����
Jreal=J;

% if (Jreal(1,1)+Jreal(2,2)<Jreal(3,3))||(Jreal(3,3)+Jreal(2,2)<Jreal(1,1))||(Jreal(3,3)+Jreal(1,1)<Jreal(2,2))
%     break
% end 
% Jreal
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��������� ������
%������ ������ � ������� ������ �����
% H=750;
H=340;
% H=8000;
r_earth=6371;
r = r_earth+H;
%����������
incl = 60*pi/180;
%������� ����������� ���� � �������� ����������
Omega0=0;
omega=0;
%��������� �������� ������
u0=0;
%�������� ������
p=r;
%����������� �������� ������������ �� ������
w02=(398600.4415/r^3)^.5;
w0=[0;w02;0];
%������� �������� �� ������������, ��������� � �������,
%� ������������, ��������� � ������
F(1,1)=cos(Omega0)*cos(omega)-sin(Omega0)*sin(omega)*cos(incl);
F(1,2)=-cos(Omega0)*sin(omega)-sin(Omega0)*cos(omega)*cos(incl);
F(1,3)=sin(Omega0)*sin(incl);
F(2,1)=sin(Omega0)*cos(omega)+cos(Omega0)*sin(omega)*cos(incl);
F(2,2)=-sin(Omega0)*sin(omega)+cos(Omega0)*cos(omega)*cos(incl);
F(2,3)=-cos(Omega0)*sin(incl);
F(3,1)=sin(omega)*sin(incl);
F(3,2)=cos(omega)*sin(incl);
F(3,3)=cos(incl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��������� ���������� ������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%"��������"
mu=7.812e6;
B0=mu/p^3;
%���� ���������� ������ ������������ �����
lambda0=-71.88*DEG_TO_RAD;
delta0=10.26*DEG_TO_RAD;
%��������� ������ ������ ������������ �����
k0(1)=cos(lambda0)*sin(delta0);
k0(2)=sin(lambda0)*sin(delta0);
k0(3)=cos(delta0);

%%%%%%%%%%%%%%%%%%%%%%%
%��������� ������ �������
t0 = [2015 1 1 0 0 0];
%��������� � ������
t0scal=datenum(t0);
%%%%%%%%%%%%%%%%%%%%%%%