function [Bmag, dB] = Bmag_orb(t,u)

global incl F

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������� ������-������ (���������) ��������
%������ - ��������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_sat=[cos(u) cos(incl)*sin(u) sin(incl)*sin(u)]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������� �������� �� �������, ��������� � �������, � �����������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G=[cos(u) -sin(u) 0;
sin(u) cos(u) 0;
0 0 1;];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������ ��������� ���� - ��������� ������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bmag=Inclined_field(t,r_sat);

%������� ��� ����� - ��� ������ - ����������� - ���������
Bmag=G*F'*Bmag;
dB=Bmag;
Bmag=[Bmag(2);Bmag(3);Bmag(1)];
end

