function [out_alpha, out_AAdiag, out_w, out_qr,Control] = attitude_calculation(W0, QR0, alpha, t1, t2)
%������ ��� "���+�����������"

%������� ��������� ������ � �������� ����������
% clear all
% close all
% clc

%�������� ����������
definitions

%�������� ���������� ����������
control_parameters

%����������� ������� ���������� ������ �����
%������ ��-���������. ��������� ��������: @controlSM - ����������
%����������, @contro�Damp - ������������� �������� ��������
control=@controlPD;

%����������: ��������� ������ ��� ������������� � ��-����������,
%������������ ��� ����������� ����������
global mm
%����������� ��������
global w02
%��� ��������������
global h
h=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��������� ������ - ������������ ���������� � ������� �� ����� � ������ ��������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% qr=angle2quat(60.0*DEG_TO_RAD,130.0*DEG_TO_RAD,100.0*DEG_TO_RAD,'YZX')';
% qr=angle2quat(alpha0(2),alpha0(3),alpha0(1),'YZX')';
qr = QR0;
A=quat2dcm(qr');
%��������� �������� � ��������
% w=[10;10;10]*w02;
w=W0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����� ����� ��������������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=3*3600/h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��������� �������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������ �����������
ww=zeros(3,time+1);
AAdiag=zeros(3,time+1);
% AAdiag_kal=zeros(3,time+1);
tt=zeros(1,time);
QR = zeros(4,time+1);
%��
B_measurement=zeros(3,time); %��������� ���������� ����
% Est=zeros(7,time);           %������ ������� ���������
Est=zeros(10,time);           %������ ������� ���������
% State=zeros(7,time);         %�������� �������� ������� ���������
State=zeros(10,time);         %�������� �������� ������� ���������
dw=0.1*pi/180*ones(3,1); %�������� ���� ����
%������ �������� ���������� ����������
nom=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������� ������ ��������� � ����������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
steps_est=1; %����� ����� ��� ����������� �������� ��������
steps_ctrl=5;%����� ����� ��� ����������
steps_both=steps_est+steps_ctrl;%��������� ����������

%������ ��������� ������ � ��������� �������
ww(:,1)=w;
AAdiag(1,1)=A(1,1);
AAdiag(2,1)=A(2,2);
AAdiag(3,1)=A(3,3);
% alpha_new = zeros(t2-t1+1);
alpha_new(1) = acos(A(1,1));
tt(1)=0;
mm=zeros(3,1);
%������������� ���������� ������� �������
Kalman_init; 
%�������� ���� ��� ������������� ������ ���
K = [cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1]; %����� -> �����
% K = [-sin(alpha) cos(alpha) 0; 0 0 1; cos(alpha) sin(alpha) 0]; %����� -> �����
for i=1:t2-t1+1  %i=t1:t2
% i
%     A(:,i)=alpha0';
    %������� ����� � �������� ������
%     t=h*(i-1)+tt(1);
    t=h*(t1+i-1)+tt(1);
    u =  w02*t + u0;
    
    %������ ���������� ���� � ���
    Bmag=Bmag_orb(t,u);
    %������ ���������� ���� � ���
    Bx = A*Bmag;

    %��������� ��� ����������, � ����������� �� �������� ����. ����������
    %�������� �� ������ ������ ������� �������
    if (rem(i,steps_both)<=steps_est)&&(rem(i,steps_both)>0)
        mm_last=mm;
        mm=zeros(3,1);
    else
%         A1 = A_est*K';
%         mm = control(w_est,Bx_est,A1);
        A1 = A*K';
        mm = control(w,Bx,A1);
    end
    
    if max(abs(mm))>0.3
        mm=0.3*mm/max(abs(mm));
    end

    %��������� ���������� �� ��������� �������� ��������
%     mm = control(w,Bx,A);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������������� ��������� �������� ������� �����-����� 4-�� �������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [qr,w] = RK4(w,qr,u,t);
    A=quat2dcm(qr'); 
    
    %������ ������� 
    if (rem(i,steps_both)<=steps_est)&&(rem(i,steps_both)>0)&&(i>1)
        %���������� ���������� ���� � ��� - ������������� ���������
        %������� ������������ ������������ 100���+���������� ������ ���
        %200��� � 100��� �������� ���� ���������
%         B_measurement(:,nom)=normrnd(Bx+b_shift*1e-9*ones(3,1),100*1e-9*ones(3,1));
        B_measurement(:,nom)=normrnd(Bx+b_shift*1e-9*[1;0;0],100*1e-9*ones(3,1));
        b_sigma=100/norm (Bx)*1e-9;
%         Kalman_parameters{2,1}=diag([ones(3,1)*b_sigma^2;ones(3,1)*w_sigma^2]);
        Kalman_parameters{2,1}=diag([ones(3,1)*b_sigma^2]);
        W_measurement(:,nom)=normrnd(w+dw,w_sigma*ones(3,1));
        %��������� ��������� ������� ��������� �������� - ����� ����������
        %��� ���, ���� ��� ���� �������� ���������
        if (rem(i,steps_both)==1)
            last_interval=steps_ctrl;
        else
            last_interval=1;
        end
        %���������� ���������� � ��������� ������
        State(:,nom+1)=[qr;w;dw];
        %�������������� ������������ �������, ����� �������� �������������
        %� ��
        mm=mm_last;  
        %������ ������� 
%         [Est(:,nom+1), P] = Kalman_mag(Est(:,nom), P, last_interval, B_measurement(:,nom),Bmag,Kalman_parameters,u,t);
%         [Est(:,nom+1), P] = Kalman_mag_dus(Est(:,nom), P, last_interval, [B_measurement(:,nom);W_measurement(:,nom)],Bmag,Kalman_parameters,u,t);
%         [Est(:,nom+1), P] = Kalman_mag_dus_dw(Est(:,nom), P, last_interval, [B_measurement(:,nom);W_measurement(:,nom)],Bmag,Kalman_parameters,u,t);
        mm=zeros(3,1);
        %��������� ������ ������� ��������, ��������, ���� ��� ����������
        A_est=quat2dcm(Est(1:4,nom+1)');
        w_est=Est(5:7,nom+1);
        %������ ������� ���������� ���� � ��� �� ������� ������� �������
        Bmag=Bmag_orb(t+h,u+w02*h);
        Bx_est=quatrotate(Est(1:4,nom+1)',Bmag')';
        nom=nom+1;
    end
    
    %���������� ���������� ��� ���������� ��������
    ww(:,i+1)=w;
    AAdiag(1,i+1)=A(1,1);
    AAdiag(2,i+1)=A(2,2);
    AAdiag(3,i+1)=A(3,3);
    
    alpha_new(i+1) = acos(A(1,1));
    QR(:,i+1)=qr;
    Control(:,i)=mm;
    
    tt(:,i+1)=t;
end

out_alpha = alpha_new(1:t2-t1+2);
out_AAdiag = AAdiag(:,1:t2-t1+2);

out_qr = QR(:,1:t2-t1+2);
out_w = ww(:,1:t2-t1+2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure('Color',[1 1 1])
% subplot(2,1,1);
% plot(tt/3600,ww(1,:)*180/pi,'-k',tt/3600,ww(2,:)*180/pi,':k',tt/3600,ww(3,:)*180/pi,'-.k');
% legend('\omega_1','\omega_2','\omega_3','Location','Best')
% xlabel('time, hours');
% ylabel('\Omega, deg/s');
% 
% subplot(2,1,2);
% plot(tt/3600,acos(AAdiag(1,:))*180/pi,'-k',tt/3600,acos(AAdiag(2,:))*180/pi,':k',tt/3600,acos(AAdiag(3,:))*180/pi,'-.k');
% legend('\gamma_1_1','\gamma_2_2','\gamma_3_3','Location','Best')
% xlabel('time, hours');
% ylabel('direction cosines, deg');
% ylim([0 180])
% % title('��������� �������� ��')
% 
% figure('Color',[1 1 1])
% tt1=1:nom;
% tt1=tt1*steps_both;
% %����� ��������, �� ������� ������������ ���������� ��������������������
% %������
% N=800;
% %���������� ������ ������ ������� ������� ��������
% Err_Omega1=(State(5,1:nom)'-Est(5,1:nom)')*180/pi;
% Err_Omega2=(State(6,1:nom)'-Est(6,1:nom)')*180/pi;
% Err_Omega3=(State(7,1:nom)'-Est(7,1:nom)')*180/pi;
% subplot(2,1,1);
% plot(tt1/3600,Err_Omega1,'-k',tt1/3600,Err_Omega2,':k',tt1/3600,Err_Omega3,'-.k');
% hold on
% %���������� �������������������� ����������
% sigma1=sqrt(sum(Err_Omega1(end-N:end).^2)/N);
% sigma2=sqrt(sum(Err_Omega2(end-N:end).^2)/N);
% sigma3=sqrt(sum(Err_Omega3(end-N:end).^2)/N);
% s='\Delta\omega_1; \sigma';
% leg1=sprintf('%s=%6.2e deg/s',s,sigma1);
% s='\Delta\omega_2; \sigma';
% leg2=sprintf('%s=%6.2e deg/s',s,sigma2);
% s='\Delta\omega_3; \sigma';
% leg3=sprintf('%s=%6.2e deg/s',s,sigma3);
% legend(leg1,leg2,leg3)
% xlabel('time, hours');
% ylabel('\Delta\Omega, deg/s');
% ylim([-1e-2 1e-2])
% % title('������ ����������� ������� ��������')
% 
% % figure('Color',[1 1 1])
% %���������� ����� ���������� �� �����������
% [psi phi theta]=quat2angle(Est(1:4,1:nom)','YXZ');
% [psi_true phi_true theta_true]=quat2angle(State(1:4,1:nom)','YXZ');
% %���������� ������� ����� ��������� ������ � �� ��������
% delt1=(psi-psi_true)*180/pi;
% delt2=(phi-phi_true)*180/pi;
% delt3=(theta-theta_true)*180/pi;
% subplot(2,1,2);
% plot(tt1/3600,delt1,'-k',tt1/3600,delt2,':k',tt1/3600,delt3,'-.k');
% %���������� �������������������� ����������
% sigma1=sqrt(sum(delt1(end-N:end).^2)/N);
% sigma2=sqrt(sum(delt2(end-N:end).^2)/N);
% sigma3=sqrt(sum(delt3(end-N:end).^2)/N);
% s='\Delta\alpha; \sigma';
% leg1=sprintf('%s=%6.2e deg',s,sigma1);
% s='\Delta\beta; \sigma';
% leg2=sprintf('%s=%6.2e deg',s,sigma2);
% s='\Delta\gamma; \sigma';
% leg3=sprintf('%s=%6.2e deg',s,sigma3);
% legend(leg1,leg2,leg3)
% xlabel('time, hours');
% ylabel('Euler angles estimation error, deg');
% ylim([-5 5])
end
% % title('������ ����������� ����������')

% figure('Color',[1 1 1])
% plot(tt1/3600,(Est(8,1:nom)-dw(1)*ones(1,nom))*180/pi,'-k',tt1/3600,(Est(9,1:nom)-dw(2)*ones(1,nom))*180/pi,':k',tt1/3600,(Est(10,1:nom)-dw(1)*ones(1,nom))*180/pi,'-.k');
% ylim([-0.1 0.1])
% xlabel('�����, ����');
% ylabel('������  \Delta\omega, ^\circ/c');