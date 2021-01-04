function dq=adaptcontroller(t,q,alpha,beta)

% State Vector = [th1,th2,th1dot,th2dot,p1,p2,p3]
%% Gains
K=[120 0;0 120]; 
gamma1 = 0.500; gamma2 = 0.150; gamma3 = 0.6;
%% Desired State Qunatities
time = rem(t,4);
th1_des = polyval(alpha,time);
alpha_d = polyder(alpha);
th1dot_des = polyval(alpha_d,time);
alpha_dd = polyder(alpha_d);
th1ddot_des = polyval(alpha_dd,time);
 

th2_des = polyval(beta,time);
beta_d = polyder(beta);
th2dot_des = polyval(beta_d,time);
beta_dd = polyder(beta_d);
th2ddot_des = polyval(beta_dd,time);

%% Parameters
 lbd = 3.5;
qr_dot = [th1dot_des;th2dot_des] - lbd*([q(1);q(2)] - [th1_des;th2_des]);
qr_ddot = [th1ddot_des;th2ddot_des] - lbd*([q(3);q(4)] - [th1dot_des;th2dot_des]);
ep = [q(1);q(2)] - [th1_des;th2_des];
epdot = [q(3);q(4)] - [th1dot_des;th2dot_des];
s = epdot + lbd*ep;

%% Determining Y0
p1 = 0; p2 = 0; p3 = 0;
M0=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C0=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];
Y0 = M0*qr_ddot + C0*qr_dot;
%% Determining Y1
p1 = 1; p2 = 0; p3 = 0;
M1=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C1=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];
Y1 = M1*qr_ddot + C1*qr_dot - Y0;

%% Determining Y2
p1 = 0; p2 = 1; p3 = 0;
M2=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C2=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];
Y2 = M2*qr_ddot + C2*qr_dot - Y0;

%% Determining Y3
p1 = 0; p2 = 0; p3 = 1;
M3=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C3=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];
Y3 = M3*qr_ddot + C3*qr_dot - Y0;

%% Defining M_hat and C_hat
p1=q(5); p2=q(6); p3=q(7);
M_hat=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C_hat=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];


%% Matrices and other necessary equations for solving ODE
p1 = 3.4; p2 = 0.4; p3 = 0.3;
M=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];

tau= M_hat*qr_ddot + C_hat*qr_dot - K*s;

%% Defining ODE equations
dq(1:2)=q(3:4);
dq(3:4)=inv(M)*(tau-C*q(3:4));
dq(5) = (-1/gamma1)*(s'*Y1);
dq(6) = (-1/gamma2)*(s'*Y2);
dq(7) = (-1/gamma3)*(s'*Y3);
dq=dq(:);

end