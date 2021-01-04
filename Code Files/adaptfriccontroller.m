function dq=adaptfriccontroller(t,q,alpha,beta)

% State Vector = [th1,th2,th1dot,th2dot,p1,p2,p3,v1,v2,c1,c2]
%% Gains
K=[90 0;0 90]; 
gamma1 = 0.500; gamma2 = 0.150; gamma3 = 0.6; gamma4 = 1.7; gamma5 = 1.4; gamma6 = 2.5; gamma7 = 2; 
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
p1 = 0; p2 = 0; p3 = 0; v1 = 0; v2 = 0; c1 = 0; c2 = 0;
M0=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C0=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];
f0=[v1*q(3) + c1*sign(q(3)); v2*q(4) + c2*sign(q(4))];
Y0 = M0*qr_ddot + C0*qr_dot + f0;
%% Determining Y1
p1 = 1; p2 = 0; p3 = 0; v1 = 0; v2 = 0; c1 = 0; c2 = 0;
M1=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C1=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];
f1=[v1*q(3) + c1*sign(q(3)); v2*q(4) + c2*sign(q(4))];
Y1 = M1*qr_ddot + C1*qr_dot + f1 - Y0;

%% Determining Y2
p1 = 0; p2 = 1; p3 = 0; v1 = 0; v2 = 0; c1 = 0; c2 = 0;
M2=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C2=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];
f2=[v1*q(3) + c1*sign(q(3)); v2*q(4) + c2*sign(q(4))];
Y2 = M2*qr_ddot + C2*qr_dot + f2 - Y0;

%% Determining Y3
p1 = 0; p2 = 0; p3 = 1; v1 = 0; v2 = 0; c1 = 0; c2 = 0;
M3=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C3=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];
f3=[v1*q(3) + c1*sign(q(3)); v2*q(4) + c2*sign(q(4))];
Y3 = M3*qr_ddot + C3*qr_dot + f3 - Y0;

%% Determining Y4
p1 = 0; p2 = 0; p3 = 1; v1 = 1; v2 = 0; c1 = 0; c2 = 0;
M4=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C4=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];
f4=[v1*q(3) + c1*sign(q(3)); v2*q(4) + c2*sign(q(4))];
Y4 = M4*qr_ddot + C4*qr_dot + f4 - Y0;

%% Determining Y5
p1 = 0; p2 = 0; p3 = 1; v1 = 0; v2 = 1; c1 = 0; c2 = 0;
M5=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C5=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];
f5=[v1*q(3) + c1*sign(q(3)); v2*q(4) + c2*sign(q(4))];
Y5 = M5*qr_ddot + C5*qr_dot + f5 - Y0;

%% Determining Y6
p1 = 0; p2 = 0; p3 = 1; v1 = 0; v2 = 0; c1 = 1; c2 = 0;
M6=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C6=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];
f6=[v1*q(3) + c1*sign(q(3)); v2*q(4) + c2*sign(q(4))];
Y6 = M6*qr_ddot + C6*qr_dot + f6 - Y0;
%% Determining Y7
p1 = 0; p2 = 0; p3 = 1; v1 = 0; v2 = 0; c1 = 0; c2 = 1;
M7=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C7=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];
f7=[v1*q(3) + c1*sign(q(3)); v2*q(4) + c2*sign(q(4))];
Y7 = M7*qr_ddot + C7*qr_dot + f7 - Y0;

%% Defining M_hat and C_hat
p1=q(5); p2=q(6); p3=q(7); v1=q(8); v2=q(9); c1=q(10); c2=q(11);
M_hat=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C_hat=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];
f_hat=[v1*q(3) + c1*sign(q(3)); v2*q(4) + c2*sign(q(4))];

%% Matrices and other necessary equations for solving ODE
p1 = 3.4; p2 = 0.4; p3 = 0.3; v1 = 0.4; v2 = 0.2; c1 = 0.2; c2 = 0.1;
M=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];
f=[v1*q(3) + c1*sign(q(3)); v2*q(4) + c2*sign(q(4))];

tau= M_hat*qr_ddot + C_hat*qr_dot + f_hat - K*s;

%% Defining ODE equations
dq(1:2)=q(3:4);
dq(3:4)=inv(M)*(tau-C*q(3:4)-f);
dq(5) = (-1/gamma1)*(s'*Y1);
dq(6) = (-1/gamma2)*(s'*Y2);
dq(7) = (-1/gamma3)*(s'*Y3);
dq(8) = (-1/gamma4)*(s'*Y4);
dq(9) = (-1/gamma5)*(s'*Y5);
dq(10) = (-1/gamma6)*(s'*Y6);
dq(11) = (-1/gamma7)*(s'*Y7);
dq=dq(:);

end