function dq = pdcontroller(t,q,alpha,beta)

% State Vector = [th1;th2;th1dot;th2dot]
%% Gains
Kp=[350 0;0 300]; Kd=diag([40,50],0); 

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
p1 = 3.4; p2 = 0.4; p3 = 0.3;

%% Matrices from EOM
M=[p1+2*p3*cos(q(2)) p2+p3*cos(q(2)); p2+p3*cos(q(2)) p2];
C=[-p3*q(4)*sin(q(2)) -p3*sin(q(2))*(q(3) + q(4)); p3*q(3)*sin(q(2)) 0];
tau=M*[th1ddot_des;th2ddot_des] + C*[th1dot_des;th2dot_des] - Kp*([q(1);q(2)]-[th1_des;th2_des]) - Kd*([q(3);q(4)]-[th1dot_des;th2dot_des]);

%% Differential Equations
dq(1:2)=q(3:4);
dq(3:4)=inv(M)*(tau-C*q(3:4));
dq=dq(:);

end