clc; clear;
tspan=[0:0.1:30];
t=rem(tspan,4);

%% Desired State Quantities

alpha=flip([-0.293572873100118;-2.79480142732306e-13;-6.31088724176809e-30;1.56911520813689e-14;-2.04498001865826;3.11179352991225;-1.95854046380799;0.652725763443537;-0.121486339735189;0.0119709556256578;-0.000487732209207617]);
th1_des = polyval(alpha,t);
alpha_d = polyder(alpha);
th1dot_des = polyval(alpha_d,t);
alpha_dd = polyder(alpha_d);
th1ddot_des = polyval(alpha_dd,t);
 
beta=flip([1.35868301631567;4.21588689884326e-13;-1.26217744835362e-29;-4.23365046723726e-14;3.10970978903363;-4.40264764951139;2.32812051434579;-0.548250042474179;0.0438464762609300;0.00313436879623240;-0.000511556707107722]);
th2_des = polyval(beta,t);
beta_d = polyder(beta);
th2dot_des = polyval(beta_d,t);
beta_dd = polyder(beta_d);
th2ddot_des = polyval(beta_dd,t);

%% Parameters
p1 = 3.4*ones(length(tspan),1);
p2 = 0.4*ones(length(tspan),1);
p3 = 0.3*ones(length(tspan),1);
v1 = 0.4*ones(length(tspan),1);
v2 = 0.2*ones(length(tspan),1);
c1 = 0.2*ones(length(tspan),1);
c2 = 0.1*ones(length(tspan),1);
%% Solving the ODE
% State Vector = [th1;th2;th1dot;th2dot;p1;p2;p3;v1;v2;c1;c2]
% intial condition
q0=[-0.28;1.358;0.01;-0.08;3.2;0.35;0.22;0.34;0.19;0.21;0.09];


[t,q]=ode45(@(t,q)adaptfriccontroller(t,q,alpha,beta),tspan,q0);

error = q - [th1_des' th2_des' th1dot_des' th2dot_des' p1 p2 p3 v1 v2 c1 c2];
lbd=5;
for i=1:length(t)
    p1 = q(i,5); p2 = q(i,6); p3 = q(i,7); v1 = q(i,8); v2 = q(i,9); c1 = q(i,10); c2 = q(i,11);
M_hat=[p1+2*p3*cos(q(i,2)) p2+p3*cos(q(i,2)); p2+p3*cos(q(i,2)) p2];
C_hat=[-p3*q(i,4)*sin(q(i,2)), (-p3*sin(q(i,2))*(q(i,3) + q(i,4))); p3*q(i,3)*sin(q(i,2)), 0];
f_hat=[v1*q(i,3) + c1*sign(q(i,3)); v2*q(i,4) + c2*sign(q(i,4))];
K=[90 0;0 90];
qr_dot = [th1dot_des(i);th2dot_des(i)] - lbd*([q(i,1);q(i,2)] - [th1_des(i);th2_des(i)]);
qr_ddot = [th1ddot_des(i);th2ddot_des(i)] - lbd*([q(i,3);q(i,4)] - [th1dot_des(i);th2dot_des(i)]);
ep = [q(i,1);q(i,2)] - [th1_des(i);th2_des(i)];
epdot = [q(i,3);q(i,4)] - [th1dot_des(i);th2dot_des(i)];
s = epdot + lbd*ep;
tau(:,i)= M_hat*qr_ddot + C_hat*qr_dot + f_hat - K*s;
end

figure(1)
subplot(2,1,1);
plot(t,th1_des);
hold on;
plot(t,q(:,1));
legend('Theta1-Desired','Actual Theta1');
title("Theta1");
subplot(2,1,2);
plot(t,error(:,1));
title("Error in Theta1");

figure(2)
subplot(2,1,1);
plot(t,th2_des);
hold on;
plot(t,q(:,2));
legend('Theta2-Desired','Actual Theta2');
title("Theta2");
subplot(2,1,2);
plot(t,error(:,2));
title("Error in Theta2");

figure(3)
subplot(2,1,1);
plot(t,th1dot_des);
hold on;
plot(t,q(:,3));
legend('Theta1Dot-Desired','Actual Theta1Dot');
title("Theta1Dot");
subplot(2,1,2);
plot(t,error(:,3));
title("Error in Theta1Dot");

figure(4)
subplot(2,1,1);
plot(t,th2dot_des);
hold on;
plot(t,q(:,4));
legend('Theta2Dot-Desired','Actual Theta2Dot');
title("Theta2Dot");
subplot(2,1,2);
plot(t,error(:,4));
title("Error in Theta2Dot");

figure(5)
subplot(3,1,1);
plot(t,error(:,5));
title("Error in p1");
subplot(3,1,2);
plot(t,error(:,6));
title("Error in p2");
subplot(3,1,3);
plot(t,error(:,7));
title("Error in p3");

figure(6)
subplot(4,1,1);
plot(t,error(:,8));
title("Error in v1");
subplot(4,1,2);
plot(t,error(:,9));
title("Error in v2");
subplot(4,1,3);
plot(t,error(:,10));
title("Error in c1");
subplot(4,1,4);
plot(t,error(:,11));
title("Error in c2");

figure(7)
subplot(2,1,1);
plot(t,tau(1,:));
title("tau_1");
subplot(2,1,2);
plot(t,tau(2,:));
title("tau_2");