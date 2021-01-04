clear; clc;
% There are 2 joint variables: th1-->Theta1; th2-->Theta2;
%% Inverse Kinematics 
L1 = 0.38; L2 = 0.24;
angles1 = inv_kine(0.48,0.1);
angles2 = inv_kine(0.38,0);
angles3 = inv_kine(0.48,-0.1);
angles4 = inv_kine(0.58,0);

%% Trajectory Generation
% Hermite Polynomial Interpolation
% A = co-efficients of alphas; alpha=[alpha0;alpha1;...;alpha10]; b =
% values of conditions 

A = [1 0 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0 0;
     0 0 2 0 0 0 0 0 0 0 0;
     0 0 0 6 0 0 0 0 0 0 0;
     1 1 1 1 1 1 1 1 1 1 1;
     1 2 4 8 16 32 64 128 256 512 1024;
     1 3 9 27 81 243 729 2187 6561 19683 59049;
     1 4 16 64 256 1024 4096 16384 65536 262144 1048576;
     0 1 8 3*4^2 4*4^3 5*4^4 6*4^5 7*4^6 8*4^7 9*4^8 10*4^9;
     0 0 2 24 12*4^2 20*4^3 30*4^4 42*4^5 56*4^6 72*4^7 90*4^8;
     0 0 0 6 96 60*4^2 120*4^3 210*4^4 336*4^5 504*4^6 720*4^7];
 
 b1 = [angles1(1,1); 0; 0; 0; angles2(1,1); angles3(1,1); angles4(1,1); angles1(1,1); 0; 0; 0];
 
 alpha = flip(A\b1);
 
 tspan=[0:0.1:30]';
 t=rem(tspan,4);
 th1_des = polyval(alpha,t);
 alpha_d = polyder(alpha);
 th1dot_des = polyval(alpha_d,t);
 alpha_dd = polyder(alpha_d);
 th1ddot_des = polyval(alpha_dd,t);
 alpha_ddd = polyder(alpha_dd);
 j1_des = polyval(alpha_ddd,t);
 
 figure(1)
 subplot(4,1,1);
 plot(tspan,th1_des);
 title('Theta1-Desired');
 subplot(4,1,2);
 plot(tspan,th1dot_des);
 title('ThetaDot1-Desired');
 subplot(4,1,3);
 plot(tspan,th1ddot_des);
 title('ThetaDDot1-Desired');
 subplot(4,1,4);
 plot(tspan,j1_des);
 title('Jerk1-Desired');
 
 b2 = [angles1(1,2); 0; 0; 0; angles2(1,2); angles3(1,2); angles4(1,2); angles1(1,2); 0; 0; 0];
 
 beta = flip(A\b2);
 
 th2_des = polyval(beta,t);
 beta_d = polyder(beta);
 th2dot_des = polyval(beta_d,t);
 beta_dd = polyder(beta_d);
 th2ddot_des = polyval(beta_dd,t);
 beta_ddd = polyder(beta_dd);
 j2_des = polyval(beta_ddd,t);
 figure(2)
 subplot(4,1,1);
 plot(tspan,th2_des);
 title('Theta2-Desired');
 subplot(4,1,2);
 plot(tspan,th2dot_des);
 title('ThetaDot2-Desired');
 subplot(4,1,3);
 plot(tspan,th2ddot_des);
 title('ThetaDDot2-Desired');
 subplot(4,1,4);
 plot(tspan,j2_des);
 title('Jerk2-Desired');
 
%% Forward Kinematics for x vs. y plot of the links
% For link 1: x1 and for link 2: x2
x1 = [L1*cos(th1_des) L1*sin(th2_des)];
x2 = [L1*cos(th1_des) + L2*cos(th1_des + th2_des) , L1*sin(th1_des) + L2*sin(th1_des + th2_des)];

figure(3)
subplot(2,1,1);
plot(x1(:,1),x1(:,2)); 
title('Trajectory of 1st Link in Space');
subplot(2,1,2);
plot(x2(:,1),x2(:,2)); 
title('Trajectory of 2nd Link in Space');

 %% Torque Calculation
 
 p1 = 3.4; p2 = 0.4; p3 = 0.3;
 
 for i=1:length(t)
 M=[p1+2*p3.*cos(th2_des(i)) p2+p3.*cos(th2_des(i)); p2+p3.*cos(th2_des(i)) p2];
 C=[-p3.*th2dot_des(i).*sin(th2_des(i)) -p3.*sin(th2_des(i)).*(th1dot_des(i) + th2dot_des(i)); p3.*th1dot_des(i).*sin(th2_des(i)) 0];
 T(:,i)=M*[th1ddot_des(i);th2ddot_des(i)] + C*[th1dot_des(i);th2dot_des(i)];
 end
 
 figure(4)
 subplot(2,1,1);
 plot(tspan,T(1,:));
 title('Desired Torque for Joint 1');
 subplot(2,1,2);
 plot(tspan,T(2,:));
 title('Desired Torque for Joint 2');