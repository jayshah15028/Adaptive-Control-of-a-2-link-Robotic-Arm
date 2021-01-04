function angles = inv_kine(x,y)

% L1 and L2 are the link lengths.
L1 = 0.38; L2 = 0.24;
th1 = zeros(2,1); th2 = zeros(2,1);
th2(1) = acos((x^2 + y^2 - L1^2 - L2^2)/(2*L1*L2));
th2(2) = -acos((x^2 + y^2 - L1^2 - L2^2)/(2*L1*L2));

for i=1:2
    A = [L1+L2*cos(th2(i)) -L2*sin(th2(i)); L2*sin(th2(i)) L1+L2*cos(th2(i))];
    a = inv(A)*[x;y];
    th1(i) = atan2(a(2),a(1));
end

angles = [th1 th2];

end