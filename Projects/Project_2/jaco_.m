clc;
clear;
Y1 = [1 1 1]';
for i = 1:5
y_int = Y1;%[1 1 1]';%Y1(1,:)';
syms y1 y2 y3;
y_var = [y1,y2,y3];
g_1 = y2 - 2*y1 + y1*y2 -y1 +1;
g_2 = y3 + y1 -2*y2 + y2*y3 -y1*y2;
g_3 = -1.8*y3 -y2*y3 + y2 + 0.2;
jaco = jacobian([g_1,g_2,g_3],[y1 y2 y3]);
jaco = double(subs(jaco,{y1, y2, y3},y_int'));
g1 = @(y_var) y_var(2) -2 * y_var(1) + y_var(1)*y_var(2) - y_var(1) + 1;
g2 = @(y_var) y_var(3) + y_var(1) -2*y_var(2) + y_var(2)*y_var(3) -y_var(1)*y_var(2);
g3 = @(y_var) -1.8*y_var(3) -y_var(2)*y_var(3) + y_var(2) + 0.2;
G = [g1(y_int');g2(y_int');g3(y_int')];
Y1 = y_int - inv(jaco) * G;
fprintf('%d iteration \n',i);
disp(Y1);
end
