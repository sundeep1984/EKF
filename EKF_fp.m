% fixed point kalman filtering exercise - uses the fixed point toolbox


clear all
close all
tk=pi/4;
wk=1; % in rad/s
num_iter=5000;
xk_hat = fi([cos(tk) sin(tk) wk])';

yk_hat = fi([cos(tk) sin(tk)])';
Tsc=16e-3;%s
wb = 1;

H = fi([1 0 0; 0 1 0]);

P = fi(eye(3));
Q = fi(0.0001*eye(3));
R = fi(0.0001*eye(2));

outputs = [];
x_1temp=[0;0;0];
for i=1:num_iter
    
    x_1temp(1) = xk_hat.data(1)*round(cos(wb*Tsc*xk_hat.data(3))*2^16)/2^16 - xk_hat.data(2)*round(sin(wb*Tsc*xk_hat.data(3))*2^16)/2^16;
    x_1temp(2) = xk_hat.data(2)*round(cos(wb*Tsc*xk_hat.data(3))*2^16)/2^16 + xk_hat.data(2)*round(sin(wb*Tsc*xk_hat.data(3))*2^16)/2^16; 
    x_1temp(3) = xk_hat.data(3);
    
    x_1 = fi(x_1temp);
    %x_1 = [xk_hat(1)*cos(wb*Tsc*round(xk_hat.data(3)*2^16)/2^16)-xk_hat(2)*sin(wb*Tsc*round(xk_hat.data(3)*2^16)/2^16);xk_hat(1)*cos(wb*Tsc*round(xk_hat.data(3)*2^16)/2^16)+xk_hat(2)*sin(wb*Tsc*round(xk_hat.data(3)*2^16)/2^16); xk_hat(3)];
    %x_1=x_1+n_x;
    
    y_1=fi([1 0 0; 0 1 0])*x_1;
    
    F = fi([cos(wb*Tsc*xk_hat.data(3))  -sin(wb*Tsc*xk_hat.data(3))  -wb*Tsc*(xk_hat.data(1)*cos(wb*Tsc*xk_hat.data(3))+xk_hat.data(2)*sin(wb*Tsc*xk_hat.data(3)));
         cos(wb*Tsc*xk_hat.data(3))  sin(wb*Tsc*xk_hat.data(3))  -wb*Tsc*(xk_hat.data(1)*cos(wb*Tsc*xk_hat.data(3))-xk_hat.data(2)*sin(wb*Tsc*xk_hat.data(3)));
         0 0 1]);
     
     P_t = F.data*P.data*F.data' + Q.data;
     
     K = fi(P_t*H.data'*inv(H.data*P_t*H.data' + R.data));
     xk_hat = x_1 +K*(yk_hat-y_1);
     
     P = P_t - K*H*P_t;
     outputs = [outputs y_1];
end

out=sqrt(outputs.data(1,:).^2+outputs.data(2,:).^2);

figure
plot((1:num_iter)*Tsc,outputs.data(1,:))
hold on
plot((1:num_iter)*Tsc,outputs.data(2,:),'r')
ylabel('X and Y components');
xlabel('Time(s)')
figure
plot((1:num_iter)*Tsc,out);
xlabel('Time(s)')
