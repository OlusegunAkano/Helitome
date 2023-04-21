%neo-hookean

function F = projectcode(x)
% options = odeset('RelTol',1e-12,'AbsTol',1e-12);
%P = x(1);
[theta,U] = ode45(@equation,[0 pi],[3.4 0 0 4.8235] );
ans = U(end,1)
hold on
grid on
figure(8)
plot(U(:,1),U(:,3))
xlabel('\rho')
ylabel('\eta')
title('\gamma = 0.3')
end 

function dydx = equation(theta,U)
    gamma = 0.6;
    P = 4.2627;
    R = 1 + (gamma* (cos(theta)));

    Sn_1 = (((U(2)^2 + U(4)^2)^2)*((U(1)^2)*(R^2))) + (3 *(gamma^4)*(R^4));
    Sn_2 = (((U(2)^2 + U(4)^2)^2)*((U(1)^2)*(R^2))) - ((gamma^4)*(R^4));

    Vn_1 = (2*((U(2)*R) + (U(1)*gamma*sin(theta))))*((gamma^4)*(R^5));
    Vn_2 = (((U(2)^2 + U(4)^2)^2)*((U(1)^2)*(R^2))) - ((gamma^4)*(R^4));
    Vn_3 = (((U(2)^2 + U(4)^2)^2)*((U(1)^4)*(gamma^2))) - ((gamma^4)*(R^4));

    lambda1 = 1/gamma * sqrt(U(2)^2 + U(4)^2);
    lambda2 = U(1)/R;

    V1 = (((Vn_2*U(1)*U(2)*((U(2)^2 + U(4)^2)^2)*gamma*R*sin(theta))+(Vn_3*((U(2)^2 + U(4)^2)^3)*(R^2)))*(lambda1^6)*(lambda2^2)*(gamma^4)*(R^4))-(Vn_1*U(2)*((U(2)^2 + U(4)^2)^2)*(lambda1^6)*(lambda2^2)*(gamma^4)*(R^4)) - (0.5*P*U(1)*U(4)*((U(2)^2 + U(4)^2)^(1.5))*(lambda1^11)*(lambda2^5)*(gamma^10)*(R^10));
    V2 = (((Vn_2*U(1)*U(2)*((U(4)^2 + U(4)^2)^2)*gamma*R*sin(theta)))*(lambda1^6)*(lambda2^2)*(gamma^4)*(R^4))-(Vn_1*U(4)*((U(2)^2 + U(4)^2)^2)*(lambda1^6)*(lambda2^2)*(gamma^4)*(R^4)) + (0.5*P*U(1)*U(2)*((U(2)^2 + U(4)^2)^(1.5))*(lambda1^11)*(lambda2^5)*(gamma^10)*(R^10));

    S22 = ((U(1)*(R^2)*(U(2)^2 + U(4)^2)*((Sn_1*U(2)^2) + (Sn_2*U(4)^2)))*(lambda1^6)*(lambda2^2)*(gamma^4)*(R^4));
    S24 = ((U(1)*(R^2)*(U(2)^2 + U(4)^2)*((Sn_1*U(2)*U(4)) - (Sn_2*U(2)*U(4))))*(lambda1^6)*(lambda2^2)*(gamma^4)*(R^4));
    S42 = ((U(1)*(R^2)*(U(2)^2 + U(4)^2)*((Sn_1*U(2)*U(4)) - (Sn_2*U(2)*U(4))))*(lambda1^6)*(lambda2^2)*(gamma^4)*(R^4));
    S44 = ((U(1)*(R^2)*(U(2)^2 + U(4)^2)*((Sn_1*U(4)^2) + (Sn_2*U(2)^2)))*(lambda1^6)*(lambda2^2)*(gamma^4)*(R^4));

            mass = [1 0 0 0; 
            0 S22 0 S24; 
            0 0 1 0;
            0 S42 0 S44];

    dydx = mass\ [U(2); V1; U(4); V2];
end