%mooney-rivlin
%initial guesses M gamma, rho, x - P etatheta


function F = projectcodemr(x)
% options = odeset('RelTol',1e-12,'AbsTol',1e-12);
%P = x(1);
[theta,U] = ode45(@equation,[0 pi],[3.1 0 0 1.1756 ]);
F = U(end,1)
hold on
grid on
figure(2)
plot(U(:,1),U(:,3))
xlabel('\rho')
ylabel('\eta')
title('\gamma = 0.4')
end 


function dydx = equation(theta,U)
        gamma = 0.4;
        M = 0.3; P = 5.8384;
        R = 1 + (gamma* (cos(theta)));

        Sm_1 = (((U(2)^2 + U(4)^2)^2)*((U(1)^2)*(R^2))) + (3 *((gamma)^4)*(R^4));%checked
        Sm_2 = (((U(2)^2 + U(4)^2)^2)*((U(1)^2)*(R^2))) - (((gamma)^4)*(R^4));%checked
        Sm_3 = ((-3*((U(2)^2 + U(4)^2)^(-2)))*((U(1)^(-2))*(R^(-2)))) - (((gamma)^(-4))*(R^(-4)));%checked
        Sm_4 = (((U(2)^2 + U(4)^2)^(-2))*((U(1)^(-2))*(R^(-2)))) - (((gamma)^(-4))*(R^(-4)));%checked

        Vm_1 = (2*((U(2)*R) + (U(1)*gamma*sin(theta))))*((gamma^4)*(R^5));%checked
        Vm_2 = (((U(2)^2 + U(4)^2)^2)*((U(1)^2)*(R^2))) - (((gamma)^4)*(R^4));%checked
        Vm_3 = ((-2)*((U(2)*R) + (U(1)*gamma*sin(theta))))*((gamma^(-4))*(R^(-3)));%checked
        Vm_4 = (((U(2)^2 + U(4)^2)^(-2))*((U(1)^(-2))*(R^(-2)))) - (((gamma)^(-4))*(R^(-4)));%checked
        Vm_5 = (((U(2)^2 + U(4)^2))*((U(1)^4)*(gamma^2))) - (((gamma)^4)*(R^4));%checked
        Vm_6 = (((U(2)^2 + U(4)^2)^(-1))*((U(1)^(-4))*(gamma^(-2)))) - (((gamma)^(-4))*(R^(-4)));%checked

        lambda1 = 1/gamma * sqrt(U(2)^2 + U(4)^2);
        lambda2 = U(1)/R;

        V1 = (((Vm_2*U(1)*U(2)*((U(2)^2 + U(4)^2)^2)*gamma*R*sin(theta))+(Vm_5*((U(2)^2 + U(4)^2)^3)*(R^2)))*(lambda1^4))-((Vm_1*U(2)*((U(2)^2 + U(4)^2)^2))*(lambda1^4))-M*((Vm_4*U(1)*U(2)*((U(2)^2 + U(4)^2)^2)*gamma*R*sin(theta)) + (Vm_6*((U(2)^2 + U(4)^2)^3)*(R^2)))*(lambda1^8)*(lambda2^4)*(gamma^8)*(R^8) + M*(Vm_3* U(2)*((U(2)^2 + U(4)^2)^2))*(lambda1^8)*(lambda2^4)*(gamma^8)*(R^8) - 0.5* P* U(1)*U(4)*((U(2)^2 + U(4)^2)^(1.5))*(lambda1^9)*(lambda2^3)*(gamma^6)*(R^6);%checked
        V2 = ((Vm_2*U(1)*U(4)*((U(2)^2 + U(4)^2)^2)*gamma*R*sin(theta))*(lambda1^4))-((Vm_1*U(4)*((U(2)^2 + U(4)^2)^2))*(lambda1^4))-(M*((Vm_4*U(1)*U(4)*((U(2)^2 + U(4)^2)^2)*gamma*R*sin(theta)))*(lambda1^8)*(lambda2^4)*(gamma^8)*(R^8)) +(M*(Vm_3*U(4)*((U(2)^2 + U(4)^2)^2))*(lambda1^8)*(lambda2^4)*(gamma^8)*(R^8)) + (0.5*P*U(1)*U(2)*((U(2)^2 + U(4)^2)^(1.5))*(lambda1^9)*(lambda2^3)*(gamma^6)*(R^6));%checked

        S22 = ((U(1)*(R^2)*(U(2)^2 + U(4)^2)*((Sm_1*(U(2)^2)) + (Sm_2*(U(4)^2))))*(lambda1^4)) - M*(U(1)*(R^2)*(U(2)^2 + U(4)^2)*((Sm_3*(U(2)^2)) + (Sm_4*(U(4)^2))))*(lambda1^8)*(lambda2^4)*(gamma^8)*(R^8);%checked
        S24 = ((U(1)*(R^2)*(U(2)^2 + U(4)^2)*((Sm_1*U(2)*U(4)) - (Sm_2*U(2)*U(4))))*(lambda1^4)) - (M*(U(1)*(R^2)*(U(2)^2 + U(4)^2)*((Sm_3*U(2)*U(4)) - (Sm_4*U(2)*U(4))))*(lambda1^8)*(lambda2^4)*(gamma^8)*(R^8));%checked
        S42 = S24;
        S44 = ((U(1)*(R^2)*(U(2)^2 + U(4)^2)*((Sm_1*U(4)^2) + (Sm_2*U(2)^2)))*(lambda1^4)) - (M*(U(1)*(R^2)*(U(2)^2 + U(4)^2)*((Sm_3*U(4)^2) + (Sm_4*U(2)^2)))*(lambda1^8)*(lambda2^4)*(gamma^8)*(R^8));%checked

        mass = [1 0 0 0; 
                0 S22 0 S24; 
                0 0 1 0;
                0 S42 0 S44];
        
        dydx = mass\ [U(2); V1; U(4); V2];
end