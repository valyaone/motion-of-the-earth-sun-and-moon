clear all;
%Пока всё в СИ
%g=6,67*10^(-11);
%Msun=1.9891*10^30; % масса Солнца кг
gmsun=13.267297*10^(19);% в CИ
vx0z=0; % начальная скорость земли по x
vy0z=29270; % скорость в перигелии  по оси y м/c
x0z=1.52098232*10^11; % расстояние Земли в апогеи м
y0z=0;
% MEarth=5.97*10^24;
gmearth=39.8199*10^13;
vx0m=0; % начальная скорость Луны по x
vy0m=29270+1023; % скорость в  по оси y м/c
x0m=1.52098232*10^11+3.84467*10^8; % расстояние Луны м
y0m=0;
Mmoon=7.35*10^22;
gmmoon=49.0245*10^11;

% tau=400;
tau=1000;
M=28*24*60*60/tau;
%M=10000

%ДАННЫЕ ЗЕМЛИ
U(1,1)=x0z;
U(1,2)=y0z;
U(1,3)=vy0z;
U(1,4)=vx0z;

%ДАННЫЕ ЛУНЫ
U(1,5)=x0m;
U(1,6)=y0m;
U(1,7)=vy0m;
U(1,8)=vx0m;

U;
figure(1);
method=1;

if method==1 %ERK1
    for m=2:M
        znam=((U(m-1,1)-U(m-1,5))^2+(U(m-1,2)-U(m-1,6))^2)^(3/2);
        w(m-1,1)=U(m-1,4);
        w(m-1,2)=U(m-1,3);
        w(m-1,3)=-gmsun*U(m-1,2)/((U(m-1,1)^2+U(m-1,2)^2)^(3/2))-gmmoon*(U(m-1,2)-U(m-1,6))/znam;
        w(m-1,4)=-gmsun*U(m-1,1)/((U(m-1,1)^2+U(m-1,2)^2)^(3/2))-gmmoon*(U(m-1,1)-U(m-1,5))/znam;
        w(m-1,5)=U(m-1,8);
        w(m-1,6)=U(m-1,7);
        w(m-1,7)=-gmsun*U(m-1,6)/((U(m-1,5)^2+U(m-1,6)^2)^(3/2))+gmearth*(U(m-1,2)-U(m-1,6))/znam;
        w(m-1,8)=-gmsun*U(m-1,5)/((U(m-1,5)^2+U(m-1,6)^2)^(3/2))+gmearth*(U(m-1,1)-U(m-1,5))/znam;
        for i=1:8
            U(m,i)=U(m-1,i)+tau*w(m-1,i);
        end
    end
end

if method==2 %ERK2
    a2=2/3;
    b1=1/4;
    b2=3/4;
    c1=0;
    c2=0;
    for m=2:M
        znam=((U(m-1,1)-U(m-1,5))^2+(U(m-1,2)-U(m-1,6))^2)^(3/2);
        w1(m-1,1)=U(m-1,4);
        w1(m-1,2)=U(m-1,3);
        w1(m-1,3)=-gmsun*U(m-1,2)/((U(m-1,1)^2+U(m-1,2)^2)^(3/2))-gmmoon*(U(m-1,2)-U(m-1,6))/znam;
        w1(m-1,4)=-gmsun*U(m-1,1)/((U(m-1,1)^2+U(m-1,2)^2)^(3/2))-gmmoon*(U(m-1,1)-U(m-1,5))/znam;
        w1(m-1,5)=U(m-1,8);
        w1(m-1,6)=U(m-1,7);
        w1(m-1,7)=-gmsun*U(m-1,6)/((U(m-1,5)^2+U(m-1,6)^2)^(3/2))+gmearth*(U(m-1,2)-U(m-1,6))/znam;
        w1(m-1,8)=-gmsun*U(m-1,5)/((U(m-1,5)^2+U(m-1,6)^2)^(3/2))+gmearth*(U(m-1,1)-U(m-1,5))/znam;
        
        Uwn=U+tau*a2*w1;
        
        znam=((Uwn(m-1,1)-Uwn(m-1,5))^2+(Uwn(m-1,2)-Uwn(m-1,6))^2)^(3/2);
        w2(m-1,1)=Uwn(m-1,4);
        w2(m-1,2)=Uwn(m-1,3);
        w2(m-1,3)=-gmsun*Uwn(m-1,2)/((Uwn(m-1,1)^2+Uwn(m-1,2)^2)^(3/2))-gmmoon*(Uwn(m-1,2)-Uwn(m-1,6))/znam;
        w2(m-1,4)=-gmsun*Uwn(m-1,1)/((Uwn(m-1,1)^2+Uwn(m-1,2)^2)^(3/2))-gmmoon*(Uwn(m-1,1)-Uwn(m-1,5))/znam;
        w2(m-1,5)=Uwn(m-1,8);
        w2(m-1,6)=Uwn(m-1,7);
        w2(m-1,7)=-gmsun*Uwn(m-1,6)/((Uwn(m-1,5)^2+Uwn(m-1,6)^2)^(3/2))+gmearth*(Uwn(m-1,2)-Uwn(m-1,6))/znam;
        w2(m-1,8)=-gmsun*Uwn(m-1,5)/((Uwn(m-1,5)^2+Uwn(m-1,6)^2)^(3/2))+gmearth*(Uwn(m-1,1)-Uwn(m-1,5))/znam;
        for i=1:8
            U(m,i)=U(m-1,i)+tau*(w1(m-1,i)*b1+b2*w2(m-1,i));
        end
    end
end

if method==3 %ERK3
    a2=1/2;
    a3=3/4; 
    b1=2/9;
    b2=3/9;
    b3=4/9;  
    for m=2:M
        znam=((U(m-1,1)-U(m-1,5))^2+(U(m-1,2)-U(m-1,6))^2)^(3/2);
        w1(m-1,1)=U(m-1,4);
        w1(m-1,2)=U(m-1,3);
        w1(m-1,3)=-gmsun*U(m-1,2)/((U(m-1,1)^2+U(m-1,2)^2)^(3/2))-gmmoon*(U(m-1,2)-U(m-1,6))/znam;
        w1(m-1,4)=-gmsun*U(m-1,1)/((U(m-1,1)^2+U(m-1,2)^2)^(3/2))-gmmoon*(U(m-1,1)-U(m-1,5))/znam;
        w1(m-1,5)=U(m-1,8);
        w1(m-1,6)=U(m-1,7);
        w1(m-1,7)=-gmsun*U(m-1,6)/((U(m-1,5)^2+U(m-1,6)^2)^(3/2))+gmearth*(U(m-1,2)-U(m-1,6))/znam;
        w1(m-1,8)=-gmsun*U(m-1,5)/((U(m-1,5)^2+U(m-1,6)^2)^(3/2))+gmearth*(U(m-1,1)-U(m-1,5))/znam;
        
        Uwn=U+tau*a2*w1;
        
        znam=((Uwn(m-1,1)-Uwn(m-1,5))^2+(Uwn(m-1,2)-Uwn(m-1,6))^2)^(3/2);
        w2(m-1,1)=Uwn(m-1,4);
        w2(m-1,2)=Uwn(m-1,3);
        w2(m-1,3)=-gmsun*Uwn(m-1,2)/((Uwn(m-1,1)^2+Uwn(m-1,2)^2)^(3/2))-gmmoon*(Uwn(m-1,2)-Uwn(m-1,6))/znam;
        w2(m-1,4)=-gmsun*Uwn(m-1,1)/((Uwn(m-1,1)^2+Uwn(m-1,2)^2)^(3/2))-gmmoon*(Uwn(m-1,1)-Uwn(m-1,5))/znam;
        w2(m-1,5)=Uwn(m-1,8);
        w2(m-1,6)=Uwn(m-1,7);
        w2(m-1,7)=-gmsun*Uwn(m-1,6)/((Uwn(m-1,5)^2+Uwn(m-1,6)^2)^(3/2))+gmearth*(Uwn(m-1,2)-Uwn(m-1,6))/znam;
        w2(m-1,8)=-gmsun*Uwn(m-1,5)/((Uwn(m-1,5)^2+Uwn(m-1,6)^2)^(3/2))+gmearth*(Uwn(m-1,1)-Uwn(m-1,5))/znam;
        
        Uwn3=U+tau*a3*w2;
        
        znam=((Uwn3(m-1,1)-Uwn3(m-1,5))^2+(Uwn3(m-1,2)-Uwn3(m-1,6))^2)^(3/2);
        w3(m-1,1)=Uwn3(m-1,4);
        w3(m-1,2)=Uwn3(m-1,3);
        w3(m-1,3)=-gmsun*Uwn3(m-1,2)/((Uwn3(m-1,1)^2+Uwn3(m-1,2)^2)^(3/2))-gmmoon*(Uwn3(m-1,2)-Uwn3(m-1,6))/znam;
        w3(m-1,4)=-gmsun*Uwn3(m-1,1)/((Uwn3(m-1,1)^2+Uwn3(m-1,2)^2)^(3/2))-gmmoon*(Uwn3(m-1,1)-Uwn3(m-1,5))/znam;
        w3(m-1,5)=Uwn3(m-1,8);
        w3(m-1,6)=Uwn3(m-1,7);
        w3(m-1,7)=-gmsun*Uwn3(m-1,6)/((Uwn3(m-1,5)^2+Uwn3(m-1,6)^2)^(3/2))+gmearth*(Uwn3(m-1,2)-Uwn3(m-1,6))/znam;
        w3(m-1,8)=-gmsun*Uwn3(m-1,5)/((Uwn3(m-1,5)^2+Uwn3(m-1,6)^2)^(3/2))+gmearth*(Uwn3(m-1,1)-Uwn3(m-1,5))/znam;
        for i=1:8
            U(m,i)=U(m-1,i)+tau*(w1(m-1,i)*b1+b2*w2(m-1,i)+b3*w3(m-1,i));
        end
    end
end

if method==4 %ERK4
    a2=1/2;
    a3=1/2; 
    a4=1;
    b1=1/6;
    b2=1/3;
    b3=1/3; 
    b4=1/6;
    for m=2:M
        znam=((U(m-1,1)-U(m-1,5))^2+(U(m-1,2)-U(m-1,6))^2)^(3/2);
        w1(m-1,1)=U(m-1,4);
        w1(m-1,2)=U(m-1,3);
        w1(m-1,3)=-gmsun*U(m-1,2)/((U(m-1,1)^2+U(m-1,2)^2)^(3/2))-gmmoon*(U(m-1,2)-U(m-1,6))/znam;
        w1(m-1,4)=-gmsun*U(m-1,1)/((U(m-1,1)^2+U(m-1,2)^2)^(3/2))-gmmoon*(U(m-1,1)-U(m-1,5))/znam;
        w1(m-1,5)=U(m-1,8);
        w1(m-1,6)=U(m-1,7);
        w1(m-1,7)=-gmsun*U(m-1,6)/((U(m-1,5)^2+U(m-1,6)^2)^(3/2))+gmearth*(U(m-1,2)-U(m-1,6))/znam;
        w1(m-1,8)=-gmsun*U(m-1,5)/((U(m-1,5)^2+U(m-1,6)^2)^(3/2))+gmearth*(U(m-1,1)-U(m-1,5))/znam;
        
        Uwn2=U+tau*a2*w1;
        
        znam=((Uwn2(m-1,1)-Uwn2(m-1,5))^2+(Uwn2(m-1,2)-Uwn2(m-1,6))^2)^(3/2);
        w2(m-1,1)=Uwn2(m-1,4);
        w2(m-1,2)=Uwn2(m-1,3);
        w2(m-1,3)=-gmsun*Uwn2(m-1,2)/((Uwn2(m-1,1)^2+Uwn2(m-1,2)^2)^(3/2))-gmmoon*(Uwn2(m-1,2)-Uwn2(m-1,6))/znam;
        w2(m-1,4)=-gmsun*Uwn2(m-1,1)/((Uwn2(m-1,1)^2+Uwn2(m-1,2)^2)^(3/2))-gmmoon*(Uwn2(m-1,1)-Uwn2(m-1,5))/znam;
        w2(m-1,5)=Uwn2(m-1,8);
        w2(m-1,6)=Uwn2(m-1,7);
        w2(m-1,7)=-gmsun*Uwn2(m-1,6)/((Uwn2(m-1,5)^2+Uwn2(m-1,6)^2)^(3/2))+gmearth*(Uwn2(m-1,2)-Uwn2(m-1,6))/znam;
        w2(m-1,8)=-gmsun*Uwn2(m-1,5)/((Uwn2(m-1,5)^2+Uwn2(m-1,6)^2)^(3/2))+gmearth*(Uwn2(m-1,1)-Uwn2(m-1,5))/znam;
                
        Uwn3=U+tau*a3*w2;
        
        znam=((Uwn3(m-1,1)-Uwn3(m-1,5))^2+(Uwn3(m-1,2)-Uwn3(m-1,6))^2)^(3/2);
        w3(m-1,1)=Uwn3(m-1,4);
        w3(m-1,2)=Uwn3(m-1,3);
        w3(m-1,3)=-gmsun*Uwn3(m-1,2)/((Uwn3(m-1,1)^2+Uwn3(m-1,2)^2)^(3/2))-gmmoon*(Uwn3(m-1,2)-Uwn3(m-1,6))/znam;
        w3(m-1,4)=-gmsun*Uwn3(m-1,1)/((Uwn3(m-1,1)^2+Uwn3(m-1,2)^2)^(3/2))-gmmoon*(Uwn3(m-1,1)-Uwn3(m-1,5))/znam;
        w3(m-1,5)=Uwn3(m-1,8);
        w3(m-1,6)=Uwn3(m-1,7);
        w3(m-1,7)=-gmsun*Uwn3(m-1,6)/((Uwn3(m-1,5)^2+Uwn3(m-1,6)^2)^(3/2))+gmearth*(Uwn3(m-1,2)-Uwn3(m-1,6))/znam;
        w3(m-1,8)=-gmsun*Uwn3(m-1,5)/((Uwn3(m-1,5)^2+Uwn3(m-1,6)^2)^(3/2))+gmearth*(Uwn3(m-1,1)-Uwn3(m-1,5))/znam;
        
        Uwn4=U+tau*a4*w3;
        
        znam=((Uwn4(m-1,1)-Uwn4(m-1,5))^2+(Uwn4(m-1,2)-Uwn4(m-1,6))^2)^(3/2);
        w4(m-1,1)=Uwn4(m-1,4);
        w4(m-1,2)=Uwn4(m-1,3);
        w4(m-1,3)=-gmsun*Uwn4(m-1,2)/((Uwn4(m-1,1)^2+Uwn4(m-1,2)^2)^(3/2))-gmmoon*(Uwn4(m-1,2)-Uwn4(m-1,6))/znam;
        w4(m-1,4)=-gmsun*Uwn4(m-1,1)/((Uwn4(m-1,1)^2+Uwn4(m-1,2)^2)^(3/2))-gmmoon*(Uwn4(m-1,1)-Uwn4(m-1,5))/znam;
        w4(m-1,5)=Uwn4(m-1,8);
        w4(m-1,6)=Uwn4(m-1,7);
        w4(m-1,7)=-gmsun*Uwn4(m-1,6)/((Uwn4(m-1,5)^2+Uwn4(m-1,6)^2)^(3/2))+gmearth*(Uwn4(m-1,2)-Uwn4(m-1,6))/znam;
        w4(m-1,8)=-gmsun*Uwn4(m-1,5)/((Uwn4(m-1,5)^2+Uwn4(m-1,6)^2)^(3/2))+gmearth*(Uwn4(m-1,1)-Uwn4(m-1,5))/znam;
        for i=1:4
            U(m,i)=U(m-1,i)+tau*(w1(m-1,i)*b1+b2*w2(m-1,i)+b3*w3(m-1,i)+b4*w4(m-1,i));
        end
    end
    %plot(U(:,2),U(:,1));
end

moon = plot(U(:,5)-U(:,1),U(:,6)-U(:,2),'k');
hold on;
earth1 = plot(0,0,'.');
xlim([-4.2 * 10^8 4.2 * 10^8]);
ylim([-4.2 * 10^8 4.2 * 10^8]);
xlabel('X');
ylabel('Y');
title('Вращение Луны вокруг Земли');
legend('Луна','Земля');

figure(2)
plot(U(:,1),U(:,2));
hold on;
sun = plot(0,0,'y.');
xlim([-1.6 * 10^11 1.6 * 10^11]);
ylim([-1.6 * 10^11 1.6 * 10^11]);
xlabel('X');
ylabel('Y');
title('Вращение Земли вокруг Солнца');
legend('Земля','Солнце');

