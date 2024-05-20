ZL = 150;
f0 = 100e6;
Zo = 50;
gamma_yL  = z2gamma(1/ZL, 1/Zo);

beta = 2 * pi / (8); % 2π/λ * λ/8 = π/4
GZL=z2gamma(ZL, Zo);
GZL_t = GZL * exp(-1j * beta); %GZL_transformado


%clf;                                                                         %%clear figure
abaco           = figure();
rfckt_handler   = rfckt.passive;                                             %%elemento de rfckt a dibujar
[h1, hsm]       = circle (rfckt_handler, f0, 'Gamma', abs(gamma_yL),'R',1);
hold on;

[h1,hsm] = circle (rfckt_handler,f0,'R',1);

x1 = h1.XData-0.5;
y1 = h1.YData+0.5;

plot(x1,y1);


%%busco las intersecciones
%[ptA, ptB] = imped_match_find_circle_intersections_helper([0 0], abs(gamma_yL), [0.5 0], 0.5);
[ptA, ptB] = imped_match_find_circle_intersections_helper([0 0], abs(gamma_yL), [0.5 0.5], 0.5);
%plot(ptA(1),ptA(2), 'k.','MarkerSize',16);
plot(ptB(1),-ptB(2), 'k.','MarkerSize',16);

ptA_mod = norm(ptA);                ptB_mod = norm(ptB);                                %%Compute magnitudes
ptA_ang = atan2(ptA(2), ptA(1));    ptB_ang = atan2(ptB(2), ptB(1));                    %%Compute angles
gamma_A = ptA_mod*exp(1i*ptA_ang);  gamma_B = ptB_mod*exp(-1i*ptB_ang);                  %%With it compute gammas
%text (real(gamma_A) + .1, imag(gamma_A)-.1, 'P_A','FontSize', 12, 'FontWeight','Bold'); %%graphing points
text (real(gamma_B) + .1, imag(gamma_B)+.1, 'y_1','FontSize', 12, 'FontWeight','Bold'); 

%%grafico yL
plot(real(gamma_yL), imag(gamma_yL), 'k.','MarkerSize',16);
text (real(gamma_yL) + .05, imag(gamma_yL)-.08, 'y_L','FontSize', 12, 'FontWeight','Bold');

z_aux=gamma2z(gamma_B)/Zo;
[h1,hsm] = circle (rfckt_handler,f0,'R',real(z_aux));



[pt1, pt2] = imped_match_find_circle_intersections_helper([0 0.5], abs(gamma_yL), [1-1.25/2 0], 1.25/2);
%[pt1, pt2] = imped_match_find_circle_intersections_helper([0 0.5], abs(gamma_yL),[max(x1)-min(x1), max(y1)-min(y1)], 0.5);
plot(pt1(1),pt1(2), 'k.','MarkerSize',16);
plot(pt2(1),pt2(2), 'k.','MarkerSize',16);

pt1_mod = norm(pt1);                pt2_mod = norm(pt2);                                %%Compute magnitudes
pt1_ang = atan2(pt1(2), pt1(1));    pt2_ang = atan2(pt2(2), pt2(1));                    %%Compute angles
gamma_1 = pt1_mod*exp(1i*pt1_ang);  gamma_2 = pt2_mod*exp(1i*pt2_ang);                  %%With it compute gammas
text (real(gamma_1) + .1, imag(gamma_1)-.1, 'y_{2a}','FontSize', 12, 'FontWeight','Bold'); %%graphing points
text (real(gamma_2) + .1, imag(gamma_2)+.1, 'y_{2b}','FontSize', 12, 'FontWeight','Bold');


%Me quedo con el P1 porque tiene menos ROE. Ahora me tengo que mover a
%coeficiente de reflexión contante hasta el círculo unitario. entonces
%tengo que crear un círculo centrado en 0 de ratio P1. 


[h1, hsm]       = circle (rfckt_handler, f0, 'Gamma', abs(gamma_1)); %Circulo gamma constante de p1
%ahora tengo que intersectar con el circulo unitario 

[pt1a, pt2a] = imped_match_find_circle_intersections_helper([0 0], abs(gamma_1), [0.5 0], 0.5);
plot(pt1a(1),pt1a(2), 'k.','MarkerSize',16);
plot(pt2a(1),pt2a(2), 'k.','MarkerSize',16);

pt1a_mod = norm(pt1a);                pt2a_mod = norm(pt2a);                                %%Compute magnitudes
pt1a_ang = atan2(pt1a(2), pt1a(1));    pt2a_ang = atan2(pt2a(2), pt2a(1));                    %%Compute angles
gamma_1a = pt1a_mod*exp(1i*pt1a_ang);  gamma_2a = pt2a_mod*exp(1i*pt2a_ang);                  %%With it compute gammas
text (real(gamma_1a) + .1, imag(gamma_1a)-.1, 'y_{3a}','FontSize', 12, 'FontWeight','Bold'); %%graphing points
text (real(gamma_2a) + .1, imag(gamma_2a)+.1, 'y_{3b}','FontSize', 12, 'FontWeight','Bold');

%La admitancia y2a se mueve a gamma_2a constante hasta y3b (o y3a)
%stub_1=y2a-y1

y2a=gamma2z(gamma_1)/Zo;
y1=z_aux;

%stub_1=y2a-y1;

%el stub2 agrega una susceptancia de (-y3b)
y3a=(gamma2z(gamma_2a)/Zo); %me tengo que quedar con la parte imaginaria!


%ahora dibujo los stubs 

%y1=gamma2z(gamma_B)/Zo
ys1=y2a-y1; %stub 1
ys2=1-y3a;  %stub 2
gamma_s1=z2gamma(ys1*Zo,Zo);
gamma_s2=z2gamma(ys2*Zo,Zo);
[h1, hsm]       = circle (rfckt_handler, f0, 'X', imag(ys1)); %reactancia que agrega el stub1
plot(real(gamma_s1), imag(gamma_s1), 'k.','MarkerSize',16);
text (real(gamma_s1)-.1,  imag(gamma_s1)+.1, 'y_{S1}','FontSize', 12, 'FontWeight','Bold');

[h1, hsm]       = circle (rfckt_handler, f0, 'X', imag(ys2));
plot(real(gamma_s2), imag(gamma_s2), 'k.','MarkerSize',16);
text (real(gamma_s2)+.1,  imag(gamma_s2)-.1, 'y_{S2}','FontSize', 12, 'FontWeight','Bold');

%cálculo de lambdas: 
theta1=180-atan2(imag(gamma_s1),real(gamma_s1))*180/pi;
l_stub1=theta1*0.25/180
theta2=180-atan2(imag(gamma_s2),real(gamma_s2))*180/pi;
l_stub2=theta2*0.25/180


%[h1, hsm]       = circle (rfckt_handler, f0, 'X', imag(y3b));
