%% Modelos de Crecimiento con Agentes Heterogeneos
%% (El Modelo de Ayagari 1994)
% Mauricio Tejada
% 
% ILADES - Universidad Alberto Hurtado
% 
% Noviembre 2018
%% El Modelo
% Las familias resuelven el siguiente problema:
% 
% $$\max E_{0}\sum_{t=0}^{\infty}\beta^{t}\frac{c_{t}^{1-\sigma}}{1-\sigma}$$
% 
% sujetos a la restricción presupuestaria:
% 
% $$k_{t+1}=\begin{cases}\left[1+(1-\tau)(q_{t}-\delta)\right]k_{t}+(1-\tau)w_{t}-c_{t} 
% & si\,\epsilon_{t}=1\\\left[1+(1-\tau)(q_{t}-\delta)\right]k_{t}+s_{t}-c_{t} 
% & si\,\epsilon_{t}=0\end{cases}$$
% 
% donde $\epsilon = 1$corresponde a empleo y $\epsilon=0$ corresponde a desempleo. 
% El empleo evolucona en el tiempo según: 
% 
% $$\pi(\epsilon'|\epsilon)=\Pr\left[\epsilon_{t+1}=\epsilon'|\epsilon_{t}=\epsilon\right]=\left[\begin{array}{cc}p_{uu} 
% & p_{ue}\\p_{eu} & p_{ee}\end{array}\right]$$
% 
% La representación recursiva del problema es:
% 
% $$\begin{array}{rll}v(k,\epsilon;K) & =\max\left\{ u(c)+\beta E\left[v(k',\epsilon';K)|\epsilon\right]\right\} 
% \\s.a\\ & k'+c=\left[1+(1-\tau)(q(K)-\delta)\right]k+(1-\tau)\epsilon w(K)+(1-\epsilon)s_t\\ 
% & k\geq0\\ & \pi(\epsilon'|\epsilon)=\left[\begin{array}{cc}p_{uu} & p_{ue}\\p_{eu} 
% & p_{ee}\end{array}\right]\\ & k_{0},\epsilon_{0}\,\,\,dados.\end{array}$$
% 
% Las firmas resuleven el siguiente problema:
% 
% $$\max K_{t}^{\alpha}L_{t}^{1-\alpha}-w_{t}(K_{t})L_{t}-q_{t}(K_{t})K_{t}$$
% 
% Finalmente, el gobierno tiene su presupuesto equilibrado:
% 
% $$(1-L_{t})s_{t}=\tau\left[L_{t}w(K_{t})+K_{t}\left(q(K_{t})-\delta\right)\right]$$
%% Solución del Modelo
% La parametrización del modelo es la siguiente:

alpha = 0.36;
beta  = 0.97;
delta = 0.005;
sigma = 1.5;
rep   = 0.25;
pp    = [0.5,    0.5;
         0.0435, 0.9565];
     
     
% Opciones globales
flag_print_conv = 0;
     
% Opciones para iterar la función valor:     
N_points_valueiteration = 100;
kmax = 200;
saved_init = 1;
crit = 1e-5;
step = 1;
maxiter = 10000;

% Opciones para iterar la distribución invariante:
N_points_distinv = 3*N_points_valueiteration;
step_dist = 0.5;

% Opciones para iterar el capital agregado
crit_k = 0.1;
step_k = 0.01;
count_k = 1;
Kiter = [];

%% Paso 1:
% Calcular el empleo en estado estacionario $L_{ss}$ usando: $L_{t}=p_{ue}(1-L_{t-1})+p_{ee}L_{t-1}\longrightarrow   
% L_{ss}=\frac{p_{ue}}{1+p_{ue}-p_{ee}}$

L_ss = pp(1,2)/(1+pp(1,2)-pp(2,2));
U_ss = 1 - L_ss;
%disp(L_ss)
%disp(U_ss)

%% Paso 2
% Realizar una conjetura inicial sobre el capital agregado de la economía:

K_0=(alpha/(1/beta-1+delta))^(1/(1-alpha))*L_ss; % Estado estacionario de %una economía de AR
s_0 = rep*(1-alpha)*(K_0/L_ss)^alpha;
%disp(K_0);
%K_0 = 33.2869; % Equilibrio

diff_crit_k=1;

while diff_crit_k>crit_k
    
    Kiter = [Kiter; K_0];
    
    disp('==========================');
    disp('Iteración');
    disp(count_k);

    %% Paso 3
    % Calcular el salario y la tasa de interés bajo la conjetura:

    w_0=(1-alpha)*(K_0/L_ss)^alpha;
    q_0=alpha*(L_ss/K_0)^(1-alpha);
    %disp(w_0);
    %disp(q_0);
    %% Paso 4
    % Calcular la tasa de impuesto que satisface la restricción presupuestaria del 
    % gobierno bajo la conjetura.
    
    tau_0 = ((1-L_ss)*s_0)/(L_ss*w_0 + K_0*(q_0-delta));
    %disp(tau_0);
    %% Paso 5
    % Calcular las funciones de política para las familias bajo la conjetura:

    kgrid = linspace(0,kmax,N_points_valueiteration);

    % Inicializar la función valor. Dos opciones: (1) valores inciales
    % guardados de la iteración anterios o (2) valores iniciales poco
    % informativos v=0.
    if saved_init == 1
        load initvalfun;
    else
        val_emp_0 = zeros(N_points_valueiteration,1);
        val_des_0 = zeros(N_points_valueiteration,1);
    end

    diff_crit = 1;
    opts=optimset('Diagnostics','off','Display','off');

    utf = @(x) (x.^(1-sigma)-1)/(1-sigma);

    val_emp_1 = zeros(N_points_valueiteration,1);
    val_des_1 = zeros(N_points_valueiteration,1); 
    kdeci_emp = zeros(N_points_valueiteration,1);
    kdeci_des = zeros(N_points_valueiteration,1);

    c_emp = @(k,kf) (1 + (1-tau_0)*(q_0-delta))*k + (1-tau_0)*w_0 - kf;
    c_des = @(k,kf) (1 + (1-tau_0)*(q_0-delta))*k + s_0 - kf;

    count = 1;

    while diff_crit>crit

        Tv_emp = @(k, kf) utf(c_emp(k,kf)) + ... 
            beta*(pp(2,1)*interp1(kgrid,val_des_0,kf,'spline') + pp(2,2)*interp1(kgrid,val_emp_0,kf,'spline')); 
        Tv_des = @(k, kf) utf(c_des(k,kf)) + ...
            beta*(pp(1,1)*interp1(kgrid,val_des_0,kf,'spline') + pp(1,2)*interp1(kgrid,val_emp_0,kf,'spline'));

        for i = 1:N_points_valueiteration
            maxk_emp = (1 + (1-tau_0)*(q_0-delta))*kgrid(i) + (1-tau_0)*w_0;
            maxk_des = (1 + (1-tau_0)*(q_0-delta))*kgrid(i) + s_0;

            kdeci_emp(i) = fmincon(@(x) -Tv_emp(kgrid(i),x), maxk_emp*0.5, [], [], [], [], 0, maxk_emp, [], opts);
            val_emp_1(i) = Tv_emp(kgrid(i), kdeci_emp(i));

            kdeci_des(i) = fmincon(@(x) -Tv_des(kgrid(i),x), maxk_des*0.5, [], [], [], [], 0, maxk_des, [], opts);
            val_des_1(i) = Tv_des(kgrid(i), kdeci_des(i));
        end

        % Criterio para evaluar convergencia
        val_1 = [val_emp_1; val_des_1];
        val_0 = [val_emp_0; val_des_0];
        diff_crit = abs(max(val_1-val_0));

        if flag_print_conv==1
            disp(diff_crit);
        end

        % Actualizamos conjetura
        val_emp_0 = step*val_emp_1 + (1-step)*val_emp_0;
        val_des_0 = step*val_des_1 + (1-step)*val_des_0;

        count = count+1;
        if count==maxiter
            disp('Máximo número de iteraciones alcanzado, no hubo convergencia');
            break;
        end

    end

    save initvalfun val_emp_0 val_des_0;
    disp('Convergencia alcanzada en el problema de las familias!!');
    %% Paso 6
    % Dadas las funciones de política, calcular la distribución de capital de estado 
    % estacionario tanto para agentes empleados como desempleados. 
    % 
    % Para esto es necesario primero invertir la función de política: $k=g_{0}^{-1}(k',\epsilon)$
    % 
    % Invertir la función de política

    fpol_emp = @(kf,k) kf - interp1(kgrid,kdeci_emp,k,'spline');
    fpol_des = @(kf,k) kf - interp1(kgrid,kdeci_des,k,'spline');

    kdeci_emp_inv = zeros(N_points_valueiteration,1);
    kdeci_des_inv = zeros(N_points_valueiteration,1);

    for i = 1:N_points_valueiteration
        kdeci_emp_inv(i) = fzero(@(x) fpol_emp(kgrid(i),x),kgrid(i));
        kdeci_des_inv(i) = fzero(@(x) fpol_des(kgrid(i),x),kgrid(i));    
    end
    %% 
    % Para hallar la distribución invariante iteramos:
    % 
    % $$F_{0,j+1}(k',\epsilon')  =\sum_{\epsilon\in\{0,1\}}\pi(\epsilon'|\epsilon)F_{0,j}(g_{0}^{-1}(k',\epsilon),\epsilon)$$
    % 
    % o lo que es lo mismo:
    % 
    % $$\begin{array}{rcl}F_{0,j+1}(k',1) & =p_{ue}F_{0,j}(g_{0}^{-1}(k',0),0)+p_{ee}F_{0,j}(g_{0}^{-1}(k',1),1)\\F_{0,j+1}(k',0) 
    % & =p_{uu}F_{0,j}(g_{0}^{-1}(k',0),0)+p_{eu}F_{0,j}(g_{0}^{-1}(k',1),1)\end{array}$$

    % Definimos el grid de puntos en k' (usamos un grid más fino que el anterior)
    kprima = linspace(0,kmax,N_points_distinv)';

    % Inicializamos la distribución invariante
    F_emp = ones(N_points_distinv,1);
    F_des = ones(N_points_distinv,1);
    F_emp(kprima<K_0) = 0; 
    F_des(kprima<K_0) = 0; 

    % Generamos las funciones de política invertidas para evaluarlas
    % posteriormente
    f_kdeci_emp_inv = @(x) interp1(kgrid,kdeci_emp_inv,x,'spline');
    f_kdeci_des_inv = @(x) interp1(kgrid,kdeci_des_inv,x,'spline');

    diff_crit =1;
    while diff_crit>crit

        k_emp_i = f_kdeci_emp_inv(kprima);
        k_des_i = f_kdeci_des_inv(kprima);

        f_F_emp = interp1(kprima, F_emp, k_emp_i, 'spline');
        f_F_des = interp1(kprima, F_des, k_des_i, 'spline');  

        % Restricciones sobre la función de distribución acumulada
        f_F_emp(k_emp_i<=0) = 0;
        f_F_des(k_des_i<=0) = 0;
        f_F_emp(k_emp_i>=kmax) = L_ss;
        f_F_des(k_des_i>=kmax) = 1-L_ss;

        % Evolución de la distribución
        F_emp_new = pp(1,2)*f_F_des + pp(2,2)*f_F_emp;
        F_des_new = pp(1,1)*f_F_des + pp(2,1)*f_F_emp;

        % Criterio para evaluar convergencia
        F_new = [F_emp_new; F_des_new];
        F_0   = [F_emp; F_des];

        diff_crit = abs(max(F_new-F_0));

        if flag_print_conv==1
            disp(diff_crit);
        end

        % Actualizamos conjetura
        F_emp_tem = step_dist*F_emp_new + (1-step_dist)*F_emp;
        F_des_tem = step_dist*F_des_new + (1-step_dist)*F_des;

        % Normalization
        F_emp = F_emp_tem*(L_ss/f_F_emp(kmax));
        F_des = F_des_tem*((1-L_ss)/f_F_des(kmax));

    end
    disp('Convergencia alcanzada en la distribución invariante!!');
    %% 
    % Paso 7
    % 
    % Dada la distribución invariante, calcular el stock de capital $K_1$ consistente 
    % con dicha distribución:
    % 
    % $$K_{1} =\int_{0}^{\infty}kf_{0}(k,0)dk+\int_{0}^{\infty}kf_{0}(k,1)dk$$
    % 
    % Usamos integración numérica:

    f_emp = diff(F_emp)./diff(kprima);
    f_des = diff(F_des)./diff(kprima);

    integrando_emp = @(k) k.*interp1(kprima(2:end), f_emp, k, 'spline');
    integrando_des = @(k) k.*interp1(kprima(2:end), f_des, k, 'spline');

    K_1 = integral(integrando_emp,0,kmax) + integral(integrando_des,0,kmax);

    disp(' ');    
    disp('Conjetura vs Actualización');
    disp([K_0, K_1, K_0-K_1]);

    diff_crit_k = abs(K_1-K_0);

    K_0 = step_k*K_1+(1-step_k)*K_0;

    count_k = count_k+1;

end
save results;

%% Ahora podemos presentar los resultados gráficamente:

% Distribución de capital en la economía
plot(kprima(2:end),[f_emp f_des]);
title('Distribución de Capital');
xlabel('k');
ylabel('Densidad');
xlim([-1 200]);
ylim([0 0.025]);
% Convergencia 
plot(Kiter);
title('Convergencia del Capital Agregado');
xlabel('Iteraciones');
ylabel('Capital');

%% 
%