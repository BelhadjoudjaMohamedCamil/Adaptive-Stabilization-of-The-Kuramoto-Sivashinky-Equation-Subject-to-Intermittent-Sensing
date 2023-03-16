%Adaptive Stabilization Of The Kuramoto-Sivashinsky Equation Subject
%To Intermittent Sensing

%American Control Conference (ACC) 2023, San Diego.

%M.C. Belhadjoudja, M. Maghenem, E. Witrant, C. Prieur

%  We simulate the system of two KS equations :
%     wt + wwx + lambda1 wxx + wxxxx = 0, x in [0,Y]
%     vt + vvx + lambda1 vxx + vxxxx = 0, x in [Y,1]
%     w(0) = U1(t), w(Y) = v(Y) = 0, v(1) = U2(t),
%     wx(0) = wx(Y) = vx(Y) = vx(1) = 0;
%  U1 and U2 are the control inputs.

clc; clear all; close all;

%Declaration of constants
Y = 0.5; %Length of the first space interval.
N = 9; %The number of collocation points for each PDE is (N+1).
dx = Y/N; %Space step.
x_vecw = 0:dx:Y; %Space vector for first PDE.
x_vecv = Y:dx:2*Y; %Space vector for second PDE.
tf = 0.008; %Final time of simulation.
dt = 0.0000001; %Time step.
t_vec = 0:dt:tf; %Time vector.
Lambda1 = (4*pi^2)/(Y^2)+50; %Destabilizing parameter (anti-diffusion).

%Declaration of the RBF matrices (A,D1,...,D4). We use multiquadrics
%RBF (MQ), namely : Phi(r) = sqrt(r^2+c^2), where c is a shape parameter.
c = 0.4;
A = zeros(N+1,N+1); D1 = A; D2 = A; D4 = A; D3 = A;
for i = 1:(N+1)
    ki = i-1;
    for j = 1:(N+1)
        kj = j-1;
        x = abs(ki-kj)*dx;
        A(i,j) = (c^2 + x^2)^(1/2);
        D1(i,j) = x/(c^2 + x^2)^(1/2);
        D2(i,j) = 1/(c^2 + x^2)^(1/2) - x^2/(c^2 + x^2)^(3/2);
        D3(i,j) = (3*x^3)/(c^2 + x^2)^(5/2) - (3*x)/(c^2 + x^2)^(3/2);
        D4(i,j) = (18*x^2)/(c^2 + x^2)^(5/2) - 3/(c^2 + x^2)^(3/2) - (15*x^4)/(c^2 + x^2)^(7/2);
    end
end

%Initialization of the different signals.
W = zeros(N+1,length(t_vec)); Wx = W; V = W; Vx = Wx; %The PDEs states
U1 = zeros(length(t_vec),1); U2 = U1; %The controls
V1 = zeros(length(t_vec),1); V2 = V1; %The Lyapunov function candidates
hattheta1 = 0*ones(length(t_vec),1); hattheta2 = 0*ones(length(t_vec),1); %adaptation parameters
wxxx = zeros(length(t_vec),1); vxxx = wxxx; %Third derivatives at 0 and 1

%Initial condition
W(:,1) = -3*(cos(4*pi*x_vecw)-1);
V(:,1) = -3*(cos(4*pi*x_vecv)-1);

%Initial value of the Lyapunov function candidates (Riemannian sums)
 for idx = 1:length(x_vecw)
        V1(1) = (1/2)*dx*W(idx,1)^2+V1(1);
 end
 for idx = 1:length(x_vecv)
        V2(1) = (1/2)*dx*V(idx,1)^2+V2(1);
 end

%The time intervals I1 = [0,t1)U[t2,t3)U[t4,t5)U[t6,t7)U[t8,t9),
%I2 = [t1,t2)U[t3,t4)U[t5,t6)U[t7,t8)U[t9,tf).
t1 = 0.001;
t2 = 0.002;
t3 = 0.0028;
t4 = 0.0039;
t5 = 0.005;
t6 = 0.0055;
t7 = 0.0065;
t8 = 0.007;
t9 = 0.0076;
t10 = tf;

%Control parameters
beta = 3;
sigma = 100; %Desired minimum decay rate
increment = 0.01; %Increment on hattheta1 and hattheta2

%Running the control algorithm
for idt = 1:(length(t_vec)-1)
   
    %Computing first derivatives
    if idt == 1
        Wx(:,idt) = 12*pi*sin(4*pi*x_vecw); %Exact derivative of initial condition
        Vx(:,idt) = 12*pi*sin(4*pi*x_vecv);
    else
        Wx(3:(end-3),idt) = (1/dx)*(W(3:(end-3),idt)-W(2:(end-4),idt)); %Euler backward scheme
        Vx(3:(end-3),idt) = (1/dx)*(V(3:(end-3),idt)-V(2:(end-4),idt));
    end  

    %Computing the third derivatives at 0 and 1 to design the controls.
    wxxx(idt) = (1/dx^3)*(W(4,idt)-3*W(3,idt)+3*W(2,idt)-W(1,idt)); %Euler forward scheme
    vxxx(idt) = (1/dx^3)*(V(end,idt)-3*V((end-1),idt)+3*V((end-2),idt)-V((end-3),idt)); %Euler backward scheme
   
    %Updating hattheta1 and hattheta1
    if ((t_vec(idt)>=t2 && t_vec(idt)<t3) && (V1(t2/dt) <= exp(-sigma*t2)*V1(1)))||((t_vec(idt)>=t4 && t_vec(idt)<t5) && (V1(t4/dt) <= exp(-sigma*(t4-t2))*V1(t2/dt)))||((t_vec(idt)>=t6 && t_vec(idt)<t7) && (V1(t6/dt) <= exp(-sigma*(t6-t4))*V1(t4/dt)))||((t_vec(idt)>=t8 && t_vec(idt)<t9) && (V1(t8/dt) <= exp(-sigma*(t8-t6))*V1(t6/dt)))
        hattheta1(idt) = hattheta1(idt-1);
        hattheta2(idt) = hattheta2(idt-1);
    elseif ((t_vec(idt)>=t3 && t_vec(idt)<t4) && (V2(t3/dt) <= exp(-sigma*(t3-t1))*V2(t1/dt)))||((t_vec(idt)>=t5 && t_vec(idt)<t6) && (V2(t5/dt) <= exp(-sigma*(t5-t3))*V2(t3/dt)))||((t_vec(idt)>=t7 && t_vec(idt)<t8) && (V2(t7/dt) <= exp(-sigma*(t7-t5))*V2(t5/dt)))||((t_vec(idt)>=t9 && t_vec(idt)<t10) && (V2(t9/dt) <= exp(-sigma*(t9-t7))*V2(t7/dt)))
        hattheta2(idt) = hattheta2(idt-1);
        hattheta1(idt) = hattheta1(idt-1);
    elseif t_vec(idt)>= 0 && t_vec(idt)<t2
        hattheta1(idt) = hattheta1(1);
        hattheta2(idt) = hattheta2(1);
    elseif ((t_vec(idt)>=t2 && t_vec(idt)<t3))||((t_vec(idt)>=t4 && t_vec(idt)<t5))||((t_vec(idt)>=t6 && t_vec(idt)<t7))||((t_vec(idt)>=t8 && t_vec(idt)<t9))
        hattheta1(idt) = hattheta1(idt-1)+increment;
        hattheta2(idt) = hattheta2(idt-1);
    elseif ((t_vec(idt)>=t3 && t_vec(idt)<t4))||((t_vec(idt)>=t5 && t_vec(idt)<t6))||((t_vec(idt)>=t7 && t_vec(idt)<t8))||((t_vec(idt)>=t9 && t_vec(idt)<t10))
        hattheta2(idt) = hattheta2(idt-1)+increment;
        hattheta1(idt) = hattheta1(idt-1);
    end

    %Parameters used to update the control inputs.
    L1 = (1/3)*(1+3*hattheta1(idt))*V1(idt)^(2/3);
    k1 = -(1+3*hattheta1(idt))*beta*V1(idt)^(1/3);
    L2 = (1/3)*(1+3*hattheta2(idt))*V2(idt)^(2/3);
    k2 = (1+3*hattheta2(idt))*beta*V2(idt)^(1/3);
   
    %Updating the control inputs
    if (t_vec(idt)>=0 && t_vec(idt)<t1)||(t_vec(idt)>=t2 && t_vec(idt)<t3)||(t_vec(idt)>t4 && t_vec(idt)<t5)||(t_vec(idt)>t6 && t_vec(idt)<t7)||(t_vec(idt)>t8 && t_vec(idt)<t9)
        U2(idt+1) = 0;
        if abs(wxxx(idt)) >= L1
            U1(idt+1) = -sign(wxxx(idt))*V1(idt)^(1/3);
        else
            U1(idt+1) = k1;
        end
    elseif t_vec(idt) >= 0
        U1(idt+1) = 0;
        if abs(vxxx(idt)) >= L2
            U2(idt+1) = sign(vxxx(idt))*V2(idt)^(1/3);
        else
            U2(idt+1) = k2;
        end    
    end
   
    %Uncomment the following two lines to obtain the open-loop simulations.
    %U1(idt+1)=0;
    %U2(idt+1)=0;
   
    %Updating the state W.
    Mw = A + (dt/2)*(diag(W(:,idt))*[zeros(2,N+1);D1(3:(end-2),:);zeros(2,N+1)]+diag(Wx(:,idt))*A+Lambda1*[zeros(2,N+1);D2(3:(end-2),:);zeros(2,N+1)]+[zeros(2,N+1);D4(3:(end-2),:);zeros(2,N+1)]);
    Nw = [zeros(2,N+1);A(3:(end-2),:);zeros(2,N+1)]-(dt/2)*(Lambda1*[zeros(2,N+1);D2(3:(end-2),:);zeros(2,N+1)]+[zeros(2,N+1);D4(3:(end-2),:);zeros(2,N+1)]);
    Gw = [U1(idt+1);U1(idt+1);zeros(N-1,1)];
    W(:,idt+1) = A*inv(Mw)*Nw*inv(A)*W(:,idt)+A*inv(Mw)*Gw;
   
    %Updating the state V.
    Mv = A + (dt/2)*(diag(V(:,idt))*[zeros(2,N+1);D1(3:(end-2),:);zeros(2,N+1)]+diag(Vx(:,idt))*A+Lambda1*[zeros(2,N+1);D2(3:(end-2),:);zeros(2,N+1)]+[zeros(2,N+1);D4(3:(end-2),:);zeros(2,N+1)]);
    Nv = [zeros(2,N+1);A(3:(end-2),:);zeros(2,N+1)]-(dt/2)*(Lambda1*[zeros(2,N+1);D2(3:(end-2),:);zeros(2,N+1)]+[zeros(2,N+1);D4(3:(end-2),:);zeros(2,N+1)]);
    Gv = [zeros(N-1,1);U2(idt+1);U2(idt+1)];
    V(:,idt+1) = A*inv(Mv)*Nv*inv(A)*V(:,idt)+A*inv(Mv)*Gv;
   
    %Updating the Lyapunov function candidates.
    for idx = 1:length(x_vecw)
        V1(idt+1) = (1/2)*dx*W(idx,idt+1)^2+V1(idt+1);
    end
    for idx = 1:length(x_vecv)
        V2(idt+1) = (1/2)*dx*V(idx,idt+1)^2+V2(idt+1);
    end

end

%Final value of hattheta1 and hattheta2.
hattheta1(end) = hattheta1(end-1);
hattheta2(end) = hattheta2(end-1);

%The total Lyapunov function candidate.
Lya = V1 + V2;

%Plotting the results.
V = V(2:end,:);
x_vecv = x_vecv(2:end);
figure(1) %Plot the PDEs
[tt,xx] = meshgrid(t_vec,[x_vecw,x_vecv]);
mesh(xx,tt,[W;V])
%colormap hot
xlabel({'Space'},'Interpreter','latex')
ylabel({'Time'},'Interpreter','latex')
zlabel({'$u(t,x)$'},'Interpreter','latex')
figure(2) %Plot the Lyapunov function candidates
plot(t_vec,V1,'linewidth',3)
hold on
plot(t_vec,V2,'linewidth',3)
plot(t_vec,Lya,'linewidth',3)
y = ylim;
plot([t1 t1],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t2 t2],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t3 t3],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t4 t4],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t5 t5],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t6 t6],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t7 t7],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t8 t8],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t9 t9],[y(1) y(2)],'--','color',[0,0,0]+0.5)
legend({'$V_{1}$','$V_{2}$','$W$'},'Interpreter','latex')
xlabel({'Time'},'Interpreter','latex')
ylabel({'Lyapunov function candidates'},'Interpreter','latex')
grid on
figure(3) %Plot the control inputs
plot(t_vec,U1,'linewidth',3)
hold on
plot(t_vec,U2,'linewidth',3)
y = ylim;
plot([t1 t1],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t2 t2],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t3 t3],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t4 t4],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t5 t5],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t6 t6],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t7 t7],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t8 t8],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t9 t9],[y(1) y(2)],'--','color',[0,0,0]+0.5)
grid on
legend({'$u_{1}$','$u_{2}$'},'Interpreter','latex')
xlabel({'Time'},'Interpreter','latex')
ylabel({'Control inputs'},'Interpreter','latex')
figure(5) %Plot hattheta1 and hattheta2
plot(t_vec,hattheta1,'linewidth',3)
hold on
plot(t_vec,hattheta2,'linewidth',3)
y = ylim;
plot([t1 t1],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t2 t2],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t3 t3],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t4 t4],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t5 t5],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t6 t6],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t7 t7],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t8 t8],[y(1) y(2)],'--','color',[0,0,0]+0.5)
plot([t9 t9],[y(1) y(2)],'--','color',[0,0,0]+0.5)
grid on
legend({'$\hat{\theta}_{1}$','$\hat{\theta}_{2}$'},'Interpreter','latex')
xlabel({'Time'},'Interpreter','latex')
ylabel({'Adaptation parameters'},'Interpreter','latex')

