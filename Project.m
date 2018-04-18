clear , close, clc
%% defining varialbles
ax=-pi;
bx=-ax;
ay=ax;
by=bx;
N=100;
dx=(bx-ax)/(N+1);
dy=(by-ay)/(N+1);
T=1;
dt=T/(N+1);
r=dt/(dx)^2;
%% laying out the axis
x=zeros(N+2,1);
for ii=1:N+1
    x(ii+1)=dx*ii;
end

y=zeros(N+2,1);
for ii=1:N+1
    y(ii+1)=dy*ii;
end

t=zeros(N+2,1);
for ii=1:N+1
    t(ii+1)=dt*ii;
end

u=zeros(N+2,N+2,N+2);
%% inputing the boundary conditions 

% u(ax,y,t)=(by-y)^2*cos(pi*y/by)
for ii=1:N+2
    u(1,:,:)=(by-y(ii))^2*cos(pi*y(ii)/by);
end

% u(bx,y,t)=y*(by-y)^2
for ii=1:N+2
    u(N+2,:,:)=y(ii)*(by-y(ii))^2;
end

% f(y)=(by-y)^2*cos(pi*y/by)
f_ay=(by-y)^2*cos(pi*y/by);
%g(y)=y*(by-y)^2
g_ay=ay*(by-ay)^2;

%u(x,ay,t)=f_ay+(x-ax)/(bx-ax)*(g_ay-f_ay)

for ii=1:N+2
    u(:,1,:)=f_ay+(x(ii)-ax)/(bx-ax)*(g_ay-f_ay);
end

%% 