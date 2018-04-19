clear , close, clc
%% defining varialbles
ax=-pi;
bx=-ax;
ay=ax;
by=bx;
N=10;
M=5;
dx=(bx-ax)/(N+1);
dy=(by-ay)/(N+1);
T=1;
dt=T/(M+1);
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

t=zeros(2*M+3,1);
for ii=1:2*(M+1)
    t(ii+1)=dt/2*ii;
end

u=zeros(N+2,N+2,2*M+3);
%% inputing the boundary conditions 

% u(ax,y,t)=(by-y)^2*cos(pi*y/by)
for ii=1:N+2
    u(1,ii,:)=(by-y(ii))^2*cos(pi*y(ii)/by);
end

% u(bx,y,t)=y*(by-y)^2
for ii=1:N+2
    u(N+2,ii,:)=y(ii)*(by-y(ii))^2;
end

% f(y)=(by-y)^2*cos(pi*y/by)
f_ay=(by-ay)^2*cos(pi*ay/by);
%g(y)=y*(by-y)^2
g_ay=ay*(by-ay)^2;

%u(x,ay,t)=f_ay+(x-ax)/(bx-ax)*(g_ay-f_ay)

for ii=1:N+2
    u(ii,1,:)=f_ay+(x(ii)-ax)/(bx-ax)*(g_ay-f_ay);
end

%% defining algorithm parameters
B=[r*ones(N^2+N,1),2*(1-r)*ones(N^2+N,1),r*ones(N^2+N,1)];

d=[-N 0 N];
B=spdiags(B,d,N^2+N,N^2+N);
% % A_n=2*(r+1)*eye(N)-r*diag(ones(N-1,1),1)-r*diag(ones(N-1,1),-1);
% % tmp=repmat({A_n},N,1);
% % A=blkdiag(tmp{:});
% % A=sparse(A);
% % a=diag(A);
% % b=diag(A,1);
% % c=diag(A,-1);
a=2*(r+1)*ones((N^2)+N,1);
b=zeros(N+1,1);
b(1:N)=-r*ones(N,1);
b=repmat(b,N,1);
c=zeros(N+1,1);
c(2:N)=-r*ones(N-1,1);
c(N+1)=-2*r;
c=repmat(c,N,1);

for ii=1:M+1
    un=u(2:N+1,2:N+2,(2*ii-1));
    un=reshape(un',N^2+N,1);
    f=B*un;
    f=reshape(f,N,N+1);
    f(:,1)=f(:,1)+r*u(2:N+1,1,2*ii);
    f(1,:)=f(1,:)+r*u(1,2:N+2,2*ii-1);
    f(N,:)=f(N,:)+r*u(N+2,2:N+2,2*ii-1);
    f=reshape(f',N^2+N,1);
    u_n_half=tridiag(a,b,c,f);
    u_n_half=reshape(u_n_half,N,N+1);
    u(2:N+1,2:N+2,2*ii)=u_n_half;
    
    u_n_half=reshape(u_n_half',N^2+N,1);
    f=B*u_n_half;
    f=reshape(f,N,N+1);
    f(:,1)=f(:,1)+r*u(2:N+1,1,2*ii);
    f(1,:)=f(1,:)+r*u(1,2:N+2,2*ii+1);
    f(N,:)=f(N,:)+r*u(N+2,2:N+2,2*ii+1);
    f=reshape(f',N^2+N,1);
    u_n_1=tridiag(a,b,c,f);
    u_n_1=reshape(u_n_1,N,N+1);
    u(2:N+1,2:N+2,2*ii+1)=u_n_1;
    
end

%% ploting result
fr=figure;
[xx,yy]=meshgrid(x,y);
filename='testanimated.gif';
for ii=1:2*M+3
    time=t(ii);
    s=num2str(time);
    cla;
    surf(xx,yy,u(:,:,ii));
    xlabel('x axis')
    ylabel('y axis')
    title(['u(x,y) for t=' num2str(time) ' sec' ])
    drawnow
    frame=getframe(fr);
    im=frame2im(frame);
    [imind,cm]=rgb2ind(im,256);
    
    if ii==1
        imwrite(imind,cm,filename,'gif','loopcount',inf);
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
end


