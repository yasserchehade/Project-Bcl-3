clear , close, clc
%% defining varialbles
ax=-pi;
bx=-ax;
ay=ax;
by=bx;
N=100;
M=500;
dx=(bx-ax)/(N+1);
dy=(by-ay)/(N+1);
T=1;
dt=T/(M+1);
r=dt/(dx)^2;
%% laying out the axis
x=zeros(N+2,1);
x(1)=-pi;
for ii=1:N+1
    x(ii+1)=-pi+dx*ii;
end

y=zeros(N+2,1);
y(1)=-pi;
for ii=1:N+1
    y(ii+1)=-pi+dy*ii;
end

t=zeros(M+2,1);
for ii=2:(M+2)
    t(ii)=dt*(ii-1);
end

u=zeros(N+2,N+2,M+2);
%% inputing the boundary conditions 

% u(ax,y,t)=(by-y)^2*cos(pi*y/by)
for ii=1:N+2
    u(1,ii,:)=(by-y(ii))^2*cos(pi*y(ii)/by);
end
u_ax_y_t=u(1,2:N+2,1);
% u(bx,y,t)=y*(by-y)^2
for ii=1:N+2
    u(N+2,ii,:)=y(ii)*(by-y(ii))^2;
end
u_bx_y_t=u(N+2,2:N+2,1);
% f(y)=(by-y)^2*cos(pi*y/by)
f_ay=(by-ay)^2*cos(pi*ay/by);
%g(y)=y*(by-y)^2
g_ay=ay*(by-ay)^2;

%u(x,ay,t)=f_ay+(x-ax)/(bx-ax)*(g_ay-f_ay)

for ii=1:N+2
    u(ii,1,:)=f_ay+(x(ii)-ax)/(bx-ax)*(g_ay-f_ay);
end
u_x_ay_t=u(2:N+1,1,1);

%% explicit method
for nn=2:M+2
    for jj=2:N+1
        for kk=2:N+1
            u(jj,kk,nn)=r*u(jj-1,kk,nn-1)+(1-4*r)*u(jj,kk,nn-1)+r*u(jj+1,kk,nn-1)+r*u(jj,kk-1,nn-1)+r*u(jj,kk+1,nn-1);
        end
        u(jj,N+2,nn)=r*u(jj-1,N+2,nn-1)+(1-4*r)*u(jj,N+2,nn-1)+r*u(jj+1,N+2,nn-1)+2*r*u(jj,N+1,nn-1);
    end
end

% %% defining algorithm parameters
% B=[r*ones(N^2+N,1),2*(1-r)*ones(N^2+N,1),r*ones(N^2+N,1)];
% 
% d=[-N 0 N];
% B=spdiags(B,d,N^2+N,N^2+N);
% 
% B2=[r*ones(N+1,1),2*(1-r)*ones(N+1,1),r*ones(N+1,1)];
% B2=repmat(B2,N,1);
% d=[-N 0 N];
% B2=spdiags(B2,d,N^2+N,N^2+N);
% d=r*diag(ones(N,1),-N);
% B2(N^2-N+1:N^2+N,N^2-N+1:N^2+N)=B2(N^2-N+1:N^2+N,N^2-N+1:N^2+N)+d;
% 
% a=2*(r+1)*ones(N+1,1);
% b=zeros(N+1,1);
% b(1:N)=-r*ones(N,1);
% c=zeros(N+1,1);
% c(2:N)=-r*ones(N-1,1);
% c(N+1)=-2*r;
% 
% a2=2*(r+1)*ones(N,1);
% b2=zeros(N,1);
% b2(1:N-1)=-r*ones(N-1,1);
% c2=zeros(N,1);
% c2(2:N)=-r*ones(N-1,1);
% 
% u_n_half=zeros(N^2+N,1);
% u_n_1=zeros(N^2+N,1);
% 
% 
% for ii=1:M+1
%    un=u(2:N+1,2:N+2,ii);
%    un=reshape(un',N^2+N,1);
%    f=B*un;
%    f=reshape(f,N+1,N);
%    f=f';
%    f(:,1)=f(:,1)+u_x_ay_t;
%    f(1,:)=f(1,:)+u_ax_y_t;
%    f(N,:)=f(N,:)+u_bx_y_t;
%    f=reshape(f',N^2+N,1);
%    for jj=1:N
%        u_n_half(1+(jj-1)*(N+1):jj*(N+1))=tridiag(a,c,b,f(1+(jj-1)*(N+1):jj*(N+1)));
%    end
%    u_n_half=reshape(u_n_half,N+1,N);
%    u_n_half=reshape(u_n_half',N^2+N,1);
%    f=B2*u_n_half;
%    f=reshape(f,N,N+1);   
%    f(:,1)=f(:,1)+u_x_ay_t;
%    f(1,:)=f(1,:)+u_ax_y_t;
%    f(N,:)=f(N,:)+u_bx_y_t;
%    f=reshape(f,N^2+N,1);
%    for jj=1:N+1
%        u_n_1(1+(jj-1)*N:jj*N)=tridiag(a,c,b,f(1+(jj-1)*N:jj*N));
%    end
%    u(2:N+1,2:N+2,ii+1)=reshape(u_n_1,N,N+1);   
% end

% ploting result
% [xx,yy]=meshgrid(x',y);
% for ii=1:M+2
%     figure
%     surf(yy,xx,u(:,:,ii));
%     xlabel('x axis')
%     ylabel('y axis')
% %     title(['u(x,y) for t=' num2str(time) ' sec' ])
% end
fr=figure;
[xx,yy]=meshgrid(x',y);
filename='testanimated.gif';
for ii=1:M+2
    time=t(ii);
    s=num2str(time);
    cla;
    surf(yy,xx,u(:,:,ii));
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




