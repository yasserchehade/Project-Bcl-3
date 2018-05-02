clear , close, clc
%% defining varialbles
ax=-pi;
bx=-ax;
ay=ax;
by=bx;
N=100;
M=5000;
dx=(bx-ax)/(N+1);
dy=(by-ay)/(N+1);
T=1;
dt=T/(M+1);
r=dt/(dx)^2;
Q=50;
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
    u(ii,1,:)=f_ay+((x(ii)-ax)/(bx-ax)*(g_ay-f_ay));
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
% U=u;
%% cranck nichilson
% RHS=zeros(N+2,N+2);
% B=sparse(N^2+N,N^2+N);
% B_diag=(2+4*r)*eye(N)+diag(-r*ones(N-1,1),1)+diag(-r*ones(N-1,1),-1);
% B_sup=-r*eye(N);
% B_sup1=-2*r*eye(N);
% 
% for ii=1:N+1
%     B((ii-1)*(N)+1:(ii-1)*(N)+(N),(ii-1)*(N)+1:(ii-1)*(N)+(N))=B_diag;
% end
% 
% for ii=2:N+1
%     B((ii-2)*(N)+1:(ii-2)*(N)+(N),(ii-1)*(N)+1:(ii-1)*(N)+(N))=B_sup;
%     B((ii-1)*(N)+1:(ii-1)*(N)+(N),(ii-2)*(N)+1:(ii-2)*(N)+(N))=B_sup;
% end
% 
% B(N^2+1:N^2+N,N^2-N+1:N^2)=B_sup1;
% 
% 
% for nn=2:M+1
%     for jj=2:N+1
%         for kk=2:N+1
%             RHS(jj,kk)=r*u(jj-1,kk,nn)+(2-4*r)*u(jj,kk,nn)+r*u(jj+1,kk,nn)+r*u(jj,kk-1,nn)+r*u(jj,kk+1,nn);
%         end
%         RHS(jj,N+2)=r*u(jj-1,N+2,nn)+(2-4*r)*u(jj,N+2,nn)+r*u(jj+1,N+2,nn)+2*r*u(jj,N,nn);
%     end
%     RHS=RHS(2:N+1,2:N+2);
%     RHS(1,:)=RHS(1,:)+r*u(1,2:N+2,1);
%     RHS(N,:)=RHS(N,:)+r*u(N+2,2:N+2,1);
%     RHS(:,1)=RHS(:,1)+r*u(2:N+1,1,1);
%     rhs=reshape(RHS,N^2+N,1);
%     
%     un=reshape(u(2:N+1,2:N+2,nn),N^2+N,1);
%     un1=B\rhs;
%     u(2:N+1,2:N+2,nn+1)=reshape(un1,N,N+1);
%     
% end

%%
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




