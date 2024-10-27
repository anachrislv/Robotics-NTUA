clear all; 
close all;
l0 = 18;
l1 = 20.0;  %% in cm 
l2 = 35.0;
l3 = 30.0;
 
dt = 0.01;    
    
xd0 = 20.0;	%starting point
xd1 =  20.0; %ending point
yd0 = 5.00; %starting point
yd1 = 20.00;  %ending point

xd(1) = xd0; 
yd(1) = 5.00; 

Tf = 15.0;
u_max = 1.5; %maximum velocity only for the y-coordinate (no motion in the x axis)
a = 0.3; %acceleration
kmax=Tf/dt + 1; 
Tf=2*Tf; %back and forth motion
t=0:dt:Tf;
t2=t;
t2(end)=[]
time=dt;
u_y=0; %initially
u_x=0; %for every t (no motion in the x axis)
u_y1(1)=0;
katw_orio=kmax./3;
panw_orio=2*kmax./3;
yd1(1)=5

%equations for the forth motion
for k=2:kmax;
   if k < katw_orio +1
       %accelerating motion
       yd1(k)=yd1(k-1) + u_y*dt+ 0.5 * a * dt*dt;
       u_y=u_y+a*dt;
       u_y1(k)=u_y;
   elseif (k> (katw_orio ) ) && (k < panw_orio+1)
       %constant velocity
       yd1(k) = yd1(k-1) + u_max*dt;
       u_y=u_max;
       u_y1(k)=u_y;
       
   elseif k> panw_orio
       %decelerating motion
       yd1(k)=yd1(k-1)+u_y*dt - 0.5 * a * dt*dt;
       u_y=u_y-a*dt;
       u_y1(k)=u_y;
       
   end

   
   
end  



yd2(kmax)=20;
yd2(1)=5
k=kmax-1;
u_y=0;
u_y2(kmax)=0;
%equations for the back motion

while k>1
    if k> (panw_orio)
        %accelerating
        yd2(k)=yd2(k+1) - u_y*dt - 0.5 * a * dt *dt;
        u_y=u_y+a*dt;
        u_y2(k)=u_y;
    elseif k> (katw_orio) && k < (panw_orio+1)
        %constant velocity
        yd2(k) = yd2(k+1) - u_max*dt;
        u_y=u_max;
        u_y2(k)=u_y;
    elseif k < (katw_orio+1)
        %decelerating
        yd2(k)=yd2(k+1) - u_y*dt + 0.5 * a * dt *dt;
        u_y=u_y-a*dt;
        u_y2(k)=u_y;

    end
    k=k-1;
end

yd2=flip(yd2);
i=kmax+1;

for k=1:2*kmax-1
    if k < kmax+1
        yd(k)=yd1(k);
        u_y(k)=u_y1(k);
    else 
        yd(k)=yd2(k-kmax+1); %merging the 2 motion arrays
        u_y(k)=u_y2(k-kmax+1);%merging the 2 velocity arrays
    end 


end


for i=1:2*kmax-1;
    xd(i)=20;
    zd(i)=20;
    u_x(i)=0;
    u_z(i)=0;
end



%plotting motion in the y axis
fig1=figure;
plot(t,yd); 
ylabel('y axis motion (cm)'); 
xlabel('time t (sec)'); 

%plotting motion in the x axis
fig2=figure;
plot(t,xd); 
ylabel('x axis motion (cm)'); 
xlabel('time t (sec)'); 


%plotting motion in the z axis
fig3=figure;
plot(t,zd); 
ylabel('z axis motion (cm)'); 
xlabel('time t (sec)'); 

%plotting velocity in the y axis
fig4=figure;
plot(t,u_y); 
ylabel('y axis velocity (cm)'); 
xlabel('time t (sec)'); 


%plotting velocity in the x axis
fig5=figure;
plot(t,u_x); 
ylabel('x axis velocity (cm)'); 
xlabel('time t (sec)'); 



%plotting velocity in the z axis
fig6=figure;
plot(t,u_z); 
ylabel('z axis velocity (cm)'); 
xlabel('time t (sec)'); 

 %q3
N = xd(:).^2 + yd(:).^2 + (zd(:)-l0).^2 -(l3)^2 - (l2)^2 - (l1)^2;
D = 2*l3*l2;
C3 = N./D;
S3 = sqrt(1-(C3).^2);
q3 = atan2(S3,C3);


%q2
n1 = cos(q3).*l3 + l2;
n2 = sin(q3).*l3;
den = sqrt((n1).^2 + (n2).^2);
b = atan2(n2,n1);
v = yd(:)./den;
q2 = asin(v) - b;


n3 = cos(q2+q3).*l3 + cos(q2).*l2;
d3 = l1;
b = atan2(d3,n3);
r = sqrt((n3).^2 + (d3).^2);
w = xd(:)./r;
q1 = acos(w) - b;

xd1 = zeros(2*kmax,1);
yd1 = zeros(2*kmax,1);
zd1 = l0.*ones(2*kmax,1);

xd2 = (-1)*l1.*sin(q1);
yd2 = l1.*cos(q1);
zd2 = l0.*ones(2*kmax,1);

xd3 = l2.*cos(q2).*cos(q1) - l1.*sin(q1);
yd3 = l2.*sin(q2);
zd3 = l2.*sin(q1).*cos(q2) + l1.*cos(q1) + l0;

fig7 = figure; 
q_d1=diff(q1(:,1))
u_q2=diff(q2)./diff(t2);
plot(t2, u_q2);  
ylabel('q2 velocity '); 
xlabel('time t (sec)');  

fig8 = figure; 
q_d1=diff(q1(:,1));
u_q1=diff(q1)./diff(t2);
plot(t2, u_q1);  
ylabel('q1 velocity '); 
xlabel('time t (sec)');  

fig9 = figure; 
q_d1=diff(q1(:,1))
u_q3=diff(q3)./diff(t2);
plot(t2, u_q3);  
ylabel('q3 velocity '); 
xlabel('time t (sec)');  

fig10 = figure;
plot(t, q1);
ylabel('q1 '); 
xlabel('time t (sec)'); 


fig11 = figure;
plot(t, q2);
ylabel('q2 '); 
xlabel('time t (sec)'); 


fig12 = figure;
plot(t, q3);
ylabel('q3 '); 
xlabel('time t (sec)');


fig14 = figure; 
axis([-5 70 -20 40]) 
axis on 
hold on 
xlabel('x '); 
ylabel('z '); 
loops = 5;
plot(xd,zd,'rs'); 
dtk=100;
plot([0],[0],'o'); 
for n=1:loops;  
  for tk=1:dtk:3001;    
    pause(1);	%% pause motion to view successive robot configurations    
     plot([0,xd1(tk)],[0,zd1(tk)]);	
     plot([xd1(tk)],[zd1(tk)],'o');
     plot([xd1(tk),xd2(tk)],[zd1(tk),zd2(tk)]);	
     plot([xd2(tk)],[zd2(tk)],'o'); 
     plot([xd2(tk),xd3(tk)],[zd2(tk),zd3(tk)]);	
     plot([xd3(tk)],[zd3(tk)],'o');
     plot([xd3(tk),xd(tk)],[zd3(tk),zd(tk)]);
    plot([xd(tk)],[zd(tk)],'g+');  
  end
end

fig13 = figure; 
axis([-5 70 -20 40]) %
axis on 
hold on 
xlabel('x (cm)'); 
ylabel('y (cm)'); 
loops = 5;
plot(xd,yd,'rs'); 
dtk=100; % plot robot position every dtk samples, to animate its motion 
plot([0],[0],'o'); 
for n=1:loops;  
  for tk=1:dtk:3001;    %%% 	
    pause(1);	%% pause motion to view successive robot configurations    
     plot([0,xd1(tk)],[0,yd1(tk)]);	
     plot([xd1(tk)],[yd1(tk)],'o');
     plot([xd1(tk),xd2(tk)],[yd1(tk),yd2(tk)]);	
     plot([xd2(tk)],[yd2(tk)],'o'); 
     plot([xd2(tk),xd3(tk)],[yd2(tk),yd3(tk)]);	
     plot([xd3(tk)],[yd3(tk)],'o');
     plot([xd3(tk),xd(tk)],[yd3(tk),yd(tk)]);
    plot([xd(tk)],[yd(tk)],'g+');  
  end
end
