N=10000; %number of particles
dt=1e-3;    %time step
v_avg=0;    %average of the velocities
sig=400;    %standard deviation for the velocities
tmax=500;   %number of steps
kb=1.38064852e-23;  %boltzmann constant

vels=zeros(N,2);    %zeros array to store xy components of the velocity
pos=zeros(N,2);     %zeros array to store xy coordinates of starting position
intervals=pos(1,1)+6/(sqrt(N)-1);   %calculating the spacing between the particles
n=0;    %simple counter 
for i=1:sqrt(N) 
    for j=1:sqrt(N)
        n=n+1;  %step the counter by one
        pos(n,1)=(i-1)*intervals;   %set nth x coordinate to i
        pos(n,2)=(j-1)*intervals;   %seth nth y coordinate to j
        vels(n,1)=sig.*randn;       %set nth vx, assigning a random value 
        vels(n,2)=sig.*randn;       %set nth vy, assigning a random value 
    end
    
end

%below two arrays are created to store position of particle at each step
%(will make the creation of the animation easier)
pos_x_t=zeros(N,tmax);
pos_x_t(:,1)=pos(:,1);  %first column is starting x coordinate of each particle
pos_y_t=zeros(N,tmax);  
pos_y_t(:,1)=pos(:,2);  %first column is starting y coordinate of each particle

vels_x_t=zeros(N,tmax);

xroom_imp=[7 -7];
yroom_imp=[7 -7];
for t=1:tmax   %step for loop a tmax number of steps
    vels=zeros(N,2);    %reset value of velocities to 0
    for r=1:N           %this for loop cycles through the N particles composing the system
        vels(r,1)=sig.*randn;   %assign random vx to rth particle
        vels_x_t(r,t+1)=vels(r,1);
        vels(r,2)=sig.*randn;   %assign random vy to rth particle
        pos_x_t(r,t+1)=pos_x_t(r,t)+vels(r,1)*dt;   %use simple equation of motion to calculate x position at following step
        pos_y_t(r,t+1)=pos_y_t(r,t)+vels(r,2)*dt;   %use simple equation of motion to calculate y position at following step
        if pos_x_t(r,t+1)>=xroom_imp(1)
            pos_x_t(r,t+1)=-pos_x_t(r,t);
            vels(r,1)=sig.*randn;
        elseif pos_x_t(r,t+1)<=xroom_imp(2)
            pos_x_t(r,t+1)=-pos_x_t(r,t);
            vels(r,1)=sig.*randn;
        end
        if pos_y_t(r,t+1)>=yroom_imp(1)
            pos_y_t(r,t+1)=-pos_y_t(r,t);
            vels(r,2)=sig.*randn;
        elseif pos_y_t(r,t+1)<=yroom_imp(2)
            pos_y_t(r,t+1)=-pos_y_t(r,t);
            vels(r,2)=sig.*randn;
        end
    end
end

xroom=[-7 7 7 -7 -7];
yroom=[-7 -7 7 7 -7];
x=xroom;
y=yroom;

set(gcf,'nextplot','replacechildren'); 
v = VideoWriter('2d Ideal gas');
v.FrameRate=30;


time=0;
open(v);
for k = 1:tmax 
   time=time+dt
   c = linspace(5,10,length(pos_x_t(:,k)));
   
   
   scatter(pos_x_t(:,k),pos_y_t(:,k),10,c,"filled");
   hold on
   
   plot(x,y,"k")
   
   
   axis([-8 8 -8 8])
   
   
   frame = getframe(gcf);
   writeVideo(v,frame);
   hold off
end

close(v);

set(gcf,'nextplot','replacechildren'); 
vid = VideoWriter('2d Ideal gas velocities');
vid.FrameRate=30;


time=0;
open(vid);
for g = 1:tmax 
   hist(sort(vels_x_t(:,g)),10)
   
   axis([-2000 2000 0 20])
   frame = getframe(gcf);
   writeVideo(vid,frame);
   hold off
end

close(vid);

Q1=zeros(N,tmax);
Q2=zeros(N,tmax);
Q3=zeros(N,tmax);
Q4=zeros(N,tmax);

NQ1=zeros();
NQ2=zeros();
NQ3=zeros();
NQ4=zeros();

for t=1:tmax
    for r=1:N
        if pos_x_t(r,t)>=0&&pos_y_t(r,t)>=0
            Q1(r,t)=1;
        elseif pos_x_t(r,t)>=0&&pos_y_t(r,t)<0
            Q2(r,t)=1;
        elseif pos_x_t(r,t)<0&&pos_y_t(r,t)<0
            Q3(r,t)=1;
        elseif pos_x_t(r,t)<0&&pos_y_t(r,t)>0
            Q4(r,t)=1;
        end
    end
   NQ1(t)=sum(Q1(:,t));
   NQ2(t)=sum(Q2(:,t));
   NQ3(t)=sum(Q3(:,t));
   NQ4(t)=sum(Q4(:,t));
end

Pr1=NQ1/N;
Pr2=NQ2/N;
Pr3=NQ3/N;
Pr4=NQ4/N;

S1=zeros();
S2=zeros();
S3=zeros();
S4=zeros();
time=zeros();
for t=1:tmax
    S1(t)=-kb*Pr1(t)*log(Pr1(t));
    S2(t)=-kb*Pr2(t)*log(Pr2(t));
    S3(t)=-kb*Pr3(t)*log(Pr3(t));
    S4(t)=-kb*Pr4(t)*log(Pr4(t));
    
    time(t+1)=time(t)+dt;
end

S1_max=max(S1);
S2_max=max(S2);
S3_max=max(S3);
S4_max=max(S4);

subplot(2,2,1)
plot(time(1:end-1),S1*1e23,"r")
hold on
plot([time(1) time(end)],[S1_max*1e23 S1_max*1e23],"k","LineStyle","--")
legend({'S(t)','S_{max}'},'Location','southwest',"fontsize",5)
axis([0 0.6 0 0.6])
xlabel("Time [s]")
ylabel("S [*10^{-24}]")
title('Entropy 1st Quadrant')

subplot(2,2,2)
plot(time(1:end-1),S2*1e23,"r")
hold on
plot([time(1) time(end)],[S2_max*1e23 S2_max*1e23],"k","LineStyle","--")
legend({'S(t)','S_{max}'},'Location','southwest',"fontsize",5)
axis([0 0.6 0 0.6])
xlabel("Time [s]")
ylabel("S [*10^{-24}]")
title('Entropy 2nd Quadrant')

subplot(2,2,3)
plot(time(1:end-1),S3*1e23,"r")
hold on
plot([time(1) time(end)],[S3_max*1e23 S3_max*1e23],"k","LineStyle","--")
legend({'S(t)','S_{max}'},'Location','southwest',"fontsize",5)
axis([0 0.6 0 0.6])
xlabel("Time [s]")
ylabel("S [*10^{-24}]")
title('Entropy 3rd Quadrant')

subplot(2,2,4)
plot(time(1:end-1),S4*1e23,"r")
hold on
plot([time(1) time(end)],[S4_max*1e23 S4_max*1e23],"k","LineStyle","--")
legend({'S(t)','S_{max}'},'Location','southwest',"fontsize",5)
axis([0 0.6 0 0.6])
xlabel("Time [s]")
ylabel("S [*10^{-24}]")
title('Entropy 4th Quadrant')
saveas(gcf,'Entropy.png')
close(gcf)


subplot(2,2,1)
scatter(pos_x_t(:,1),pos_y_t(:,1),1,"filled","r")
hold on
plot(x,y,"k")
axis([-8 8 -8 8])
xlabel("X")
ylabel("Y")
title('Position t=0ms')

subplot(2,2,2)
scatter(pos_x_t(:,10),pos_y_t(:,10),1,"filled","r")
hold on
plot(x,y,"k")
xlabel("X")
ylabel("Y")
axis([-8 8 -8 8])
title('Position t=10ms')

subplot(2,2,3)
scatter(pos_x_t(:,100),pos_y_t(:,100),1,"filled","r")
hold on
plot(x,y,"k")
xlabel("X")
ylabel("Y")
axis([-8 8 -8 8])
title('Position t=100ms')


subplot(2,2,4)
scatter(pos_x_t(:,500),pos_y_t(:,500),1,"filled","r")
hold on
plot(x,y,"k")
xlabel("X")
ylabel("Y")
axis([-8 8 -8 8])
title('Position t=500ms')
saveas(gcf,'Position.png')


