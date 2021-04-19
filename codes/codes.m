%% 确定海豚及鱼群的坐标
clc,clear;
dolph_=rand(1,3)*5 %确定海豚初始坐标doph_
dolph=dolph_  %建立海豚坐标随时间变化的矩阵
n=600;             %确定鱼群个数
r=rand(n,1)*1; %随机生成鱼群球坐标（r,theta,fi)
theta=rand(n,1)*2*pi  
fi=rand(n,1)*pi-pi/2
x=r.*cos(theta).*sin(fi);
y=r.*cos(theta);
z=r.*sin(theta);
fish=[x,y,z];      %生成鱼群坐标矩阵
scatter3(x,y,z);   %绘出三维散点图
xlabel('x轴');
xlabel('x轴')
ylabel('y轴')
zlabel('z轴')
grid on
title('鱼群的散点分布图');
hold on
plot3(dolph_(1),dolph_(2),dolph_(3),'r','markersize',50);
Ox=mean(x);
Oy=mean(y);
Oz=mean(z);
% [Ox,Oy,Oz]=[mean(x),mean(y),mean(z)] %确定海豚终点：鱼群的中心
v_dolph=0.5; %假设海豚速度
tmax=sqrt((Ox-dolph_(1))^2+(Oy-dolph_(2))^2+(Oz-dolph_(3))^2)/v_dolph;
for t=1:1:round(tmax)
    dolph=[dolph;dolph(t,:)+v_dolph*[Ox-dolph(t,1),Oy-dolph(t,2),Oz-dolph(t,3)]...
 /(sqrt((Ox-dolph(t,1))^2+(Oy-dolph(t,2))^2+(Oz-dolph(t,3))^2)*t)];
 
% 确定海豚的目标终点
% 定义距离函数d
for i=1:n
    for j=1:n
      d(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2+(z(i)-z(j))^2);  %利用欧氏距离
    end
end
%定义鱼群与海豚初始位置之间的距离
for i=1:n    df_(i)=sqrt((x(i,1)-dolph_(1))^2+(y(i,1)-dolph_(2))^2+(z(i,1)-dolph_(3))^2)';
end
%% 确定每个点的度
r=0.1  %设定半径为1
for i=1:n
    deg1(i)=0;   %设定点的初始度为0
    deg2(i)=i;   
    for j=1:n
        if d(i,j)<r
        deg1(i)=deg1(i)+1;
        deg=[deg1;deg2]';
        end
    end
end
%% 确定虚拟领导集合vlead
nnum=0;
R=3;
vlead=[];
for i=1:n
   if df_(i)<R
      nnum=nnum+1;  %确定虚拟领导的鱼群个数nnum
      vlead=[vlead,i];
   else nnum=nnum;
   end
end
%% 确定每个点的虚拟领导
% 从虚拟领导集合中找出每个点的虚拟领导
p=0.4;   %确定每个点所需虚拟领导占所有虚拟领导的百分比
num=nnum*p;  % 确定每个点的虚拟领导数
degn=sort(deg,'descend');
degnew=degn(1:num,1:2);  
%对虚拟领导集合vlead中元素对应的度进行升序排序
for i=1:n
    for j=1:num
       lead(i,j)=degnew(j,1);     
%lead(i,j)记第i个点的第j个虚拟领导点的度
       leadn(i,j)=degnew(j,2);    
%lead(i,j)记第i个点的第j个虚拟领导点的序号
    end
end
% 确定每个点的邻集内产生的虚拟领导
count=0.05*n  
%每个点至多拥有邻集内的虚拟领导个数为count个（占总数的5%）
time=0        %初始虚拟领导个数为0
leadn=0;      %每个点的领导序号leadn
for j=1:n
    for i=1:n
        if d(i,j)<r
        plus(i)=j;
        else plus=[];
        leadn=[leadn,plus'];  
        time=time+1;
        if time>count
            break
        end
    end
end
end
%% 确定虚拟领导集合中点的运动轨迹
v_fish=1 %设定鱼群游动速度
for i=1:size(vlead)
    x_lead(i)=x(vlead(i));
    y_lead(i)=y(vlead(i));
    z_lead(i)=z(vlead(i));    direct(i,:)=[x_lead-dolph(1),y_lead-dolph(2),z_lead-dolph(3)]... /sqrt((x_lead-dolph(1))^2+(y_lead-dolph(2))^2+(z_lead-dolph(3))^2);
    x_leadn(i)=x_lead(i)+v_fish*direct(1)*t;
    y_leadn(i)=y_lead(i)+v_fish*direct(2)*t;
    z_leadn(i)=z_lead(i)+v_fish*direct(3)*t;
    %虚拟领导鱼群的位置关于时间的函数
end
%% 确定其余点的运动轨迹
% 定义势函数field
for i=1:n
    for j=1:leadn
       field(i,j)=deg(j)/(d(i,j))^3-7200000*(d(i,j))^3;
    end
end
%定义每个点的附近总场之和fields
fields=zeros(n,3); %初始值为n*3的矩阵 行第几个 列xyz
for i=1:n
    for j=1:leadn  %乘以任意两点之间的方向向量
       fields(i,:)=fields(i,:)+field(i,j)*[x(j)-x(i),y(j)-y(i),z(j)-z(i)]...
       /sqrt((x(j)-x(i))^2+(y(j)-y(i))^2+(z(j)-z(i))^2);     
       x=[x,x+v_fish*fields(:,1)*t];
       y=[y,y+v_fish*fields(:,2)*t];
       z=[z,z+v_fish*fields(:,3)*t]; 
%xyz都是n*1的矩阵表示 行表示第几个
%        x=[x,[x_n]'];
%        y=[y,[y_n]'];
%        z=[z,[z_n]'];
    end
end
end
    %任意点的位置等于原来的位置+运动速度*运动方向*时间
 
%% 计算鱼群的被捕概率
eat=0;  %初始被捕鱼数eat
    for i=1:n        df(i)=sqrt((x(i)-dolph(t,1))^2+(y(i)-dolph(t,2))^2+(z(i)-dolph(t,3))^2);
        eat=eat+1;
        p_eat=eat/n;
end
