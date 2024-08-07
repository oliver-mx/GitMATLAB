%% Simulation for optimal length
%
% industrial water price:        2.4  $/m^3
% industrial electricity price:  0.14 $/kWh
%

%X0=zeros(5,20000);
if 1==0
[a,b]=size(X0);
X0(5,:)=linspace(.5,5,b);
for j=1:b
X0(1,:)=62*ones(1,b)+8*rand(1,b);
X0(3,:)=12*ones(1,b)+8*rand(1,b);
X0(2,:)=10*ones(1,b)+(X0(3,:)-10*ones(1,b)+.2).*rand(1,b);
X0(4,:)=1.01*ones(1,b)+1.5*rand(1,b);
end
end

%% calculation
%Y=zeros(1,20000);
c_b1=datetime("now");
for i=1:20000
    i,
    y=fun_1(X0(:,i),5,'Rev',1e4,1e-4); 
    Y(i)=-y;
end

%% stopped at i=12137
c_b2=datetime("now");
save DATA Y c_b1 c_b2
c_b2-c_b1, % = 140:54:14


%% plot
scatter(X0(5,:),Y,'r')












