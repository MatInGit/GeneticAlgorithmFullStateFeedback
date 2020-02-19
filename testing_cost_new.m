close all; clear all; clc
%Mateusz Leputa @ University of 
N = 1
D = 0
m1=0.01;
m2=0.01;
m3=0.01;
k1=1e4;
k2=1e5;
k3=1e6;
b1=0.1;
b2=0.4;
b3=0.4;
b4=0.4;
F1=5.6e3;
 
m=[m1,0,0;0,m2,0;0,0,m3];
b=[b4+b2+b1,-b2,-b4;-b2,b2+b3,-b3;-b4,-b3,b3+b4];
k=[k1+k3,0,-k3;0,k2,-k2;-k3,-k2,k3+k2];

A = [-b*m^-1,-k*m^-1;eye(3),[0,0,0;0,0,0;0,0,0]];
B = [m^-1;[0,0,0;0,0,0;0,0,0]];
B = [B(:,1)];
C = [[0,0,0;0,0,0;0,0,0],F1*eye(3)];
C = C(1,:)
D=0;
sys1 = ss(A,B,C,D)
[L, P] =lqr(A', C', B*(B'), 0.01);
[K,P]=lqr(A,B,[50,0,0,0,0,0;0,100,0,0,0,0;0,0,100,0,0,0;0,0,0,100,0,0;0,0,0,0,100,0;0,0,0,0,0,100],0.5) %Check 
Acl = (A-B*K)

sys2 = ss(Acl,1*B,C,D)

pop = 1000;
kill_ratio=0.8;
mutation_rate = 0.9;
itterations = 200;

cost_base = costgen1(sys1)
cost_lqr = costgen1(sys2)

%initialize random k matrixes
K_gen =[]
a = -5*max(K);
b = 5*max(K);
parfor i=1:pop
    K_gen = [K_gen;((b-a).*rand(length(K),1) + a)'];
end

cost_history=[];
K_best_gen_j=[]

for j = 1:itterations
     cost_arr = [];
    parfor i=1:pop
        Acl = (A-B*K_gen(i,:));
        sys3 = ss(Acl,B,C,D);
        cost_arr = [cost_arr;costgen1(sys3)];
    end
    %rank by cost
    [cost_arr,I] = sort(cost_arr);
    cost_history = [cost_history, cost_arr(1)];
    K_best_gen_j = [K_best_gen_j;K_gen(I(1),:)]; 
    cost_history(j)
    j
    K_gen;
    K_gen_temp = [];
    for i=1:pop
        K_gen_temp(i,:) = K_gen(I(i),:);
    end
    K_gen = K_gen_temp;
%save best k
    K_gen_temp = [];
    parfor i=1:pop*kill_ratio
        K_gen_temp(i,:) = K_gen(i,:);
    end
%Crossover & Mutation
    parfor i=(pop*kill_ratio+1):pop
        r = randi(length(K)-1);
        r1 = randi(pop*kill_ratio);
        r2 = randi(pop*kill_ratio);
        K_temp = [[K_gen(r1,1:r-1)],[K_gen(r2,r:length(K))]];
        for z = 1:length(K_temp) 
           K_temp(z)= K_temp(z)+mutation_rate*normrnd(0,1)*K_temp(z)*cost_history(j);
        end 
        K_gen_temp = [K_gen_temp;K_temp];

    end
    K_gen = K_gen_temp;
    
end       
    [cost_arr,I] = sort(cost_arr);
    cost_history = [cost_history, cost_arr(1)];
    K_best_gen_j = [K_best_gen_j;K_gen(I(1),:)];  
        Acl = (A-B*K_best_gen_j(j,:));
        sys3 = ss(Acl,B,C,D);
        
time = 0:0.00005:3;


u = sawtooth(time*(2*pi)*5,0.5);
plot(time,u)
hold
plot(time,lsim(sys2,u,time))
plot(time,lsim(sys3,u,time))

grid on 
grid minor 

legend('Base','LQR Controller','Generative K gains')