%This runs a full optimal power flow study on the caledar year.
clear;
clc;

run = 8760;
date = transpose(linspace(1,run,run));

%Nodal Equation Matrix
A_Matrix = ...
    [0	0	0	0	0	0	0	0	0	1	0	0;
     0	0	0	0	0	0	0	0	1	-1	0	0;
     0	0	0	0	0	0	1	1	-1	0	-1	0;
     0	0	0	0	0	0	0	0	0	0	0	1;
     0	0	0	0	0	0	0	0	0	0	1	-1;
     0	0	0	0	0	1	0	-1	0	0	0	0;
     0	0	0	1	1	-1	-1	0	0	0	0	0;
     0	1	1	0	-1	0	0	0	0	0	0	0;
     1	0	-1	-1	0	0	0	0	0	0	0	0;
     -1	-1	0	0	0	0	0	0	0	0	0	0;
];

%Generator Matrix
G_matrix=...
        [0	0	0
        0	0	0;
        1	0	0;
        0	0	0;
        0	0	0;
        0	1	0;
        0	0	0;
        0	0	0;
        0	0	0;
        0	0	1];

%Line Admittnace
B_matrix = ...
    [-0.0258840704253788;
    -0.236832133383858;
    -0.406295962230727;
    -0.0258768098240791;
    -0.0235423074098;
    -0.0193999667872569;
    -0.0170077745939224;
    -0.123193065708717;
    -0.0227807972458927;
    -0.0165587587554437;
    -0.0415451471113659;
    -0.0147978860331928];

%Bus angle matrix
Theta_matrix = ...
    [1	-1	0	0	0	0	0	0	0	0;
    1	0	-1	0	0	0	0	0	0	0;
    0	1	-1	0	0	0	0	0	0	0;
    0	1	0	-1	0	0	0	0	0	0;
    0	0	1	-1	0	0	0	0	0	0;
    0	0	0	1	-1	0	0	0	0	0;
    0	0	0	1	0	0	0	-1	0	0;
    0	0	0	1	0	0	0	-1	0	0;
    0	0	0	0	0	0	0	1	-1	0;
    0	0	0	0	0	0	0	0	1	-1;
    0	0	0	0	0	-1	0	1	0	0;
    0	0	0	0	0	1	-1	0	0	0];

%Line Loading Limits
Limits_matrix = ...
    [71.3
    32.4;
    34.8;
    30.1;
    21;
    77.5;
    27;
    16.92;
    32.6;
    27;
    27;
    27;];


%Input of Price, Renewable Resource and Load timeseries - hourly
PriceMX = xlsread("ProjectDataTimeseries.xlsx","InputDataHour","B2:B8761");
SolarMX = xlsread("ProjectDataTimeseries.xlsx","InputDataHour","C2:C8761");
WindMX = xlsread("ProjectDataTimeseries.xlsx","InputDataHour","D2:D8761");
LoadMX = xlsread("ProjectDataTimeseries.xlsx","InputDataHour","E1:N8761");
ResultsMX = zeros(run,26);

options = optimoptions('linprog','Display','none');

solarPrice = 30.18;
windPrice = 20.44;

for i = 1:run
    
% decision variables are three generator outputs, 12 line flows and 10 bus
% voltage magnitudes. Line flows and bus voltage have no impact on the cost
% function
f=[solarPrice;
   windPrice;
   PriceMX(i);
   zeros(12,1);
   zeros(10,1)];

%lower bounds
% Minimum generation levels, minimum power flow levels, Maximum angle
lb=[zeros(3,1);
    -1*Limits_matrix; 
    -inf(10,1)];

%upper bounds
% Maximum generation levels, Maximum power flow levels, Maximum angle
Pmax=[SolarMX(i); WindMX(i); 100];
ub=[Pmax;
    Limits_matrix; 
    inf(10,1)];

%equality constraints
Eq1=[G_matrix A_Matrix zeros(10,10)];
Eq2=[zeros(12,3) eye(12) -repmat(B_matrix,1,10).*Theta_matrix ];
Aeq=[Eq1;Eq2];

Demand = ...
    [LoadMX(i,1);
    LoadMX(i,2);
    LoadMX(i,3);
    LoadMX(i,4);
    LoadMX(i,5);
    LoadMX(i,6);
    LoadMX(i,7);
    LoadMX(i,8);
    LoadMX(i,9);
    LoadMX(i,10)];
beq=[Demand;zeros(12,1)];

[x,fval,exitflag,output,lambda] = linprog(f,[],[],Aeq,beq,lb,ub,options);
%[x,fval,exitflag,output,lambda] = linprog(f,[],[],Aeq,beq,lb,ub);
ResultsMX(i,1:25) = transpose(x);
ResultsMX(i,26) = fval;


if exitflag~=1
    disp('No feasible solutions found...')
    return
end
disp(i);

end

totalCost = 0;

for i = 1:run
    totalCost = totalCost + ResultsMX(i,26);
end

%plot results
figure
plot(date,ResultsMX(1:run,1),date,SolarMX(1:run,1))
legend("Dispatched Solar", "Available Solar")
xlabel("Time Period")
ylabel("Generation (MW)")
figure
plot(date,ResultsMX(1:run,2),date,WindMX(1:run,1))
legend("Dispatched Wind", "Available Wind")
xlabel("Time Period")
ylabel("Generation (MW)")

%add a total column to the LoadMX
for i = 1:run
    LoadMX(i,11) = sum(LoadMX(i,1:10));
end

%initialise and indicator matrix for the results
IndiMx = zeros(run,17);


%pull in the ResultsMX and flag if one of the 12 lines has hit a limit by
%comparing the lines load flow to it's limit. Store the loadflow headroom
%in the IndiMx
for i = 1:run
    
   for j = 1:12
       if abs(ResultsMX(i,j+3))==Limits_matrix(j)
           IndiMx(i,4) = 1;
           IndiMx(i,5:16) = transpose(Limits_matrix) - abs(ResultsMX(i,4:15));
       end
   end
    
    
end

i = 0;

for i = 1:run
    %check if possible power is larger than actual generation and no line limits. Indicates constrained
    %energy due to surplus renewable power
    if (SolarMX(i,1)+WindMX(i,1))>sum(ResultsMX(i,1:2)) && IndiMx(i,4) == 0
        IndiMx(i,1) = 1;
    else
        IndiMx(i,1) = 0;
    end
    
    %check if actual from generator 1 and 2 is less than possible power and
    %at least one line limit has been reached but no grid supply
    %Indicates constrained energy due to line limits
    if (SolarMX(i,1)+WindMX(i,1))>sum(ResultsMX(i,1:2)) && IndiMx(i,4) == 1 && ResultsMX(i,3)==0
        IndiMx(i,2) = 1;
    else
        IndiMx(i,2) = 0;    
    end
    
    %check if actual from generator 1 and 2 is less than possible power.
    %at least one line limit requiring grid power purchasing
    %Indicates constrained energy due to line limits
    if (SolarMX(i,1)+WindMX(i,1))>sum(ResultsMX(i,1:2)) && IndiMx(i,4) == 1 && ResultsMX(i,3)>0
        IndiMx(i,3) = 1;
    else
        IndiMx(i,3) = 0;    
    end
    
end


totalLoad = sum(LoadMX(:,11));
totalPossible = sum(SolarMX(:,1)+WindMX(:,1));
totalRenewableGen = sum(sum(ResultsMX(:,1:2)));
totalGridGen = sum(ResultsMX(:,3));

totalWastedExcess = sum(((SolarMX(:,1)+WindMX(:,1))-(ResultsMX(1:run,1)+ResultsMX(1:run,2))).*IndiMx(:,1));
totalWastedLineLimitExcess = sum(((SolarMX(:,1)+WindMX(:,1))-(ResultsMX(1:run,1)+ResultsMX(1:run,2))).*IndiMx(:,2));
totalWastedLineLimitCurtailment = sum(((SolarMX(:,1)+WindMX(:,1))-(ResultsMX(1:run,1)+ResultsMX(1:run,2))).*IndiMx(:,3));
