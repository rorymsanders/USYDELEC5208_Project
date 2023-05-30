% This script runs a optimal powerflow and BESS location selector over a one month window using the matrix form of linprog
% The Decision Variables are:
% G1 G2	G3	B8	B9	C8	D8	C9 D9 S109 S108	S98	S97	S87	S76	S73	S63	S32	S21	S35	S54	N10	N9	N8	N7	N6	N5	N4	N3	N2	N1

clear all;
clc;
close all
dur = 744;

%Chose Month asessment type

% %Summer Month
% DemandOri =  readmatrix('ProjectDataTimeseries.xlsx','Sheet','InputDataHour','Range',"E2:N745");

% %Winter Month
% DemandOri =  readmatrix('ProjectDataTimeseries.xlsx','Sheet','InputDataHour','Range',"E4346:N5089");

%Average Year
DemandOri =  readmatrix('ProjectDataTimeseries.xlsx','Sheet','AverageYear','Range',"G2:P745");

%% Equality Constraints
%Nodal Equation Matrix
A_Matrix = ...
    [0	0	0	0	0	0	0	0	0	1	0	0
    0	0	0	0	0	0	0	0	1	-1	0	0
    0	0	0	0	0	0	1	1	-1	0	-1	0
    0	0	0	0	0	0	0	0	0	0	0	1
    0	0	0	0	0	0	0	0	0	0	1	-1
    0	0	0	0	0	1	0	-1	0	0	0	0
    0	0	0	1	1	-1	-1	0	0	0	0	0
    0	1	1	0	-1	0	0	0	0	0	0	0
    1	0	-1	-1	0	0	0	0	0	0	0	0
    -1	-1	0	0	0	0	0	0	0	0	0	0
];

%Generator Matrix.
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

%C8, D8, C9, D9
Batt_matrix=...
    [0	0	0	0
    0	0	0	0
    0	0	0	0
    0	0	0	0
    0	0	0	0
    0	0	0	0
    0	0	0	0
    -1	1	0	0
    0	0	-1	1
    0	0	0	0
];

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

%Batt size consistency matrix
BattSize_matrix = [1 0;
                   0 1];

%Build the equality constraints RHS vector
Demand = [DemandOri zeros(dur,14)];
Demand1 = transpose(Demand);
Demand2 = Demand1(:);


%Construct first itteration of the Aeq matrix
Aeq = [G_matrix zeros(10,2) Batt_matrix A_Matrix zeros(10,10);
       zeros(12,3) zeros(12,2) zeros(12,4) eye(12) -repmat(B_matrix,1,10).*Theta_matrix;
       zeros(2,3) BattSize_matrix zeros(2,26)];

%extrapoloate Aeq matrix over entire period of assessment (744 intervals)
Aeq1Row = size(Aeq,1);
row = size(Aeq,1);
col = size(Aeq,2);
Aeq = repmat({Aeq}, 1, dur);
Aeq = blkdiag(Aeq{:});

%remove battery size contraint functions from initial itteration of the Ae1 maxtrix
Aeq(23:24,4:5) =  [0 0;
                   0 0];

%add the -ve battery size matix

reps = size(Aeq,1)/Aeq1Row;

FstRow = 47;
FstCol = 4;

%capture all necessary equality constraints in Aeq matrix
for i = 1:reps-1

    Aeq((FstRow+(i-1)*row):(FstRow+(i-1)*row+1),(FstCol+(i-1)*col):(FstCol+(i-1)*col+1)) = [-1 0;
                                                                                             0 -1];
end

%Build RHS vector of equality constraints
beq=Demand2;

%% Inequality Matrices

%battery power rating in units of MW
Y = 1;
%battery storage rating in units of MWh
Z = 2;
%battery charge and discharge efficiency
nc = 0.95;
nd = 0.9;
%battery starting charge
E8_0 = 0;
E9_0 = 0;

%inequality matrix formulation
BattSelc_matrix = ...
    [0	0
    -Y	0
    0	0
    0	-Y
    0	0
    -Y	0
    0	0
    0	-Y];
BattCharge_matrix = ...
    [-1	0	0	0
    1	0	0	0
    0	-1	0	0
    0	1	0	0
    0	0	-1	0
    0	0	1	0
    0	0	0	-1
    0	0	0	1];
BattEnergy_matrix = ...
    [-Z	0	nc	-1/nd	0	0
    0	0	-nc	1/nd	0	0
    0	-Z	0	0	nc	-1/nd
    0	0	0	0	-nc	1/nd];

%1st diagonal matrix
Aieq = ...
    [zeros(8,3) BattSelc_matrix BattCharge_matrix zeros(8,22);
     zeros(4,3) BattEnergy_matrix zeros(4,22)];


%diagonilisation for 48 periods
Aieq = repmat({Aieq}, 1, dur);
Aieq = blkdiag(Aieq{:});

%creating overwrite matix
BattEnergy_matrix_sub = ...
    [0	0	nc	-1/nd	0	0
    0	0	-nc	1/nd	0	0
    0	0	0	0	nc	-1/nd
    0	0	0	0	-nc	1/nd];

%including zeros in the overwrite matrix
Aieq_sub = [zeros(8,31);
            zeros(4,3) BattEnergy_matrix_sub zeros(4,22)];

%itteratively placing the overwrite matrix
for j=1:dur-1

    for i=1:(dur-j)
        R1 = (i*12+1+(j-1)*12);
        R2 = (i*12+12+(j-1)*12);
        C1 = ((j-1)*31+1); 
        C2 = ((j-1)*31+31); 
        Aieq((R1):(R2),(C1):(C2)) = Aieq_sub;
        %X = [R1, R2, C1, C2,j,i];
        %disp(X);
    end 

end

bieq =[zeros(8,1);
      -E8_0;
      E8_0;
      -E9_0;
      E9_0;];
bieq = repmat(bieq,dur,1);

%% Upper and Lower Bounds

% Chose month type
% %Summer Month
% %input price data
% QldPrice = readmatrix('ProjectDataTimeseries.xlsx','Sheet','InputDataHour','Range',"B2:B745");
% %input wind and solar data
% SolarGen = readmatrix('ProjectDataTimeseries.xlsx','Sheet','InputDataHour','Range',"C2:C745");
% WindGen = readmatrix('ProjectDataTimeseries.xlsx','Sheet','InputDataHour','Range',"D2:D745");

% %Winter Month
% %input price data
% QldPrice = readmatrix('ProjectDataTimeseries.xlsx','Sheet','InputDataHour','Range',"B4346:B5089");
% %input wind and solar data
% SolarGen = readmatrix('ProjectDataTimeseries.xlsx','Sheet','InputDataHour','Range',"C4346:C5089");
% WindGen = readmatrix('ProjectDataTimeseries.xlsx','Sheet','InputDataHour','Range',"D4346:D5089");

%Average Month
%input price data
QldPrice = readmatrix('ProjectDataTimeseries.xlsx','Sheet','AverageYear','Range',"D2:D745");
%input wind and solar data
SolarGen = readmatrix('ProjectDataTimeseries.xlsx','Sheet','AverageYear','Range',"E2:E745");
WindGen = readmatrix('ProjectDataTimeseries.xlsx','Sheet','AverageYear','Range',"F2:F745");

%line ratings
Limits_matrix = ...
    [71.3;
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
    27];


%upper bounds
Pgen=[0; 0; 100];
ub=[Pgen;
    75;
    136.2;
    inf(4,1);
    Limits_matrix; 
    inf(10,1)];
ub = repmat(ub,1,dur);
ub(1,:) = transpose(SolarGen);
ub(2,:) = transpose(WindGen);
ub = ub(:);

%lower bounds
lb=[zeros(3,1);
    0;
    0;
    -inf(4,1);
    -1*Limits_matrix; 
    -inf(10,1)];
lb = repmat(lb,dur,1);

%% Decision Variables

%Repeating decision matrix construction
%3 Generators, 6 Battery Decisions, 12 line decisions, 10 bus angle
%decisions

solarPrice = 30.18;
windPrice = 20.44;
B8Price = 1.52739726;
B9Price = 1.52739726;

f=[solarPrice;
   windPrice;
   0.15;
   B8Price;
   B9Price;
   0;
   0;
   0;
   0;
   zeros(12,1);
   zeros(10,1)];
f = repmat(f,1,dur);
f(3,:) = transpose(QldPrice);
f = f(:);

%% Optimisation 
%construct optimisation. B8 and B9 are the identified interger variables
[x,fval,exitflag,output] = intlinprog(f,[4 5],Aieq,bieq,Aeq,beq,lb,ub);

if exitflag~=1
    disp('No feasible solutions found...')
    return
end

%% Analysis

%Process outputs
X=reshape(x,31,dur);
Total_cost=fval;

%Determine BESS energy capacity
E8 = zeros(1,length(X));
E9 = zeros(1,length(X));

for i = 1:(length(X))

    if i == 1
        E8(i) = E8_0-(X(7,i)/nd)+X(6,i)*nc;
        E9(i) = E9_0-(X(9,i)/nd)+X(8,i)*nc;
    else
        E8(i) = E8_0+E8(i-1)-(X(7,i)/nd)+X(6,i)*nc;
        E9(i) = E9_0+E9(i-1)-(X(9,i)/nd)+X(8,i)*nc;
    end
    
end

%Add BESS capacity to the results matrix
X = [X(1:5, :); E8; X(6:7, :);E9;X(8:end, :)];
X = transpose(X);

%% Graphs
run = dur;
date = transpose(linspace(1,run,run));

figure
plot(date,X(:,1),date,SolarGen)
legend("Dispatched Solar", "Available Solar")
title('Dispatched Solar vs Available Solar')
xlabel('Interval') 
ylabel('MW')
xlim([0 750])
%saveas(gcf,'C:\Users\RorySanders\OneDrive - The University of Sydney (Students)\2023\ELEC5208\Project\Code\figures\1solar.png');

figure
plot(date,X(:,2),date,WindGen)
legend("Dispatched Wind", "Available Wind")
title('Dispatched Wind vs Available Wind')
xlabel('Interval') 
ylabel('MW') 
xlim([0 750])
%saveas(gcf,'C:\Users\RorySanders\OneDrive - The University of Sydney (Students)\2023\ELEC5208\Project\Code\figures\2wind.png');

figure
plot(date,E8)
legend("Battery_8")
title('CHTO Battery Capacity')
xlabel('Interval') 
ylabel('MWh') 
xlim([0 750])
ylim([0 70])
%saveas(gcf,'C:\Users\RorySanders\OneDrive - The University of Sydney (Students)\2023\ELEC5208\Project\Code\figures\3B8Cap.png');

B8Power = X(:,7:8);
B8Power(:,1) = -1.*B8Power(:,1);
figure
bar(date,B8Power,"stacked")
legend("Charge", "Discharge")
title('B_8 Charge and Discharge Cycle')
xlabel('Interval') 
ylabel('MW') 
xlim([0 750])
ylim([-35 35])
%saveas(gcf,'C:\Users\RorySanders\OneDrive - The University of Sydney (Students)\2023\ELEC5208\Project\Code\figures\4B8CD.png');

figure
plot(date,E9)
legend("Battery_9")
title('MILC Battery Capacity')
xlabel('Interval') 
ylabel('MWh') 
xlim([0 750])
ylim([0 70])
%saveas(gcf,'C:\Users\RorySanders\OneDrive - The University of Sydney (Students)\2023\ELEC5208\Project\Code\figures\5B9CAP.png');

B9Power = X(:,10:11);
B9Power(:,1) = -1.*B9Power(:,1);
figure
bar(date,B9Power,"stacked")
legend("Charge", "Discharge")
title('B_9 Charge and Discharge Cycle')
xlabel('Interval') 
ylabel('MW')
xlim([0 750])
ylim([-35 35])
%saveas(gcf,'C:\Users\RorySanders\OneDrive - The University of Sydney (Students)\2023\ELEC5208\Project\Code\figures\6B9CD.png');

%% Performance analysis

%add a total column to the LoadMX
for i = 1:dur
    DemandOri(i,11) = sum(DemandOri(i,1:10));
end

%initialise and indicator matrix for the results
IndiMx = zeros(dur,17);

%pull in the ResultsMX and flag if one of the 12 lines has hit a limit by
%comparing the lines load flow to it's limit. Store the loadflow headroom
%in the IndiMx
for i = 1:dur

   for j = 1:12
       if abs(X(i,j+11))==Limits_matrix(j)
           IndiMx(i,4) = 1;
           IndiMx(i,5:16) = transpose(Limits_matrix) - abs(X(i,12:23));
       end
   end


end

i = 0;

for i = 1:dur
    %check if possible power is larger than actual generation and no line limits. Indicates constrained
    %energy due to surplus renewable power
    if (SolarGen(i,1)+WindGen(i,1))>sum(X(i,1:2)) && IndiMx(i,4) == 0
        IndiMx(i,1) = 1;
    else
        IndiMx(i,1) = 0;
    end

    %check if actual from generator 1 and 2 is less than possible power and
    %at least one line limit has been reached but no grid supply
    %Indicates constrained energy due to line limits
    if (SolarGen(i,1)+WindGen(i,1))>sum(X(i,1:2)) && IndiMx(i,4) == 1 && X(i,3)==0
        IndiMx(i,2) = 1;
    else
        IndiMx(i,2) = 0;    
    end

    %check if actual from generator 1 and 2 is less than possible power.
    %at least one line limit requiring grid power purchasing
    %Indicates constrained energy due to line limits
    if (SolarGen(i,1)+WindGen(i,1))>sum(X(i,1:2)) && IndiMx(i,4) == 1 && X(i,3)>0
        IndiMx(i,3) = 1;
    else
        IndiMx(i,3) = 0;    
    end

end


totalLoad = sum(DemandOri(:,11))
totalPossible = sum(SolarGen(:,1)+WindGen(:,1))
totalRenewableGen = sum(sum(X(:,1:2)))
totalGridGen = sum(X(:,3))

totalWastedExcess = sum(((SolarGen(:,1)+WindGen(:,1))-(X(:,1)+X(:,2))).*IndiMx(:,1));
totalWastedLineLimitExcess = sum(((SolarGen(:,1)+WindGen(:,1))-(X(:,1)+X(:,2))).*IndiMx(:,2));
totalWastedLineLimitCurtailment = sum(((SolarGen(:,1)+WindGen(:,1))-(X(:,1)+X(:,2))).*IndiMx(:,3));
