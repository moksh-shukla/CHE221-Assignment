x = zeros(100,1);
y = zeros(100,1);
x(1) = 0.01;
for i = 2:100
    x(i) = x(i-1)+0.01;
end

p = 0.133;
T = zeros(100,1);
T1sat = Tsat(6.93051,1382.650,159.493,p);
T2sat = Tsat(7.13076,1639.856,164.162,p);

for i = 1:100
    x1 = x(i);
    x2 = 1-x(i);
    T(i) = x(i)*T1sat + (1-x(i))*T2sat;
    P1sat = Antoine(6.93051,1382.650,159.493,T(i));
    P2sat = Antoine(7.13076,1639.856,164.162,T(i));
    L12 = lambda12(T(i));
    L21 = lambda21(T(i));    
    gamma1 = exp(-log(x1+L12*x2)+x2*((L12/(x1+L12*x2))-(L21/(x2+L21*x1))));
    gamma2 = exp(-log(x2+L12*x1)-x1*((L12/(x1+L12*x2))-(L21/(x2+L21*x1))));
    P = x1*gamma1*P1sat + x2*gamma2*P2sat;
    y1 = x1*gamma1*P1sat/P;
    y2 = x2*gamma2*P2sat/P;
    epsilon = 0.0001;
    delT = 0.001;
    while(abs(1-(y1+y2))>epsilon)
        if((y1+y2)>1)
            T(i) = T(i)-delT;
            
        else
            T(i) = T(i)+delT;
        end
        
        L12 = lambda12(T(i));
        L21 = lambda21(T(i));    
        gamma1 = exp(-log(x1+L12*x2)+x2*((L12/(x1+L12*x2))-(L21/(x2+L21*x1))));
        gamma2 = exp(-log(x2+L12*x1)-x1*((L12/(x1+L12*x2))-(L21/(x2+L21*x1))));
        P = x1*gamma1*P1sat + x2*gamma2*P2sat;
        y1 = x1*gamma1*P1sat/P;
        y2 = x2*gamma2*P2sat/P;
    end
    
    y(i) = y1; 
    
end

file = "experimental.txt";
M = readmatrix(file);
Texp = M(:,1);
xexp = M(:,2);
yexp = M(:,3);

figure
plot(x,T,'*','MarkerEdgeColor','r');
xlim([0,1]);
hold on;
plot(y,T,'o','MarkerEdgeColor','b');
hold on;
plot(xexp,Texp,'-ro');
hold on;
plot(yexp,Texp,'-bx');

xlabel('x,y');
ylabel('Temp.(K)');
legend('t-cal','y-calc','t-exp','y-exp');
hold off;

figure
plot(x,y);
xlim([0,1]);
hold on;
plot(xexp,yexp,'-ro');
xlabel('x');
ylabel('y');
legend('y-calc','y-exp');
hold off;

function [psat] = Antoine(A,B,C,T)
    psat = 0.133322*(10^(A-(B/(C+(T-273)))));
end

function [T] = Tsat(A,B,C,p)
    T=(B/(A-log10(p*7.50062)))-C+273;
end

function [L12] = lambda12(T)
    v1 = (94.133/1.07)*0.001;
    v2 = (122.16/0.97)*0.001; %https://pubchem.ncbi.nlm.nih.gov/compound/3%2C5-dimethylphenol#section=Density
    R = 1.987; % all pressure units in bar, T: Kelvin
    A12 = -265.7202;
    L12 = (v2/v1)*exp(-A12/(R*T));
end

function [L21] = lambda21(T)
    v1 = (94.133/1.07)*0.001;
    v2 = (122.16/0.97)*0.001; %https://pubchem.ncbi.nlm.nih.gov/compound/3%2C5-dimethylphenol#section=Density
    R = 1.987; % all pressure units in bar, T: Kelvin    
    A21 = 1088.9075;
    L21 = (v1/v2)*exp(-A21/(R*T));
end