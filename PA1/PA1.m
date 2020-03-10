a = 2.283;   % Vanderwaal constant
b = 0.04278;
R = 0.08314; % gas constant
A = 4.22061; % antoine constant
B = 516.689;
C = 11.223;
Tc = (8*a)/(27*b*R); % critical temp
Pc = a/(27*(b^2));  %critical pressure
T = linspace(Tc-79,Tc,100); %divison of range of temp in 100 parts
Vol = linspace(0.1,10,100);
Xc = [3*b];% critical volume
Yc = [Pc];

for i = 1:100    
      Psat = (10^(A-(B/(C+T(i))))); % antoine equation
      V = [1 (-Psat*b-R*T(i))/Psat a/Psat -(a*b)/Psat]; %vanderwaal equation
      v = roots(V);
      X1(i) = v(1) ;
      X2(i) = v(3);
      Y(i) = Psat; 
end

%T = 300K Isotherm
for j = 1:100
   P1(j) = ((R*300)/(Vol(j)-b))-(a/(Vol(j)^2));
end

%T = 100K Isotherm
for j = 1:100
   P2(j) = ((R*100)/(Vol(j)-b))-(a/(Vol(j)^2));
end

%T = 150K Isotherm
for j = 1:100
   P3(j) = ((R*150)/(Vol(j)-b))-(a/(Vol(j)^2));
end

%T = 190K Isotherm
for j = 1:100
   P4(j) = ((R*190)/(Vol(j)-b))-(a/(Vol(j)^2));
end

%T = 200K Isotherm
for j = 1:100
   P5(j) = ((R*200)/(Vol(j)-b))-(a/(Vol(j)^2));
end

%T = 170K Isotherm
for j = 1:100
   P6(j) = ((R*170)/(Vol(j)-b))-(a/(Vol(j)^2));
end


plot(X1,Y,'*','MarkerEdgeColor','r','DisplayName','Dome','Linestyle','None');
xlim([-1,10]);
ylim([-1,50]);
xlabel('Volume(L/mol)');
ylabel('Pressure(bar)');
hold on;
plot(X2,Y,'*','MarkerEdgeColor','r','DisplayName','Dome','Linestyle','None');
hold on;
plot(Xc,Yc,'o','MarkerEdgeColor','b','DisplayName','Critical Point');
hold on;
plot(Vol,P1,'DisplayName','Isotherm T=300K');
hold on;
plot(Vol,P2,'DisplayName','Isotherm T=100K');
hold on;
plot(Vol,P3,'DisplayName','Isotherm T=150K');
hold on;
plot(Vol,P4,'DisplayName','Isotherm T=190K');
hold on;
plot(Vol,P5,'DisplayName','Isotherm T=200K');
hold on;
plot(Vol,P6,'DisplayName','Isotherm T=170K');
legend;
