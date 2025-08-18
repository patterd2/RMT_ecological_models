x = 0:0.1:6;
a = -0.06;
b = 0.2;
c = 0.5;

figure;
hold on;
for i = 1:100
    plot(x,polyval([unifrnd(-0.2, -0.03) unifrnd(0.1, 0.5) unifrnd(0.1 ,1)], x));
end