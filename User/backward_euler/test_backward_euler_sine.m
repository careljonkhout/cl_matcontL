A = [0 1; -1 0];

f = @(t,y) A* y;
df = @(t,y) 1;
y0 = [0 1];
t0 = 0;
tf = 20*pi;
n = 10000;


[t,y]=backward_euler(f,df,y0,t0,tf,n);
exp(1)
y(end,:)
plot(t,y);