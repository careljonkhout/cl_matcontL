f = @(t,y) y;
df = @(t,y) 1;
y0 = 1;
t0 = 0;
tf = 1;
n = 100000;


[t,y]=backward_euler(f,df,y0,t0,tf,n);
exp(1)
y(end)
plot(t,y);