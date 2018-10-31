function testadapt2


opt = contset(opt,'TestPath',mfilename('fullpath'));
opt = contset(opt,'Filename','testadapt1');

load('Data\testadapt1')
ap = [1];
ID = 2;
opt = contset(opt,'Singularities',               1);  
if(strcmp(s(ID).label ,'H '))
    data = s(ID).data;
    x  = data.x;    
    x1 = x(1:end-1);
    P0 = data.P0;
    p  = [x(end);P0(2)];
    
    h = 1e-6;
    ntst = 20;
    ncol = 4;
    %
    [x0,v0] = init_H_LC_L(@adapt22, x1, p, ap, h, ntst, ncol);
    [~, datafile] = contL(@limitcycleL,x0,v0,opt);
end


[x1,v1]=init_PD_LC(@adaptx,xlc,slc(4),40,4,1e-6);
disp('>> opt=contset(opt,''MaxNumPoints'',250);');
opt = contset(opt,'MaxNumPoints',250);
disp('>> [xlc2,vlc2,slc2,hlc2,flc2]=cont(@limitcycle,x1,v1,opt);');
[xlc2,vlc2,slc2,hlc2,flc2]=cont(@limitcycle,x1,v1,opt);
disp('>> plotcycle(xlc2,vlc2,slc2,[size(xlc2,1) 1 2]);');
figure
axes
plotcycle(xlc2,vlc2,slc2,[size(xlc2,1) 1 2]);
