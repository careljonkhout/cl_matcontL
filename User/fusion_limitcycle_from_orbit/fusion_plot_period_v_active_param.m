function fusion_plot_period_v_active_param(x,x2,a,b)
  figure
  hold on
  plot(x(end,:),x(end-1,:))
  plot(x2(end,:),x2(end-1,:))
  title(sprintf( ...
      'fusion a:%.2f n:%.2f', ...
       a,b));
  xlabel('q_{inf}')
  ylabel('period')
  legend('N=25','N=50')