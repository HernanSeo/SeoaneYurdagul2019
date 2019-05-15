clear all
clc

load zt_data_arg_logs_series.txt
load zn_data_arg_logs_series.txt
load g_data_arg_logs_series.txt
load sudden_stops.txt

zT=zt_data_arg_logs_series;
zN=zn_data_arg_logs_series;
g=g_data_arg_logs_series;
ss_t=sudden_stops;
tt_estim=1901:2004;

figure(4)
subplot(3,1,1);plot(tt_estim(1:end),zT(1:end,1),'r','LineWidth',1.5);axis('tight');
vline([ss_t],':b');
vline([1914 1931 1959 1982 1995 2001],'b');title('log(z^T): mean reverting shock')
subplot(3,1,2);plot(tt_estim(1:end),zN(1:end,1),'r','LineWidth',1.5);axis('tight');
vline([ss_t],':b');
vline([1914 1931 1959 1982 1995 2001],'b');title('log(z^N): mean reverting shock')
subplot(3,1,3);plot(tt_estim(1:end),g(1:end,1),'r','LineWidth',1.5);axis('tight');
vline([ss_t],':b');
vline([1914 1931 1959 1982 1995 2001],'b');title('log(g): permanent shock')