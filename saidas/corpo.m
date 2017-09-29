
close all; clear all;

n=71;

rho  = load('sondagem_rho.dat');
fase = load('sondagem_fase.dat');

for i=2:10:n+1
    figure(1), semilogy(rho(1,2:end), rho(i,2:end),'*-'); hold on;
    figure(2), plot(fase(1,2:end), fase(i,2:end),'*-'); hold on;
end

if (n>1)
    
figure(3),contourf(rho(1,2:end), log10(rho(2:end,1)),rho(2:end,2:end),30)
colorbar
figure(4),contourf(fase(1,2:end), log10(fase(2:end,1)),fase(2:end,2:end),30)
colorbar

end

