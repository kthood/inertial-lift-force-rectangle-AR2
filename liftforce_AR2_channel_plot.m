% liftforce_AR2_channel_plot

Re=36; L=45*10^-6; AR=2; rsph=6*10^-6; 
rho=10^3; mu = 10^-1; U = Re*mu/(rho*L);

xs = -.5*AR*L+rsph+.01*L:.05*L:.5*AR*L-rsph-.01*L;
ys = -.5*L+rsph+.01*L:.05*L:.5*L-rsph-.01*L;
[xs,ys]=meshgrid(xs,ys);
forcex = zeros(size(xs)); forcey = zeros(size(ys));
for j=1:size(xs,1)
    for l=1:size(xs,2)
        [forcex(j,l),forcey(j,l)]=getLiftForce_AR2(xs(j,l),ys(j,l),rsph,U,rho,L);
    end
end

figure('Units', 'pixels', 'Position', [100 400 400 200]);
h=10^6;
hold on
quiver(h*xs,h*ys,h*forcex,h*forcey)
set(gca,'XLim',[-.5*AR*L*h,.5*AR*L*h])
set(gca,'YLim',[-.5*L*h,.5*L*h])
xaxis=L*h*(-.5*AR:.1:.5*AR);
plot(xaxis, zeros(size(xaxis)),'k-')
plot(zeros(size(xaxis)),xaxis,'k-')
axis equal
box on
xlabel('x (\mum)')
ylabel('y (\mum)')
title(['Re = ',num2str(Re),', a = ',num2str(rsph/(10^-6)),'\mum'])
