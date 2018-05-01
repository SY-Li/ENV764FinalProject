clear all; close all; clc

% senarios:
% safe plant(psi50x>psi50s):psi50x=-2;psi50s=-1.5;d=0.01,0.04;psir=0,-1;
% unsafe plant(psi50x<psi50s):psi50x=-1.5;psi50s=-2;d=0.01,0.015,0.02;psir=0,-1;

% parameters setting
gxm=0.0086;gsm=0.47;psi50x=-2;psi50s=-1.5;d=0.04;psir=0;
% ODEs for dynamic system
f=@(t,y)[2*gxm.*(1+(y(2)./psi50x).^5).^(-1).*(psir-2*y(1)+y(2));
    2*gxm.*(1+(y(2)./psi50x).^5).^(-1).*(y(1)-y(2))-gsm.*(1+(y(2)./psi50s).^5).^(-1).*d];

% plot vector field
y1=linspace(-5,0,20);
y2=linspace(-5,0,20);
[x,y]=meshgrid(y1,y2);
u=zeros(size(x));
v=zeros(size(y));

t=0;
for i = 1:numel(x)
    Yprime = f(t,[x(i); y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end
figure
quiver(x,y,u,v,'r'); 
hold on

xlabel('\psix')
ylabel('\psil')
axis tight equal;

% plot one solution solved by numerical method
[ts,ys]=ode45(f,[0,5000],[-1.5;-1.5]);
plot(ys(:,1),ys(:,2),'k')
xlim([-5.0 0])
ylim([-5.0 0])

% plot two nullclines
nullclinex=@(psix,psil)psir-2*psix+psil;
nullclinel=@(psix,psil)2*gxm.*(1+(psil./psi50x).^5).^(-1).*(psix-psil)-gsm.*(1+(psil./psi50s).^5).^(-1).*d;
fimplicit(nullclinex,[-5 0 -5 0],'--b','LineWidth',1.5)
fimplicit(nullclinel,[-5 0 -5 0],'g','LineWidth',1.5)
hold off

% the second plot: how evaporation and leaf water supply change with leaf
% water potential
psil=linspace(-5,0);
e=gsm.*(1+(psil./psi50s).^5).^(-1).*d;
exl=2*gxm.*(1+(psil./psi50x).^5).^(-1).*(psir-psil)*0.5;
figure
plot(psil,e,'--b',psil,exl,'-r')
xlabel('\psil')
ylabel('Fluxes')
