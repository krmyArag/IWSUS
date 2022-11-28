clear;clc;
N=1000;
W=20;r=1;h=W*r;u0=3;
z=linspace(0,h,N);
zr=linspace(h,h*1.5,N);
ip=2;
% % %
% *all;
% *u, *v, *uv: U0
% *U,*V: base
% *m*: mean
% res*: res v.s. height
% Res*: weight-sum of every level
% %
%% RESu
umu=0;vmu=0;wmu=0;
ResUWall=0;ResURoad=0;ResURoof=0;
for ip=1:6
    uu=ufitu(ip,W,r,u0,z);vu=vfitu(ip,W,r,u0,z);wu=wfitu(ip,W,r,u0,z); %u0
    umu=umu+mean(abs(uu))/6;vmu=vmu+mean(abs(vu))/6;wmu=wmu+mean(abs(wu))/6;% mean wind speed
    uur=ufitu(ip,W,r,u0,zr);vur=vfitu(ip,W,r,u0,zr);wur=wfitu(ip,W,r,u0,zr); %u0
    %% wall RES
    % u,v
    resuWall=(11.8+4.2*sqrt(uu.^2+vu.^2)).^-1;
    %% road RES
    % u,v,w
    resuRoad=(11.8+4.2*sqrt(uu.^2+vu.^2+wu.^2)).^-1;
    %% roof RES
    % u,v,w
    resuRoof=(11.8+4.2*sqrt(uur.^2+vur.^2+wur.^2)).^-1;
    %% total RES
    ResuWall=0;ResuRoad=0;ResuRoof=0;
    for i=1:N-1
        dh=z(i+1)-z(i);
        dhr=zr(i+1)-zr(i);
        ResuWall=ResuWall+((resuWall(i+1)+resuWall(i))*dh*0.5)/h;
        ResuRoad=ResuRoad+((resuRoad(i+1)+resuRoad(i))*dh*0.5)/h;
        ResuRoof=ResuRoof+((resuRoof(i+1)+resuRoof(i))*dhr*0.5)/(1.5*h-h);
    end
    ResUWall=ResUWall+ResuWall/6;
    ResURoad=ResURoad+ResuRoad/6;
    ResURoof=ResURoof+ResuRoof/6;
    
end
%% MeanSpeedRES
Resum=(11.8+4.2*sqrt(mean(abs(umu))^2+mean(abs(vmu))^2+mean(abs(wmu))^2))^-1;


%% RESv
umv=0;vmv=0;wmv=0;ResV=0;
ResVWall=0;ResVRoad=0;ResVRoof=0;
for ip=1:6
    uv=ufitv(ip,W,r,u0,z);vv=vfitv(ip,W,r,u0,z);wv=wfitv(ip,W,r,u0,z); %u0
    uvr=ufitv(ip,W,r,u0,zr);vvr=vfitv(ip,W,r,u0,zr);wvr=wfitv(ip,W,r,u0,zr); %u0
    if ismember(ip,[1 3 5])
        umv=umv+mean(abs(uv))/3;vmv=vmv+mean(abs(vv))/3;wmv=wmv+mean(abs(wv))/3;
        %% wall RES
        resvWall=(11.8+4.2*sqrt(uv.^2+vv.^2)).^-1;
        %% Total RESwall
        ResvWall=0;
        for i=1:N-1
            dh=z(i+1)-z(i);
            ResvWall=ResvWall+((resvWall(i+1)+resvWall(i))*dh*0.5)/h;
        end
        ResVWall=ResVWall+ResuWall/3;
        ResVRoad=ResVRoad+ResuWall*2/9;
    elseif ismember(ip,[2 4 6])
        %         umv=umv+2*mean(abs(uv))/9;
        %         vmv=vmv+2*mean(abs(vv))/9;
        %         wmv=wmv+2*mean(abs(wv))/9;
        %% road RES
        resvRoad=(11.8+4.2*sqrt(uv.^2+vv.^2+wv.^2)).^-1;
        %% roof RES
        resvRoof=(11.8+4.2*sqrt(uvr.^2+vvr.^2+wvr.^2)).^-1;
        %% Total RESroad
        ResvRoad=0;ResvRoof=0;
        for i=1:N-1
            dh=z(i+1)-z(i);
            dhr=zr(i+1)-zr(i);
            ResvRoad=ResvRoad+((resvRoad(i+1)+resvRoad(i))*dh*0.5)/h;
            ResvRoof=ResvRoof+((resvRoof(i+1)+resvRoof(i))*dhr*0.5)/(1.5*h-h);
        end
        ResVRoad=ResVRoad+ResvRoad/9;
        ResVRoof=ResVRoof+ResvRoof/3;
    else
        disp('no ip. check!')
    end
    
end
Resvm=(11.8+4.2*sqrt(mean(abs(umv))^2+mean(abs(vmv))^2+mean(abs(wmv))^2))^-1;

%% Masson wind
[utop, ucan,wcan]=masson(h,u0,r);
%wcan=mean(abs(wmu))^2;
ResMasson=1/(11.8+4.2*sqrt(ucan^2+wcan^2))

%% Average RES for UV
RESWall=(ResUWall+ResVWall)/2
RESRoad=(ResURoad+ResVRoad)/2
RESRoof=(ResURoof+ResVRoof)/2
RESc=RESWall*2/3+RESRoad/3
RESm=(Resum+Resvm)/2

%% H-obs
data = importfile('H_obs.txt', 2, 25);
H_obs=data.VarName2;
%% FLUX
ts=50:100:2350;
data = importfile('TRf.txt', 2, 25);
Troof=data.VarName2;
data = importfile('Troad2.txt', 2, 25);
Troad=data.VarName2;
data = importfile('Teast.txt', 2, 25);
Teast=data.VarName2;
data = importfile('Twest.txt', 2, 25);
Twest=data.VarName2;
data = importfile('T_canyon2.txt', 2, 25);
Tcan=data.VarName2;
data = importfile('T2_ambient.txt', 2, 25);
Tamb=data.VarName2;

rho=1.2;cp=1004.5;

% Heteo
HWall=.5*(Twest+Teast-2*Tcan)/RESWall;
HRoad=(Troad-Tcan)/RESRoad;
HRoof=(Troof-Tamb)/RESRoof;
Hc=0.4*HWall+0.2*HRoad+0.4*HRoof;

% Masson
%Hm=(0.2*(Troad-Tcan)+0.6*(Twest+Teast-2*Tcan)+0.2*(Troof-Tamb))/RESm;
% Hm=(Troad-Tcan+Twest+Teast-2*Tcan)/RESm;
Hm=(0.2*(Troad-Tcan)+0.4*.5*(Twest+Teast-2*Tcan)+0.4*(Troof-Tcan))/ResMasson;

%% Temp daily plot
figure1=figure(455)
ls=plot(ts,Troad,'r',ts,Tcan,'b',ts,Twest,'c',ts,Teast,'g',ts,Troof,'m',ts,Tamb,'k',ts,0.5*Teast+0.5*Twest,'y')
for i=1:6
    ls(i).LineWidth=1.5;
end
xlim([0,2401])
xticks([0,600,1200,1800,2400])
xticklabels({'0000','0600','1200','1800','2400'})
ylabel('Temperature (\circC)')
xlabel('LST')
ax=ls.Parent;
ax.FontName='Arial';
ax.FontSize=12;

lg=legend('Road','canyon air','east wall','west wall','roof','ambient');
lg.Box='off';
lg.FontSize=13;
lg.Location='northwest';


%% IWSUS-MASSON-obs
figure2=figure(4444)
% subplot(121)
hold on
ls=plot(ts,Hc,'r',ts,Hm,'b');%,ts,HWall,'g',ts,HRoad,'k',ts,HRoof,'y')
ss=scatter(ts,H_obs,'.k');
ss.SizeData=200;

for i=1:2
    ls(i).LineWidth=1.5;
end

xlim([0,2401])
xticks([0,600,1200,1800,2400])
xticklabels({'0000','0600','1200','1800','2400'})
ylabel('Sensible heat flux (W/m^2)')
xlabel('LST')
ax=ls.Parent;
ax.FontName='Arial';
ax.FontSize=12;
ax.Box='on';
ax.LineWidth=1.5;
ll=line([0 2401],[0 0]);
ll.LineStyle=':';
ll.LineWidth=2;
ll.Color='k';
lg=legend('IWSUS','Masson','Obs');
lg.Location='northwest';
lg.FontSize=13;
lg.Box='off';

% annotation(figure2,'textbox',...
%     [0.402429149797571 0.804714284578961 0.0441295554666867 0.071428572563898],...
%     'String',{'(a)'},...
%     'FontSize',13,...
%     'FontName','Arial',...
%     'EdgeColor','none');

%% FARNS roof&roof+canyon
figure(23333)
% subplot(122)
ls2=plot(ts,0.8*HRoof+0.1*HRoad+0.1*HWall,'r',ts,0.4*HRoof+0.2*HRoad+0.4*HWall,'b');

for i=1:2
    ls2(i).LineWidth=1.5;
end

xlim([0,2401])
xticks([0,600,1200,1800,2400])
xticklabels({'0000','0600','1200','1800','2400'})
ylabel('H (W/m^2)')
xlabel('LST')
ax=ls.Parent;
ax.FontName='Arial';
ax.FontSize=12;
ax.Box='on';
ax.LineWidth=1.5;
ll=line([0 2401],[0 0])
ll.LineStyle=':';
ll.LineWidth=2;
ll.Color='k';
lg=legend('roof','roof+canyon');
lg.FontSize=13;
lg.Box='off';

% annotation(figure2,'textbox',...
%     [0.605668016194332 0.816142856189187 0.0388663974464664 0.0638095247631982],...
%     'String','(b)',...
%     'FontSize',13,...
%     'FontName','Arial',...
%     'FitBoxToText','off',...
%     'EdgeColor','none');

subplot(122)

%% Masson horizontal
function [utop, ucan,w]=masson(h,ua,r)
dz=10;
z0=h/10;
cd=1.85;
utop=2*ua*log(h/3/z0)/(pi*log((dz+h/3)/z0));
ucan=utop*exp(-r/2);
w=sqrt(cd)*ua;
end

%% u
%% u fitting functions
function [u]=ufitu(ipos,W,ratio,u0,z)
H=W*ratio;
switch ipos
    case 5 % en_b
        %% log
        pos1=H*(.2781*exp(-0.8123*ratio));
        [~,z1]=min(abs(z-pos1));
        
        if ratio<=0.5
            fr1=(-0.5384*exp(0.5*ratio)+0.5304);
        elseif ratio<=1
            fr1=0.1944*ratio-0.2586;
        else
            fr1=0.09337*exp(-1.087*ratio)-0.09564;
        end
        u1=u0*fr1;
        a1=u1/log(pos1+1);
        u=a1*log(z+1);
        %% increasing
        pos2=H;
        [~,z2]=min(abs(z-pos2));
        
        u2=u0*(0.2887*exp(-1.124*ratio)-0.1473);
        b2=-pos1;c2=u1;a2=(u2-c2)/(pos2+b2)^2;
        u(z1:z2)=a2*(z(z1:z2)+b2).^2+c2;
        
        %% pos3
        pos3=H*(1.712*exp(-5.94*ratio)+1.612);
        [~,z3]=min(abs(z-pos3));
        u3=1.153*u0;
        b3=-pos3;c3=u3;a3=(u2-c3)/(pos2+b3)^2;
        u(z2:z3)=a3*(z(z2:z3)+b3).^2+c3;
        %% exp
        c4=u0;
        b4=-0.1;
        a4=(u3-c4)/exp(b4*pos3);
        u(z3:end)=a4*exp(b4*z(z3:end))+c4;
    case 6 % en_y
        %% log
        pos1=0.1*H;
        [~,z1]=min(abs(z-pos1));
        fr1=-0.3515*exp(-0.9082*ratio)+0.09882;
        u1=fr1*u0;
        a1=u1/log(pos1+1);
        u=a1*log(z+1);
        %% pos2
        pos2=H*(-1.278*exp(-2.643*H)+0.8242);
        if ratio < 0.5
            pos2=H*(-0.8667*ratio+0.9167);
        end
        [~,z2]=min(abs(z-pos2));
        
        u2=u0*(-0.3537*exp(-0.2449*ratio)+0.2514);
        if ratio <0.5
            u2=u0*(-0.4925*ratio+0.1848);
        end
        if ratio>1
            u2=u0*(-0.3537*exp(-0.2449)+0.2514);
        end
        a2=(u2-u1)/(pos2-pos1);
        b2=u2-a2*pos2;
        u(z1:z2)=a2*z(z1:z2)+b2;
        %% pos3
        pos3=H;
        [~,z3]=min(abs(z-pos3));
        u3=u0.*(-0.09157.*exp(1.713*ratio)+0.3768);
        if ratio>1
            u3=u0.*(-0.09157.*exp(1.713)+0.3768);
        end
        b3=-pos2;c3=u2;a3=(u3-u2)/(pos3-pos2)^2;
        u(z2:z3)=a3.*(z(z2:z3)+b3).^2+c3;
        %% pos4
        pos4=H*(0.4612*exp(-4.901*ratio)+1.364);
        [~,z4]=min(abs(z-pos4));
        u4=0.5688*u0;
        b4=-pos3;c4=u3;a4=(u4-c4)/(pos4+b4)^2;
        u(z3:z4)=a4.*(z(z3:z4)+b4).^2+c4;
        %% exp
        c5=u0;
        b5=-0.1;
        a5=(u4-c5)/exp(b5*pos4);
        u(z4:end)=a5*exp(b5*z(z4:end))+c5;
    case 3 % mid_b
        %% log
        pos1=H*(-0.01956*exp(2.563*ratio)+0.3038);
        if ratio>1
            pos1=H*(-0.01956*exp(2.563)+0.3038);
        end
        [~,z1]=min(abs(z-pos1));
        
        u1=u0*(-0.6492*exp(-2.172*ratio)+0.1069);
        a1=u1/log(pos1+1);
        u=a1*log(z+1);
        %% pos2
        pos2=H;
        [~,z2]=min(abs(z-pos2));
        
        fa=@(x) 0.4947*x-0.327;
        fb=@(x) 0.2591*exp(-3.461*ratio)-0.1188;
        if ratio<0.5
            fr2=fa(ratio);
        else
            fr2=fb(ratio);
        end
        u2=fr2*u0;
        
        b2=-pos1;c2=u1;a2=(u2-c2)/(pos2+b2)^2;
        u(z1:z2)=a2*(z(z1:z2)+b2).^2+c2;
        %% pos3
        pos3=H*(1.069*exp(-2.746*ratio)+1.254);
        [~,z3]=min(abs(z-pos3));
        
        u3=u0*(0.4118*exp(-1.276*ratio)+0.719);
        if ratio<0.5
            u3=u0*(0.118*ratio+0.8775);
        end
        b3=-pos2;c3=u2;a3=(u3-c3)/(pos3+b3)^2;
        u(z2:z3)=a3*(z(z2:z3)+b3).^2+c3;
        %% exp
        c4=u0;
        b4=-0.2;
        a4=(u3-c4)/exp(b4*pos3);
        u(z3:end)=a4*exp(b4*z(z3:end))+c4;
    case 4 % mid_y
        %% log
        pos1=H*(0.4667*ratio+0.08333);
        if ratio>0.5
            pos1=H*(4.551*exp(-5.545*ratio)+0.03222);
        end
        [~,z1]=min(abs(z-pos1));
        
        fr1=@(x) 0.01338.*(exp(2.661.*x))-0.1092;
        u1=fr1(ratio).*u0;
        if ratio>=1
            u1=fr1(1).*u0;
        end
        a1=u1/log(pos1+1);
        u=a1*log(z+1);
        %% pos2
        pos2=H*0.8;
        [~,z2]=min(abs(z-pos2));
        
        fr2=@(x) -0.2438.*exp(-3.828*x)-0.01533;
        if ratio<0.5
            fr2=@(x) -0.04493*ratio-0.02887;
        end
        u2=u0*fr2(ratio);
        b2=-pos1;c2=u1;a2=(u2-c2)/((pos2+b2).^2);
        u(z1:z2)=a2.*(z(z1:z2)+b2).^2+c2;
        %% pos3
        fp3=@(x) 1.286.*exp(-6.621.*x)+1.054;
        pos3=H*fp3(ratio);
        if ratio<=0.25
            pos3=H*fp3(0.25);
        end
        [~,z3]=min(abs(z-pos3));
        
        if ratio<=0.5
            fr3=@(x) -0.54007.*exp(0.9084.*ratio)+0.5225;
        else
            fr3=@(x) -0.4343.*exp(-1.968.*ratio)-0.1667;
        end
        u3=u0.*fr3(ratio);
        b3=-pos3;c3=u3;a3=(u2-c3)/((pos2+b3).^2);
        u(z2:z3)=a3.*(z(z2:z3)+b3).^2+c3;
        %% pos4
        pos4=H*(3.289*exp(-5.56*ratio)+1.512);
        [~,z4]=min(abs(z-pos4));
        u4=0.6038*u0;
        b4=-pos3;c4=u3;a4=(u4-c4)/(pos4+b4)^2;
        u(z3:z4)=a4*(z(z3:z4)+b4).^2+c4;
        %% exp
        c5=u0;
        b5=-0.3;
        a5=(u4-c5)/exp(b5*pos4);
        u(z4:end)=a5*exp(b5*z(z4:end))+c5;
    case 1 % out_b
        %% log
        pos1=H*(.2781*exp(-0.8123*ratio));
        [~,z1]=min(abs(z-pos1));
        
        if ratio<=0.5
            fr1=(-0.5384*exp(0.5*ratio)+0.5304);
        elseif ratio<=1
            fr1=0.1944*ratio-0.2586;
        else
            fr1=0.09337*exp(-1.087*ratio)-0.09564;
        end
        u1=u0*fr1;
        a1=u1/log(pos1+1);
        u=a1*log(z+1);
        
        %% increasing
        pos2=H;
        [~,z2]=min(abs(z-pos2));
        
        u2=u0*(0.2887*exp(-1.124*ratio)-0.1473);
        b2=-pos1;c2=u1;a2=(u2-c2)/(pos2+b2)^2;
        u(z1:z2)=a2*(z(z1:z2)+b2).^2+c2;
        %% pos3
        pos3=H*(1.712*exp(-5.94*ratio)+1.612);
        [~,z3]=min(abs(z-pos3));
        u3=1.153*u0;
        b3=-pos3;c3=u3;a3=(u2-c3)/(pos2+b3)^2;
        u(z2:z3)=a3*(z(z2:z3)+b3).^2+c3;
        %% exp
        c4=u0;
        b4=-0.1;
        a4=(u3-c4)/exp(b4*pos3);
        u(z3:end)=a4*exp(b4*z(z3:end))+c4;
    case 2 % out_y
        %% log
        pos1=0.1*H;
        [~,z1]=min(abs(z-pos1));
        fr1=-0.3515*exp(-0.9082*ratio)+0.09882;
        u1=fr1*u0;
        a1=u1/log(pos1+1);
        u=a1*log(z+1);
        %% pos2
        pos2=H*(-1.278*exp(-2.643*H)+0.8242);
        if ratio < 0.5
            pos2=H*(-0.8667*ratio+0.9167);
        end
        [~,z2]=min(abs(z-pos2));
        
        u2=u0*(-0.3537*exp(-0.2449*ratio)+0.2514);
        if ratio <0.5
            u2=u0*(-0.4925*ratio+0.1848);
        end
        if ratio>1
            u2=u0*(-0.3537*exp(-0.2449)+0.2514);
        end
        a2=(u2-u1)/(pos2-pos1);
        b2=u2-a2*pos2;
        u(z1:z2)=a2*z(z1:z2)+b2;
        
        %% pos3
        pos3=H;
        [~,z3]=min(abs(z-pos3));
        u3=u0.*(-0.09157.*exp(1.713*ratio)+0.3768);
        if ratio>1
            u3=u0.*(-0.09157.*exp(1.713)+0.3768);
        end
        b3=-pos2;c3=u2;a3=(u3-u2)/(pos3-pos2)^2;
        u(z2:z3)=a3.*(z(z2:z3)+b3).^2+c3;
        %% pos4
        pos4=H*(0.4612*exp(-4.901*ratio)+1.364);
        [~,z4]=min(abs(z-pos4));
        u4=0.5688*u0;
        b4=-pos3;c4=u3;a4=(u4-c4)/(pos4+b4)^2;
        u(z3:z4)=a4.*(z(z3:z4)+b4).^2+c4;
        %% exp
        c5=u0;
        b5=-0.1;
        a5=(u4-c5)/exp(b5*pos4);
        u(z4:end)=a5*exp(b5*z(z4:end))+c5;
    otherwise
        disp("no position")
end

end

%% v fitting functions
function [v]=vfitu(ipos,W,ratio,u0,z)
H=W*ratio;
switch ipos
    case 5 % en_b
        %% pos1
        pos1=H*(0.4051*exp(-3.592*ratio)+0.03452);
        [~,z1]=min(abs(z-pos1));
        
        v1=u0*(-0.5917*exp(-2.86*ratio)+0.09746);
        if ratio<0.5
            v1=u0*(-0.5917*exp(-2.86*0.5)+0.09746);
        end
        a1=v1/log(pos1+1);
        v=a1*log(z+1);
        
        %% pos2
        pos2=H*(1.079*exp(0.2143*ratio)-0.003025);
        if ratio>1
            pos2=H*(0.7993*exp(-0.3397*ratio)+0.7643);
        end
        [~,z2]=min(abs(z-pos2));
        v2=u0*(-0.06966*exp(-1.58*ratio)-0.1733);
        b2=-pos2;c2=v2;a2=(v1-c2)/(pos1+b2)^2;
        v(z1:z2)=a2*(z(z1:z2)+b2).^2+c2;
        
        %% pos3
        b3=-0.1;
        a3=v2/exp(b3*pos2);
        v(z2:end)=a3*exp(b3*z(z2:end));
    case 6 % en_y
        %% pos1
        pos1=H*(1.519*exp(-4.394*ratio)+0.03125);
        if ratio<0.5
            pos1=H*(0.4*ratio);
        end
        [~,z1]=min(abs(z-pos1));
        
        v1=u0*(-0.02976*exp(2.984*ratio)+0.5034);
        if ratio>1
            v1=u0*(-0.02976*exp(2.984)+0.5034);
        end
        a1=v1/log(pos1+1);
        v=a1*log(z+1);
        
        %% pos2
        pos2=H*(-2.55*exp(0.2608*ratio)+3.722);
        
        if ratio>0.5
            pos2=H*(-2.55*exp(0.2608*0.5)+3.722);
        end
        if ratio<0.25
            pos2=H;
        end
        [~,z2]=min(abs(z-pos2));
        
        v2=0.1498*u0;
        b2=-pos1;c2=v1;a2=(v2-c2)/(pos2+b2)^2;
        v(z1:z2)=a2*(z(z1:z2)+b2).^2+c2;
        
        %% pos3
        pos3=H*(-2.391*exp(0.545*ratio)+4.673);
        if ratio>0.5
            pos3=H*(-2.391*exp(0.545*0.5)+4.673);
        end
        [~,z3]=min(abs(z-pos3));
        v3=u0*(0.1326*exp(-2.765*ratio)-0.1007);
        b3=-pos3;c3=v3;a3=(v2-c3)/(pos2+b3)^2;
        v(z2:z3)=a3*(z(z2:z3)+b3).^2+c3;
        
        %% pos4
        c4=0;
        b4=-0.1;
        a4=v3/exp(b4*pos3);
        v(z3:end)=a4*exp(b4*z(z3:end));
    case 3 % mid_b
        %% pos1
        pos1=H*(1.519*exp(-4.394*ratio)+0.03125);
        if ratio<0.5
            pos1=H*(0.0486*exp(2.811*ratio)+0.001879);
        end
        [~,z1]=min(abs(z-pos1));
        
        v1=u0*(3.516*exp(-1.421*ratio)-1.378);
        a1=v1/log(pos1+1);
        v=a1*log(z+1);
        
        %% pos2
        pos2=H*(-0.3592*exp(-3.522*ratio)+0.9516);
        [~,z2]=min(abs(z-pos2));
        
        v2=u0*(-2.381*exp(-2.057*ratio)+0.9407);
        if ratio < 0.5
            v2=u0*(-1.422*ratio+0.8006);
        end
        b2=-pos1;c2=v1;a2=(v2-c2)/(pos2+b2)^2;
        v(z1:z2)=a2*(z(z1:z2)+b2).^2+c2;
        
        %% 💔pos3
        c3=u0;
        b3=-0.3;
        a3=(v2-c3)/exp(b3*pos2);
        v(z2:end)=a3*exp(b3*z(z2:end))+c3;
        v(:)=0;
    case 4 % mid_y
        %% log
        pos1=3;
        [~,z1]=min(abs(z-pos1));
        
        fr1=-1.528*exp(-1.452*ratio)+2.203;
        if ratio<0.25
            fr1=-1.528*exp(-1.452*0.25)+2.203;
        end
        v1=u0*fr1;
        a1=v1/log(pos1+1);
        v=a1*log(z+1);
        
        %% pos2
        pos2=H;
        [~,z2]=min(abs(z-pos2));
        
        v2=u0*(-0.7184*exp(-1.351*ratio)+1.635);
        if ratio<0.2
            v2=u0*(-0.7184*exp(-1.351*0.2)+1.635);
        end
        b2=-pos1;c2=v1;a2=(v2-c2)/(pos2+b2)^2;
        v(z1:z2)=a2*(z(z1:z2)+b2).^2+c2;
        
        %% pos3
        c3=u0; % v of inlet
        b3=-0.1;
        a3=(v2-c3)/exp(b3*pos2);
        v(z2:end)=a3*exp(b3*z(z2:end))+c3;
        v(:)=0;
    case 1 % out_b
        %% pos1
        pos1=H*(0.4051*exp(-3.592*ratio)+0.03452);
        [~,z1]=min(abs(z-pos1));
        
        v1=-u0*(-0.5917*exp(-2.86*ratio)+0.09746);
        if ratio<0.5
            v1=-u0*(-0.5917*exp(-2.86*0.5)+0.09746);
        end
        a1=v1/log(pos1+1);
        v=a1*log(z+1);
        
        %% pos2
        pos2=H*(1.079*exp(0.2143*ratio)-0.003025);
        if ratio>1
            pos2=H*(0.7993*exp(-0.3397*ratio)+0.7643);
        end
        [~,z2]=min(abs(z-pos2));
        v2=-u0*(-0.06966*exp(-1.58*ratio)-0.1733);
        b2=-pos2;c2=v2;a2=(v1-c2)/(pos1+b2)^2;
        v(z1:z2)=a2*(z(z1:z2)+b2).^2+c2;
        
        %% pos3
        b3=-0.1;
        a3=v2/exp(b3*pos2);
        v(z2:end)=a3*exp(b3*z(z2:end));
    case 2 % out_y
        %% pos1
        pos1=H*(1.519*exp(-4.394*ratio)+0.03125);
        if ratio<0.5
            pos1=H*(0.4*ratio);
        end
        [~,z1]=min(abs(z-pos1));
        
        v1=-u0*(-0.02976*exp(2.984*ratio)+0.5034);
        if ratio>1
            v1=u0*(-0.02976*exp(2.984)+0.5034);
        end
        a1=v1/log(pos1+1);
        v=a1*log(z+1);
        
        %% pos2
        pos2=H*(-2.55*exp(0.2608*ratio)+3.722);
        
        if ratio>0.5
            pos2=H*(-2.55*exp(0.2608*0.5)+3.722);
        end
        if ratio<0.25
            pos2=H;
        end
        [~,z2]=min(abs(z-pos2));
        
        v2=-0.1498*u0;
        b2=-pos1;c2=v1;a2=(v2-c2)/(pos2+b2)^2;
        v(z1:z2)=a2*(z(z1:z2)+b2).^2+c2;
        
        %% pos3
        pos3=H*(-2.391*exp(0.545*ratio)+4.673);
        if ratio>0.5
            pos3=H*(-2.391*exp(0.545*0.5)+4.673);
        end
        [~,z3]=min(abs(z-pos3));
        v3=-u0*(0.1326*exp(-2.765*ratio)-0.1007);
        b3=-pos3;c3=v3;a3=(v2-c3)/(pos2+b3)^2;
        v(z2:z3)=a3*(z(z2:z3)+b3).^2+c3;
        
        %% pos4
        c4=0;
        b4=-0.1;
        a4=v3/exp(b4*pos3);
        v(z3:end)=a4*exp(b4*z(z3:end));
    otherwise
        disp("no position")
end
end

%% w fitting functions
function [w]=wfitu(ipos,W,ratio,u0,z)
H=W*ratio;
switch ipos
    case 5 % en_b
        %% poly2
        fp1=@(x) -5.284*exp(-6.056*x)+0.9939;
        pos1=H*fp1(ratio);
        if ratio<0.5
            pos1=H*fp1(0.5);
        end
        [~,z1]=min(abs(z-pos1));
        w1=u0*(0.04559*exp(2.746*ratio)+0.002081);
        if ratio>=0.5
            w1=u0*(0.3718*exp(-5.31*ratio)+0.1556);
        end
        b1=-pos1;c1=w1;a1=-c1/(b1)^2;
        w=a1*(z+b1).^2+c1;
        %% linear
        pos2=H*(0.3585*exp(-1.898*ratio)+1.057);
        [~,z2]=min(abs(z-pos2));
        
        w2=u0*(-0.307*exp(-5.844*ratio)+0.1185);
        a2=(w2-w1)/(pos2-pos1);b2=w2-a2*pos2;
        w(z1:z2)=a2*z(z1:z2)+b2;
        %% pos3 poly2
        pos3=H*(1.475*exp(-3.011*ratio)+1.318);
        [~,z3]=min(abs(z-pos3));
        
        w3=u0*(-0.426*exp(-1.736*ratio)+0.4536);
        a3=(w3-w2)/(pos3-pos2);b3=w3-a3*pos3;
        w(z2:z3)=a3*z(z2:z3)+b3;
        %% exp
        c4=u0*0.05;
        b4=-0.1;
        a4=(w3-c4)/exp(b4*pos3);
        w(z3:end)=a4*exp(b4*z(z3:end))+c4;
    case 6 % en_y
        %% poly2
        pos1=H*(-3.712*exp(-5.28*ratio)+0.9836);
        if ratio <0.5
            pos1=-3.712*exp(-5.28*0.5)+0.9836;
        end
        [~,z1]=min(abs(z-pos1));
        
        w1=u0*(-0.04245*exp(4.028*ratio)-0.002003);
        if ratio>0.5
            w1=u0*(-2.098*exp(-5.081*ratio)-0.1547);
        end
        % b1=-pos1;c1=w1;a1=-c1/b1^2;
        b1=0;c1=0;a1=w1/(pos1)^2;
        w=a1*(z+b1).^2+c1;
        %% pos2
        pos2=H*(8.236*exp(-5.263*ratio)+1.857);
        [~,z2]=min(abs(z-pos2));
        
        w2=u0*(-0.7792*exp(-2.273*ratio)+0.2835);
        if ratio<0.5
            w2=u0*(-0.7792*exp(-2.273*0.5)+0.2835);
        end
        b2=-pos2;c2=w2;a2=(w1-c2)/(pos1+b2)^2;
        w(z1:z2)=a2*(z(z1:z2)+b2).^2+c2;
        
        %% pos3
        w(z2:end)=w2;
    case 3 % mid_b
       %% pos1
        pos1=H*0.9;
        [~,z1]=min(abs(z-pos1));
        
        w1=u0*(0.5417*ratio+0.03767);
        if ratio >=0.5
            w1=u0*(9.341*exp(-8.185*ratio)+0.1524);
        end
        b1=-pos1;c1=w1;a1=-c1/(b1)^2;
        w=a1*(z+b1).^2+c1;
        %% pos2
        pos2=H*(0.5214*exp(-1.5*ratio)+1.055);
        [~,z2]=min(abs(z-pos2));
        
        w2=u0*(0.225*ratio+0.05438);
        if ratio > 0.5
            w2=u0*(0.904*exp(-5.519*ratio)+0.1201);
        end
        a2=(w2-w1)/(pos2-pos1);b2=w2-a2*pos2;
        w(z1:z2)=a2*z(z1:z2)+b2;
        %% pos3
        pos3=(1.534*exp(-2.733*ratio)+1.337)*H;
        [~,z3]=min(abs(z-pos3));
        
        w3=u0*(-0.4408*exp(-2.875*ratio)+0.5083);
        a3=(w3-w2)/(pos3-pos2);b3=w3-a3*pos3;
        w(z2:z3)=a3*z(z2:z3)+b3;
        %% pos4
        c4=u0/10;
        b4=-0.1;
        a4=(w3-c4)/exp(b4*pos3);
        w(z3:end)=a4*exp(b4*z(z3:end))+c4;
    case 4 % mid_y
        %% pos1
        pos1=H*(-0.4067*ratio+0.7757);
        if ratio >= 0.5
            pos1=H*(0.0219*exp(1.688*ratio)+0.5214);
        end
        [~,z1]=min(abs(z-pos1));
        
        w1=0.2472*ratio+0.2858;
        if ratio >= 0.5
            w1=4.77*exp(-5.175*ratio)+0.05055;
        end
        w1=w1*u0;
        b1=-pos1;c1=w1;a1=-c1/b1^2;
        w=a1*(z+b1).^2+c1;
        %% pos2
        pos2=H*(4.424*exp(-4.943*ratio)+1.005);
        [~,z2]=min(abs(z-pos2));
        
        w2=u0*(0.4383*exp(-4.059*ratio)-0.1127);
        a2=(w2-w1)/(pos2-pos1);
        b2=w2-a2*pos2;
        w(z1:z2)=a2*z(z1:z2)+b2;
        %% pos3
        pos3=H*(2.654*exp(-2.872*ratio)+1.73);
        [~,z3]=min(abs(z-pos3));
        
        w3=u0*(-1.275*exp(-0.6783*ratio)+0.9835);
        b3=-pos2;c3=w2;a3=(w3-c3)/(pos3+b3)^2;
        w(z2:z3)=a3*(z(z2:z3)+b3).^2+c3;
        %% pos4
        w(z3:end)=w3;
    case 1 % out_b
         %% poly2
        fp1=@(x) -5.284*exp(-6.056*x)+0.9939;
        pos1=H*fp1(ratio);
        if ratio<0.5
            pos1=H*fp1(0.5);
        end
        [~,z1]=min(abs(z-pos1));
        w1=u0*(0.04559*exp(2.746*ratio)+0.002081);
        if ratio>=0.5
            w1=u0*(0.3718*exp(-5.31*ratio)+0.1556);
        end
        b1=-pos1;c1=w1;a1=-c1/(b1)^2;
        w=a1*(z+b1).^2+c1;
        %% linear
        pos2=H*(0.3585*exp(-1.898*ratio)+1.057);
        [~,z2]=min(abs(z-pos2));
        
        w2=u0*(-0.307*exp(-5.844*ratio)+0.1185);
        a2=(w2-w1)/(pos2-pos1);b2=w2-a2*pos2;
        w(z1:z2)=a2*z(z1:z2)+b2;
        %% pos3 poly2
        pos3=H*(1.475*exp(-3.011*ratio)+1.318);
        [~,z3]=min(abs(z-pos3));
        
        w3=u0*(-0.426*exp(-1.736*ratio)+0.4536);
        a3=(w3-w2)/(pos3-pos2);b3=w3-a3*pos3;
        w(z2:z3)=a3*z(z2:z3)+b3;
        %% exp
        c4=u0*0.05;
        b4=-0.1;
        a4=(w3-c4)/exp(b4*pos3);
        w(z3:end)=a4*exp(b4*z(z3:end))+c4;
    case 2 % out _y
        %% poly2
        pos1=H*(-3.712*exp(-5.28*ratio)+0.9836);
        if ratio <0.5
            pos1=-3.712*exp(-5.28*0.5)+0.9836;
        end
        [~,z1]=min(abs(z-pos1));
        
        w1=u0*(-0.04245*exp(4.028*ratio)-0.002003);
        if ratio>0.5
            w1=u0*(-2.098*exp(-5.081*ratio)-0.1547);
        end
        % b1=-pos1;c1=w1;a1=-c1/b1^2;
        b1=0;c1=0;a1=w1/(pos1)^2;
        w=a1*(z+b1).^2+c1;
        %% pos2
        pos2=H*(8.236*exp(-5.263*ratio)+1.857);
        [~,z2]=min(abs(z-pos2));
        
        w2=u0*(-0.7792*exp(-2.273*ratio)+0.2835);
        if ratio<0.5
            w2=u0*(-0.7792*exp(-2.273*0.5)+0.2835);
        end
        b2=-pos2;c2=w2;a2=(w1-c2)/(pos1+b2)^2;
        w(z1:z2)=a2*(z(z1:z2)+b2).^2+c2;
        
        %% pos3
        w(z2:end)=w2;
        %% pos2
    otherwise
        disp("no position")
end
end

%% v
%% u fitting functions
function [u]=ufitv(ipos,W,r,u0,z)
H=W*r;
switch ipos
    case 5 % en_c
        u=zeros(1,length(z));
    case 6 % en_y
        %% pos1
        p1=H*(0.7662*exp(-3.718*r)+0.0975);
        [~,z1]=min(abs(z-p1));
        u1=u0*(0.7217*exp(-8.564*r)-0.1092);
        if r>=1
            u1=u0*(0.7217*exp(-8.564*1)-0.1092);
        end
        a1=u1/p1;
        u=a1*z;
        
        %% pos2
        p2=H;
        [~,z2]=min(abs(z-p2));
        u2=u0*(0.02603*exp(-3.252*r)+0.009276);
        
        b2=-p2;c2=u2;a2=(u1-c2)/(p1+b2)^2;
        u=a2*(z+b2).^2+c2;
        
        %% pos3
        c3=0;
        b3=-0.05;
        a3=(u2-c3)/exp(b3*p2);
        u(z2:end)=a3*exp(b3*z(z2:end))+c3;
        %         plot(u(z2:end),z(z2:end))
    case 3 % mid_c
        u=zeros(1,length(z));
    case 4 % mid_y
        %% pos1
        p1=H;
        [~,z1]=min(abs(z-p1));
        u1=u0*(-0.04175*exp(-9.374*r)+0.02567);
        
        b1=-H;c1=u1;a1=-c1/b1^2;
        
        u=a1*(z+b1).^2+c1;
        
        %% pos2
        c3=0;
        b3=-0.05;
        a3=(u1-c3)/exp(b3*p1);
        u(z1:end)=a3*exp(b3*z(z1:end))+c3;
    case 1 % out_c
        u=zeros(1,length(z));
    case 2 % out_y
        %% pos1
        p1=H;
        [~,z1]=min(abs(z-p1));
        u1=u0*(-0.04175*exp(-9.374*r)+0.02567);
        
        b1=-H;c1=u1;a1=-c1/b1^2;
        
        u=0.66*(a1*(z+b1).^2+c1);
        
        %% pos2
        c3=0;
        b3=-0.05;
        a3=(u1-c3)/exp(b3*p1);
        u(z1:end)=0.66*(a3*exp(b3*z(z1:end))+c3);
    otherwise
        disp("no position")
end

end

%% v fitting functions
function [v]=vfitv(ipos,W,r,u0,z)
H=W*r;
switch ipos
    case 5 % en_c
        %% pos1
        p1=H*(0.5697*exp(-1.487*r)+0.049);
        [~,z1]=min(abs(z-p1));
        v1=u0*(-0.205*exp(-1.619*r)+1.192);
        
        a1=v1/log(p1+1);
        v=a1.*log(z+1);
        %% pos2
        c2=u0;
        b2=-0.05;
        a2=(v1-c2)/exp(b2*p1);
        v(z1:end)=a2*exp(b2*z(z1:end))+c2;
    case 6 % en_y
        %% pos1
        p1=H*(0.5948*exp(-1.288*r)+0.1745);
        [~,z1]=min(abs(z-p1));
        v1=u0*(1.05);
        
        a1=v1/log(p1+1);
        v=a1.*log(z+1);
        
        %% pos2
        p2=H;
        [~,z2]=min(abs(z-p2));
        v2=u0*(-0.0534*exp(-3.028*r)+1.102);
        a2=(v2-v1)/(p2-p1);b2=v2-a2*p2;
        v(z1:z2)=a2*z(z1:z2)+b2;
        
        %% pos3
        c3=u0;
        b3=-0.05;
        a3=(v2-c3)/exp(b3*p2);
        v(z2:end)=a3*exp(b3*z(z2:end))+c3;
    case 3 % mid_c
        %% pos1
        p1=H*(0.901*exp(-4.313*r)+0.194);
        [~,z1]=min(abs(z-p1));
        v1=u0*(-0.258*exp(-1.673*r)+1.2);
        
        a1=v1/log(p1+1);
        v=a1.*log(z+1);
        
        %% pos2
        c2=u0;
        b2=-0.05;
        a2=(v1-c2)/exp(b2*p1);
        v(z1:end)=a2*exp(b2*z(z1:end))+c2;
    case 4 % mid_y
        %% pos1
        p1=H*(0.5948*exp(-1.288*r)+0.1745);
        [~,z1]=min(abs(z-p1));
        v1=u0*(0.9);
        
        a1=v1/log(p1+1);
        v=a1.*log(z+1);
        
        %% pos2
        p2=H*(2.206*exp(-3.593*r)+1.173);
        [~,z2]=min(abs(z-p2));
        v2=u0*(-0.4812*exp(-0529*r)+1.505);
        a2=(v2-v1)/(p2-p1);b2=v2-a2*p2;
        v(z1:z2)=a2*z(z1:z2)+b2;
        
        %% pos3
        c3=u0;
        b3=-0.05;
        a3=(v2-c3)/exp(b3*p2);
        v(z2:end)=a3*exp(b3*z(z2:end))+c3;
    case 1 % out_c
        %% pos1
        p1=H*(0.901*exp(-4.313*r)+0.194);
        [~,z1]=min(abs(z-p1));
        v1=u0*(-0.258*exp(-1.673*r)+1.2);
        
        a1=v1/log(p1+1);
        v=0.9*a1.*log(z+1);
        
        %% pos2
        c2=u0;
        b2=-0.05;
        a2=(v1-c2)/exp(b2*p1);
        v(z1:end)=0.9*(a2*exp(b2*z(z1:end))+c2);
    case 2 % out_y
        %% pos1
        p1=H*(0.5948*exp(-1.288*r)+0.1745);
        [~,z1]=min(abs(z-p1));
        v1=u0*(0.9);
        
        a1=v1/log(p1+1);
        v=0.9*a1.*log(z+1);
        
        %% pos2
        p2=H*(2.206*exp(-3.593*r)+1.173);
        [~,z2]=min(abs(z-p2));
        v2=u0*(-0.4812*exp(-0529*r)+1.505);
        a2=(v2-v1)/(p2-p1);b2=v2-a2*p2;
        v(z1:z2)=0.9*(a2*z(z1:z2)+b2);
        
        %% pos3
        c3=u0;
        b3=-0.05;
        a3=(v2-c3)/exp(b3*p2);
        v(z2:end)=0.9*(a3*exp(b3*z(z2:end))+c3);
    otherwise
        disp("no position")
end
end

%% w fitting functions
function [w]=wfitv(ipos,W,r,u0,z)
H=W*r;
switch ipos
    case 5 % en_c
        %% pos1
        p1=H;
        [~,z1]=min(abs(z-p1));
        w1=u0*(-0.156*exp(-3.475*r)+0.09242);
        
        b1=-p1;c1=w1;a1=-c1/b1^2;
        
        w=a1*(z+b1).^2+c1;
        
        %% pos2
        c2=0;
        b2=-0.05;
        a2=(w1-c2)/exp(b2*p1);
        w(z1:end)=a2*exp(b2*z(z1:end))+c2;
    case 6 % en_y
        %% pos1
        p1=H;
        [~,z1]=min(abs(z-p1));
        w1=u0*(-0.1621*exp(-1.648*r)+0.156);
        
        b1=-p1;c1=w1;a1=-c1/b1^2;
        
        w=a1*(z+b1).^2+c1;
        
        %% pos2
        c2=0;
        b2=-0.05;
        a2=(w1-c2)/exp(b2*p1);
        w(z1:end)=a2*exp(b2*z(z1:end))+c2;
    case 3 % mid_c
        %% pos1
        p1=0.8*H;
        [~,z1]=min(abs(z-p1));
        w1=u0*(-0.06537*exp(-2.274*r)+0.04656);
        
        b1=-p1;c1=w1;a1=-c1/b1^2;
        
        w=a1*(z+b1).^2+c1;
        
        %% pos2
        c2=0;
        b2=-0.05;
        a2=(w1-c2)/exp(b2*p1);
        w(z1:end)=a2*exp(b2*z(z1:end))+c2;
    case 4 % mid_y
        %% pos1
        p1=H*(1.222*exp(-3.644*r)+0.02068);
        [~,z1]=min(abs(z-p1));
        w1=u0*(-0.0484*exp(-0.8592*r)+0.00731);
        if r<0.5
            w1=u0*(-0.05433*r+0.003);
        end
        
        b1=-p1;c1=w1;a1=-c1/b1^2;
        
        w=a1*(z+b1).^2+c1;
        
        %% pos2
        p2=H;
        [~,z2]=min(abs(z-p2));
        w2=u0*(-1.034*exp(-0.0328*r)+1.025);
        if r<0.5
            w2=u0*(-1.034*exp(-0.0328*0.5)+1.025);
        end
        b2=-p2;c2=w2;a2=(w1-c2)/(p1+b2)^2;
        
        w(z1:z2)=a2*(z(z1:z2)+b2).^2+c2;
        
        %% pos3
        w(z1:end)=w2;
    case 1 % out_c
        %% pos1
        p1=0.8*H;
        [~,z1]=min(abs(z-p1));
        w1=u0*(-0.06537*exp(-2.274*r)+0.04656);
        
        b1=-p1;c1=w1;a1=-c1/b1^2;
        
        w=0.7*(a1*(z+b1).^2+c1);
        
        %% pos2
        c2=0;
        b2=-0.05;
        a2=(w1-c2)/exp(b2*p1);
        w(z1:end)=0.7*(a2*exp(b2*z(z1:end))+c2);
    case 2 % out _y
        %% pos1
        p1=H*(1.222*exp(-3.644*r)+0.02068);
        [~,z1]=min(abs(z-p1));
        w1=u0*(-0.0484*exp(-0.8592*r)+0.00731);
        if r<0.5
            w1=u0*(-0.05433*r+0.003);
        end
        
        b1=-p1;c1=w1;a1=-c1/b1^2;
        
        w=0.7*(a1*(z+b1).^2+c1);
        
        %% pos2
        p2=H;
        [~,z2]=min(abs(z-p2));
        w2=u0*(-1.034*exp(-0.0328*r)+1.025);
        if r<0.5
            w2=u0*(-1.034*exp(-0.0328*0.5)+1.025);
        end
        b2=-p2;c2=w2;a2=(w1-c2)/(p1+b2)^2;
        
        w(z1:z2)=0.7*(a2*(z(z1:z2)+b2).^2+c2);
        
        %% pos3
        w(z1:end)=w2*0.7;
    otherwise
        disp("no position")
end
end


%% fileread
function Tr = importfile(filename, startRow, endRow)

delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

formatSpec = '%f%f%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

fclose(fileID);

Tr = table(dataArray{1:end-1}, 'VariableNames', {'DataThiefFurban_wind_profile_parameterizationcodewprofwprof_dat','VarName2'});

end
