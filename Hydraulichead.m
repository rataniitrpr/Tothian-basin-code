clear all;
close all;
clc;

L=10000; % length
D=1000; % depth
d=800; % depth of soil layer from the bottom
Hr=50; % regional depth
HL=50; % local depth
P=365; % period
Kx1=1; % conductivity of bottom layer along x
Kz1=1; % conductivity of bottom layer along z
Kx2=1; % conductivity of top layer along x
Kz2=1; % conductivity of top layer along z
Ss1=.0001; % specific storage of bottom layer along
Ss2=.0001; % specific storage of top layer along

tao=L^2*Ss1/(Kx1*365);
Dx1=Kx1/Ss1;
Dx2=Kx2/Ss2;
bta1=(Kz1/Kx1)^.5;
bta2=(Kz2/Kx2)^.5;
An1=(Kz1/Kx1)^.5;
An2=(Kz2/Kx2)^.5;
M=7; % No of local undulations
[x1,z1] = meshgrid(0:100:L,0:10:d);
[x2,z2] = meshgrid(0:100:L,d:10:D);

mloop=21;
mloop1=21;
nloop=999999;
nloop1=21;
dz=.0001;
 
t=0.75*P; % time of simulation
alp=sin(pi*t/P)^2;

%steady-state solution
for m=1:mloop

    Nm=(m-1)*pi/L;

A1=[1 -1 0 0];
aa=Nm*d/An1;
ab=Nm*d/An2;
A2=[exp(aa) exp(-aa)  -exp(ab) -exp(-ab)];
ac=Kz1/An1;
ad=Kz2/An2;
A3=[ac*exp(aa) -ac*exp(-aa)  -ad*exp(ab) ad*exp(-ab)];
ae=Nm*D/An2;
A4=[0 0 exp(ae) exp(-ae)];
A=[A1 
    A2
    A3
    A4];
B=[0 0 0 1]';
C=A\B;
C1(m)=C(1);
C2(m)=C(2);
C3(m)=C(3);
C4(m)=C(4);

if m==1
    Ams(m)=D+Hr;
elseif m==2
    Ams(m)=-Hr;
else
    Ams(m)=0;
end

if m==1
    Am0(m)=D+Hr;
    Amt(m)=D+Hr+(alp*HL);
elseif  m==2
    Am0(m)=-Hr;
    Amt(m)=-Hr;
elseif m==M+1
    Am0(m)=0;
    Amt(m)=-alp*HL;
else
    Am0(m)=0;
    Amt(m)=0;
end

end

s1=0;
s2=0;
for m=1:mloop

    Nm=(m-1)*pi/L;
    Z1=C1(m)*exp(Nm*z1/An1)+C2(m)*exp(-Nm*z1/An1);
    Z2=C3(m)*exp(Nm*z2/An2)+C4(m)*exp(-Nm*z2/An2);

    ba=Ams(m)*cos(Nm*x1).*Z1;
    s1=s1+ba;

    bb=Ams(m)*cos(Nm*x2).*Z2;
    s2=s2+bb;
end
hst1=s1; % 0<z:d
hst2=s2; % d<z:D

% transient solution

for m=1:mloop1
    
    Nm=(m-1)*pi/L;

    s=0;
    for n=1:nloop
        
        x1=(Nm*Nm*Dx1)^.5;
        x2=(Nm*Nm*Dx2)^.5;
        if x1>x2
            x3=x2;
        else
            x3=x1;
        end
        y1=x3+((n-1)*dz);
        y2=x3+(n*dz);
        
        Kz=Kz1/Kz2;
        bta=bta1/bta2;
        
        aa=y1^2/Dx1;
        r11=(aa-Nm^2)^.5;
        ba=y2^2/Dx1;
        r12=(ba-Nm^2)^.5;
        
        bb=y1^2/Dx2;
        r21=(bb-Nm^2)^.5;
        bc=y2^2/Dx2;
        r22=(bc-Nm^2)^.5;
        
        ab=sin(r11*d/bta1);
        ac=sin(r21*(D-d)/bta2);
        g1=Kz*r11*ab*ac;
        ad=cos(r11*d/bta1);
        ae=cos(r21*(D-d)/bta2);
        g2=bta*r21*ad*ae;
        z1=g1-g2;
        
        af=sin(r12*d/bta1);
        ag=sin(r22*(D-d)/bta2);
        g3=Kz*r12*af*ag;
        ah=cos(r12*d/bta1);
        ai=cos(r22*(D-d)/bta2);
        g4=bta*r22*ah*ai;
        z2=g3-g4;
        
        if z1>=0 && z2>=0
            aa=0;
           
        elseif z1<=0 && z2<=0
            aa=0;

        elseif z1>0 && z2<0
            
            aa=abs(z1)+abs(z2);
            ab=dz*abs(z1)/aa;
            ys=y1+ab;
            ab=1;
            s=s+ab;
            Nn(m,s)=ys;
            if s>nloop1
                break
            end
            
        else 
            aa=abs(z1)+abs(z2);
            ab=dz*abs(z1)/aa;
            ys=y1+ab;            
            ab=1;
            s=s+ab;
            Nn(m,s)=ys;
            if s>nloop1
                break
            end
        end
    end
end

for m=1:mloop1
    
    Nm=(m-1)*pi/L;
    for n=1:nloop1+1
    
    a1=(Nn(m,n)^2)/Dx1;
    a2=(Nn(m,n)^2)/Dx2;
    r1=(a1-Nm^2)^.5;
    r2=(a2-Nm^2)^.5;

    ca=Nm/An1;
    cb=Nm/An2;

    f1=@(y)(exp(ca*y).*cos(r1*y/bta1));
    w1=integral(f1,0,d);
    f2=@(y)(exp(-ca*y).*cos(r1*y/bta1));
    w2=integral(f2,0,d);
    ca=Ss1*((w1*C1(m))+(w2*C2(m)));

    f3=@(y)(exp(cb*y).*sin(r2*(D-y)/bta2));
    w3=integral(f3,d,D);
    f4=@(y)(exp(-cb*y).*sin(r2*(D-y)/bta2));
    w4=integral(f4,d,D);

    cb=(w3*C3(m))+(w4*C4(m));
    cc=cos(r1*d/bta1)/sin(r2*(D-d)/bta2);
    cd=Ss2*cb*cc;

    f5=@(y)(cos(r1*y/bta1).*cos(r1*y/bta1));
    w5=integral(f5,0,d);
    
    f6=@(y)(sin(r2*(D-y)/bta2).*sin(r2*(D-y)/bta2));
    w6=integral(f6,d,D);

    ce=Ss1*w5;
    cf=Ss2*w6*cc^2;
    cg=Ams(m)-Am0(m);

    Emno(m,n)=cg*(ca+cd)/(ce+cf);
    Gmn(m,n)=(ca+cd)/(ce+cf);
    end
end

[x1,z1] = meshgrid(0:100:L,0:10:d);
[x2,z2] = meshgrid(0:100:L,d:10:D);
sa=0;
sb=0;
for m=1:mloop1
    
    Nm=(m-1)*pi/L;
    for n=1:nloop1+1
    
    a1=(Nn(m,n)^2)/Dx1;
    a2=(Nn(m,n)^2)/Dx2;
    r1=(a1-Nm^2)^.5;
    r2=(a2-Nm^2)^.5;

    ba=cos(r1*z1/bta1);
    bb=Emno(m,n)*ba.*cos(Nm*x1);
    sa=sa+bb;

    bc=cos(r1*d/bta1)/sin(r2*(D-d)/bta2);
    bd=bc*sin(r2*(D-z2)/bta2);
    be=Emno(m,n)*bd.*cos(Nm*x2);
    sb=sb+be;
    end
end

s1=0;
s2=0;
for m=1:mloop

    Nm=(m-1)*pi/L;
    Z1=C1(m)*exp(Nm*z1/An1)+C2(m)*exp(-Nm*z1/An1);
    Z2=C3(m)*exp(Nm*z2/An2)+C4(m)*exp(-Nm*z2/An2);

    ba=Am0(m)*cos(Nm*x1).*Z1;
    s1=s1+ba;

    bb=Am0(m)*cos(Nm*x2).*Z2;
    s2=s2+bb;
end

ht10=s1+sa; % at t=0
ht20=s2+sb; % at t=0

for m=1:mloop1
    
    Nm=(m-1)*pi/L;
    for n=1:nloop1+1
    
        dd=Nn(m,n);
        ea=pi/P;
        eb=dd^2*sin(2*ea*t);
        ec=2*ea*cos(2*ea*t);
        ed=2*ea*exp(-t*dd^2);
        ee=eb-ec+ed;
        ef=(2*ea)^2+(dd^4);

        if m==1
            fa=HL*ea*ee/ef;
        elseif m==M+1
            fa=-HL*ea*ee/ef;
        else
            fa=0;
        end

        eg=Emno(m,n)*exp(-t*dd^2);
        eh=fa*Gmn(m,n);

        Emnt(m,n)=eg-eh;
    end
end
 
[x1,z1] = meshgrid(0:100:L,0:10:d);
[x2,z2] = meshgrid(0:100:L,d:10:D);
sa=0;
sb=0;
for m=1:mloop1
    
    Nm=(m-1)*pi/L;
    for n=1:nloop1+1
    
    a1=(Nn(m,n)^2)/Dx1;
    a2=(Nn(m,n)^2)/Dx2;
    r1=(a1-Nm^2)^.5;
    r2=(a2-Nm^2)^.5;

    ba=cos(r1*z1/bta1);
    bb=Emnt(m,n)*ba.*cos(Nm*x1);
    sa=sa+bb;

    bc=cos(r1*d/bta1)/sin(r2*(D-d)/bta2);
    bd=bc*sin(r2*(D-z2)/bta2);
    be=Emnt(m,n)*bd.*cos(Nm*x2);
    sb=sb+be;
    end
end

s1=0;
s2=0;
for m=1:mloop

    Nm=(m-1)*pi/L;
    Z1=C1(m)*exp(Nm*z1/An1)+C2(m)*exp(-Nm*z1/An1);
    Z2=C3(m)*exp(Nm*z2/An2)+C4(m)*exp(-Nm*z2/An2);

    ba=Amt(m)*cos(Nm*x1).*Z1;
    s1=s1+ba;

    bb=Amt(m)*cos(Nm*x2).*Z2;
    s2=s2+bb;
end

ht1=s1+sa; % head at time t; bottom layer
ht2=s2+sb; % head at time t; top layer


% boundary head variation; z=D
x3=0:100:L;

h0=D+Hr*(1-cos(pi*x3/L))+(alp*HL*(1-cos(M*pi*x3/L)));

plot(x3,h0)
x=[0:100:L];
z1=[0:10:d];
z2=[d:10:D];
z=[0:10:D];
for i=1:length(z)
    for j=1:length(x)
        if i<=length(z1)
        ht(i,j)= ht1(i,j);
        else
            ht(i,j)=ht2(i-length(z1)+1,j);
        end
    end
end
hhT=ht-D;
htt=hhT/HL;
xx=x/L;
zz=z/L;

figure;

[x,z] = meshgrid(0:100:L,0:10:D);

[C,h]=contour(x,z,ht);











 