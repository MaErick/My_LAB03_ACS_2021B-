for i=1:length(scs_m.objs)
    if typeof(scs_m.objs(i))=="Block" & scs_m.objs(i).gui=="SUPER_f" then
        scs_m = scs_m.objs(i).model.rpar;
        break;
    end
end
L=[6.696;21.164];
U=[11];
sys = lincos(scs_m,L,U);

A=sys.A
B=sys.B
C=sys.C
D=sys.D

P=syslin('c',A,B,C,D)    //The plant (continuous-time)
Q_xx=diag([1 1]); //Weights on states
R_uu   = 0.5; //Weight on input
Kc=lqr(P,Q_xx,R_uu)

C=[1 0   //dy1
0 1 ];//dy2
S=C*(P/.(-Kc))
//check system stability
and(real(spec(S.A))<0)
// Check by simulation
dt=0.1;
t=0:dt:100;
u=0.1*rand(t);
//y=csim(u,t,S,[6.696;21.164]);
//clf;plot(t',y');xlabel(_("time (s)"))
//L=legend(["$dy_1$","$dy_2$"]);L.font_size=4;

G=ss2tf(S)
nyquist(G)

