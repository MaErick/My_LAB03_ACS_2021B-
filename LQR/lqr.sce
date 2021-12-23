// Form the state space model (assume full state output)
A = [  -0.2032588          0
        0.2032588  -0.052646];
B = [0.2126574    
             0];
C = eye(2,2);
P = syslin("c",A, B, C);
//The compensator weights
Q_xx=diag([20 20]); //Weights on states
R_uu   = 0.2; //Weight on input
Kc=lqr(P,Q_xx,R_uu);

//form the Plant+compensator system

C=[1 0  //dy1
0 1];//dy2
S=C*(P/.(-Kc));
//check system stability
and(real(spec(S.A))<0)
// Check by simulation
dt=0.1;
t=0:dt:30;
u=0.1*rand(t);
y=csim(u,t,S,[1;0;0;0]);
clf;plot(t',y');xlabel(_("time (s)"))
L=legend(["$dy_1$","$dy_2$"]);L.font_size=4;
G=ss2tf(S)
