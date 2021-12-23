

clf();         // close current figure
clear          // clear all pasta variables
xdel(winsid()) // close all windows
A = [  -0.2032588          0  0  0  0  0
        0.2032588  -0.052646  0  0  0  0
                0          0  0  0  0  0
                0          0  0  0  0  0
                0          0  0  0  0  0
                0          0  0  0  0  0];
B = [0.2126574  0
             0  0
             0  0
             0  0
             0  0
             0  0];
C = [1  0  0  0  0  0
     0  1  0  0  0  0
     0  0  1  0  0  0
     0  0  0  1  0  0
     0  0  0  0  1  0
     0  0  0  0  0  1];
D =[0  0
    0  0
    0  0
    0  0
    0  0
    0  0];
    

// Reducing the model
// X=[u, v, r]', U=[tau_u, tau_r]', Y=[u r]'
//los puntos de interrupción para proporcionar si desea un tratamiento lineal por partes.
ap=A(1:3,1:3);  
bp=B(1:3,:); 
cp=[C(1,1:3); 
    C(3,1:3)];
dp=[D(1,1:2);
    D(3,1:2)];
// Controllability and Observability
// Cc=[B, AB, A^2 B,..., A^(n-1) B]
Cc = cont_mat(ap, bp)
rankCc=rank(Cc)
//
// O=[C; CA; CA^2;...; CA^(n-1) ]
O = obsv_mat(ap, cp)
rankO=rank(O)
// verify if the rank of Cc is n, dimension of a
// verify if the rank of O is n, dimension of a

/*             Plot singular values of LTI the model                      */
G = syslin('c', ap, bp, cp, dp);
poles=spec(ap)
tr  = trzeros(G)

w = logspace(-3,3);
sv = svplot(G,w);

scf(0);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude (dB)");


/*                                 Scaling                                 */
d2r=%pi/180;
r2d=180/%pi;

//su = diag( [1/200] );   // scaling input
                                    // tau_u_max=300 [N] tau_r_max=100 [Nm]
//sx = diag([1/1]);  // scaling state 

//sy = diag([1/1]);        // scaling output
su = diag( [1/200,1/(200*0.8)] );   // scaling input
                                    // tau_u_max=300 [N] tau_r_max=100 [Nm]
sx = diag([1/1,1/0.5,1/(50*d2r)]);  // scaling state 

sy = diag([1/1,1/(50*d2r)]);        // scaling output


ap_s = sx*ap*inv(sx)
bp_s = sx*bp*inv(su)
cp_s = sy*cp*inv(sx)
dp_s = sy*dp*inv(su)

//g=minreal(g)
G = syslin('c', ap_s, bp_s, cp_s, dp_s);
w = logspace(-3,3);
sv = svplot(G,w);

scf(1);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude (dB)");
Z=ss2tf(G)
