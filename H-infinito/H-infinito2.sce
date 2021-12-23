// HROV linearized model
// Author: juan C. Cutipa-Luque
// Source:
// x=[u v r x y apsi] 
// tau= [tau1 tau2]
// y= [u v r x y apsi] 
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
Cc = cont_mat(ap,bp) //controbabilidad de la matriz... (ap,bp)) matrices reales de dimensiones adecuadas 
rankCc=rank(Cc) //es el rango numérico de X, es decir, el número de valores singulares de X que son mayores que la norma 
//
// O=[C; CA; CA^2;...; CA^(n-1) ]
O = obsv_mat(ap, cp)
rankO=rank(O)
// verify if the rank of Cc is n, dimension of a
// verify if the rank of O is n, dimension of a

/*             Plot singular values of LTI the model                      */
G = syslin('c', ap, bp, cp, dp); //definición de sistema lineal sad
//syslin define un sistema lineal como una lista y verifica la coherencia de los datos.
poles=spec(ap) //autovalores y autovectores de una matriz o un lápiz
tr  = trzeros(G) //ceros de transmisión y rango normal
//Llamado con un argumento de salida, trzeros (Sl) devuelve los ceros de transmisión del sistema lineal Sl.
w = logspace(-3,3); //vector espaciado logarítmicamente
sv = svplot(G,w); //sigma-plot de valor singular

scf(0);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude (dB)");


/*                                 Scaling                                 */
d2r=%pi/180;
r2d=180/%pi;

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


ms=1.7;// 0.3;%1.5;    % guarantee overshot Mp < 6dB = 20*log10(2) 
wbs=0.23;//0.05;%0.23;
ee=1e-3;//1e-4
ki=1; // used to give more accurate adjustment to the cut-off frequency wbs
      // by default set it to 1
//           --------     WT Data    ------------
mt=1.3;//1.00;    % guarantee overshot Mp < 2dB = 20*log10(1.26)
wbt=4.1;//9.1;%4.1;
ee=1e-3;//1e-4

//           --------     WS     ------------

s=poly(0,'s');
ws1=(s/ms+wbs)/(s+wbs*ee),
ws2=ws1;
ws=[ws1,0;0,ws2]
//Ws=syslin('c',ws)
Ws=blockdiag(ws1,ws2)

//           --------     WT     ------------

s=poly(0,'s');
wt1=(s+wbt/mt)/(ee*s+wbt),
wt2=wt1;
wt=[wt1,0;0,wt2]
//Wt=syslin('c',wt)
Wt=blockdiag(wt1,wt2)


//           --------     WR     ------------
s=poly(0,'s');
wr1=s/s,
wr2=wr1;
wr=[wr1,0;0,wr2]



// ------------------ Plot weighting functions

svs = svplot(Ws,w);
svt = svplot(Wt,w);
scf(2);
plot2d("ln", w, [-20*log(svs')/log(10) -20*log(svt')/log(10)])
xgrid(12)
xtitle("Singular values plot inv(Ws) and inv(Wt)","Frequency (rad/s)", "Amplitude (dB)");



















