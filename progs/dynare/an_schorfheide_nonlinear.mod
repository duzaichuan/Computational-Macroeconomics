% This file implements the stationary nonlinear model equations of
% An and Schorfheide (2007, Econometric Reviews)
% =========================================================================
% Willi Mutschler (willi@mutschler.eu)
% Version: June 13, 2023
% =========================================================================

var
c       ${c}$          (long_name='detrended consumption (C/A)')
z       ${z}$          (long_name='shifter to steady-state technology growth')
p       ${\pi}$        (long_name='gross inflation rate')
r       ${R}$          (long_name='nominal interest rate')
y       ${y}$          (long_name='detrended output (Y/A)')
ystar   ${y^*}$        (long_name='detrended flex-price output (Y^*/A)')
g       ${g}$          (long_name='government consumption process (g = 1/(1-G/Y))')
rstar   ${R^*}$        (long_name='auxiliary nominal interest rate for policy rule')
// auxiliary variables for log dev. from steady-state
chat   ${\hat{c}}$     (long_name='detrended consumption (log dev. from steady-state)')
zhat   ${\hat{z}}$     (long_name='technology growth (log dev. from steady-state)')
phat   ${\hat{\pi}}$   (long_name='inflation rate (log dev. from steady-state)')
rhat   ${\hat{R}}$     (long_name='nominal interest rate (log dev. from steady-state)')
yhat   ${\hat{y}}$     (long_name='detrended output (log dev. from steady-state)')
ghat   ${\hat{g}}$     (long_name='government consumption (log dev. from steady-state)')
// observable variables
YGR    ${YGR}$         (long_name='output growth rate (quarter-on-quarter)')
INFL   ${INFL}$        (long_name='annualized inflation rate')
INT    ${INT}$         (long_name='annualized nominal interest rate')
;

varexo
epsr   ${\varepsilon^R}$   (long_name='monetary policy shock')
epsg   ${\varepsilon^g}$   (long_name='government spending shock')
epsz   ${\varepsilon^z}$   (long_name='total factor productivity growth shock')
;

parameters
TAU     ${\tau}$           (long_name='inverse of intertemporal elasticity of subsitution')
KAPPA   ${\kappa}$         (long_name='auxiliary parameter, slope of New-Keynesian Phillips curve')
PSI1    ${\psi_1}$         (long_name='Taylor rule sensitivity parameter to inflation deviations')
PSI2    ${\psi_2}$         (long_name='Taylor rule sensitivity parameter to output deviations')
RHOR    ${\rho_R}$         (long_name='Taylor rule persistence')
RHOG    ${\rho_{g}}$       (long_name='persistence government spending process')
RHOZ    ${\rho_z}$         (long_name='persistence TFP growth rate process')
RA      ${r_{A}}$          (long_name='annualized steady-state real interest rate')
PA      ${\pi^{(A)}}$      (long_name='annualized target inflation rate')
GAMQ    ${\gamma^{(Q)}}$   (long_name='quarterly steady-state growth rate of technology')
NU      ${\nu}$            (long_name='inverse of elasticity of demand in Dixit Stiglitz aggregator')
GBAR    ${\bar{g}}$        (long_name='steady state government spending process')
;


model;
//auxiliary parameters
#GAMMA = 1+GAMQ/100;   //below eq.(38)
#BETA  = 1/(1+RA/400); //below eq.(38)
#PBAR  = 1+PA/400;     //below eq.(38)
#PHI   = TAU*(1-NU)/(NU*PBAR^2*KAPPA); //relationship slope of New Keynesian Phillips Curve eq.(32)

[name='Euler equation']
1 = BETA * (c(+1)/c)^(-TAU) * 1/(GAMMA*z(+1)) * r/p(+1);
[name='Phillips curve based on Rotemberg price setting and Dixit/Stiglitz aggregator']
1 =  PHI*(p-steady_state(p))*( (1-1/(2*NU))*p + steady_state(p)/(2*NU) )
    -PHI*BETA*( (c(+1)/c)^(-TAU)*y(+1)/y*(p(+1)-steady_state(p))*p(+1) )
    +1/NU*(1-c^TAU);
[name='market clearing']
y = c + (g-1)/g*y + PHI/2*(p-steady_state(p))^2*y;
[name='monetary policy equation']
r = rstar^(1-RHOR)*r(-1)^RHOR*exp(epsr);
[name='monetary policy feedback rule']
rstar = steady_state(r) * (p/PBAR)^PSI1 * (y/ystar)^PSI2;
[name='technology growth process']
log(z) = RHOZ*log(z(-1)) + epsz;
[name='government spending process']
log(g) = (1-RHOG)*log(GBAR) + RHOG*log(g(-1)) + epsg;
[name='flex-price output']
ystar = (1-NU)^(1/TAU)*g;

[name='consumption (log dev. from steady-state)']
chat = log(c) - log(steady_state(c));
[name='total factor productivity growth (log dev. from steady-state)']
zhat = log(z) - log(steady_state(z));
[name='inflation (log dev. from steady-state)']
phat = log(p) - log(steady_state(p));
[name='nominal interest rate (log dev. from steady-state)']
rhat = log(r) - log(steady_state(r));
[name='output (log dev. from steady-state)']
yhat = log(y) - log(steady_state(y));
[name='government spending process (log dev. from steady-state)']
ghat = log(g) - log(steady_state(g));

[name='output growth (q-on-q)']
YGR = GAMQ + 100*( yhat - yhat(-1) + zhat );
[name='annualized inflation']
INFL = PA + 400*phat;
[name='annualized nominal interest rate']
INT = PA + RA + 4*GAMQ + 400*rhat;
end;

steady_state_model;
z = 1;
p = 1+PA/400;
g = GBAR;
r = (1+GAMQ/100) / (1/(1+RA/400)) * p;
rstar = r;
c = (1-NU)^(1/TAU);
y = g*c;
ystar = y;
chat=0; zhat=0; phat=0; rhat=0; yhat=0; ghat=0;
YGR  = GAMQ;
INFL = PA;
INT  = PA + RA + 4*GAMQ;
end;

@#include "an_schorfheide_params_shocks.inc"
steady;
stoch_simul(order=1,nograph);