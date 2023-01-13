// Dynare example based on stochastic growth model, Lecture 3 

// (1) declare endogenous variables

var  c, // consumption 
     k, // capital
     y, // output
     z; // technology                                                   

// (2) declare exogenous variables (shocks)

varexo   e;                               

// (3) declare parameters

parameters alpha, beta, delta, sigma, rho, sigmaeps, 
           z_ss, k_ss, c_ss, y_ss;  

alpha    = 0.35;  // capital share
beta     = 0.95;  // discount factor
delta    = 0.1;   // depreciation rate
rho      = 0.95;  // TFP persistence
sigma    = 3.00;  // CRRA
sigmaeps = 0.01;  // TFP standard deviation

// Import steady-states as parameters
load paramfile
set_param_value('z_ss',z_ss)
set_param_value('k_ss',k_ss)
set_param_value('c_ss',c_ss)
set_param_value('y_ss',y_ss)

// (4) declare the model equations. We want approximation in logs

model; // Model is declared in log variables, i.e. c is log(C), where C is the level of consumption                                   
// Hence exp(c) = C in level

// resource constraint
exp(c) + exp(k) = exp(z)*(exp(k(-1))^alpha)+(1-delta)*exp(k(-1));
                  
// consumption Euler equation
exp(c)^(-sigma) = beta*(exp(c(+1))^(-sigma))*(alpha*exp(z(+1))*(exp(k)^(alpha-1))+1-delta);

// Production function y = z*k(-1)^alpha
exp(y) = exp(z)*exp(k(-1))^alpha;

// law of motion productivity (already in logs, so no need to transform it)
z = rho*z(-1) + e;


end;


// (5) solve the steady state, using analytical results
initval;
    c = log(c_ss);
    k = log(k_ss);
    y = log(y_ss);
    z = log(z_ss);
    e = 0;
end;

steady;

// specify variance of shocks

shocks;
var e = 100*sigmaeps^2;
end;

// (6) solve the dynamics
//stoch_simul(order=1,irf=100);
stoch_simul(order=1,irf=100);
//stoch_simul(order=1,irf=100) c k;
//stoch_simul(order=1,nograph,irf=100);



 

