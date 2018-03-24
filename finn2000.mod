/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This script reproduces Mary G. Finn DSGE model for:
"perfect competition and the effects of energy price increases
on economic activity"
(2000), Journal of Money credit and Banking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Declares the endogenous variables 13;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var y ${y}$ (long_name='Production')
    l ${l}$ (long_name='Labor')
    k ${k}$ (long_name='Capital')
    u ${u}$ (long_name='Utilization rate of capital')
    delt ${\delta}$ (long_name='depreciation rate of capital')
    i ${i}$ (long_name='Investment')
    a ${a}$ (long_name='Technical relationship')
    e ${e}$ (long_name='Energy usage')
    c ${c}$ (long_name='Consumption')
    r ${r}$ (long_name='Rental income from capital services')
    w ${w}$ (long_name='Rental income from labor services')
    z ${z}$ (long_name='Technology shock')
    p ${p}$ (long_name='Energy price')
    
;
    

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//declares the exogenous variables:

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo 
    eps_z 
    eps_p
;

parameters //12

    alph ${\alpha}$ (long_name='Preference parameter')
    bet ${\beta}$ (long_name='discount factor')
    ph ${\phi}$ (long_name='Preference parameter')
    thet ${\theta}$ (long_name='Labor output share')
    ome_0 ${\omega_0}$ (long_name='Depreciation parameter 0')
    ome_1 ${\omega_1}$ (long_name='Depreciation parameter 1')
    v0 ${\nu_0}$ (long_name='Technical parameter 0')
    v1 ${\nu_1}$ (long_name='Technical parameter 1')
    rho_z ${\rho_z}$ (long_name='autocorrelation parameter for z')
    rho_p ${\rho_p}$ (long_name='autocorrelation parameter for p')
    sig_z ${\sigma_z}$ (long_name='standard deviation for z')
    sig_p ${\sigma_p}$ (long_name='standard deviation for p')
;
         
	bet = .99; 
	alph = 2;
	ph = .3215; 
	thet = .7; 
	ome_0 = .0399;
    ome_1 = 1.2467;
    v0 = .0136;
    v1 = 1.6622;
    rho_z=.95;//taken from Finn 1995 :"a theory of the capacity utilization/inflation relationship"
    rho_p=.95;//same remark
    sig_z=.007;//same remark
    sig_p=.032;//same remark

	/*------------------------------------------------------------
    model: we have 13 variables, we need 13 equations
    ------------------------------------------------------------*/

model;

    [name='production function: $y_{t}=(z_{t}l_{t})^{\theta}(k_{t}u_{t})^{1-\theta}$', eq='\#1']
    y=(exp(z)*l)^thet*(k(-1)*u)^(1-thet);

    [name='Law of motion for capital: $k_{t+1}=(1-\delta(u_{t})k_{t}+i_{t}$', eq='\#2']
    k=(1-delt)*k(-1)+i;

    [name='Depreciation rate of capital, $\delta(u_{t})=\omega_0u_{t}^\omega_1/\omega_1$', eq='\#3']
    delt=(ome_0*u^ome_1)/ome_1;

    [name='Energy relationship, $e_{t}/k_{t}=a(u_{t})$', eq='\#4']
    a=e/k(-1);

    [name='Technical relationship, $a(u_{t})=\nu_0u^{\nu_1}/\nu_1$', eq='\#5']
    a=(v0*u^v1)/v1;

    [name='Household budget constraint, $w_{t}l_{t}+r_{t}k_{t}u_{t}=c_{t}+i_{t}+p_{t}e_{t}$' , eq='\#6']
    w*l+r*k(-1)*u-c-i-exp(p)*e=0;
    
    [name='Labor marginal productivity, $w_{t}=\partial F/z_{t}\partial l_{t} z_{t}$' , eq='\#7']
    w=thet*y/l;

    [name='Capital marginal productivity, $w_{t}=\partial F/\partial k_{t}u_{t}$' , eq='\#8']
    r=(1-thet)*y/(k(-1)*u);

    [name='Intratemporal efficiency condition governing labor supply, $(1-\phi)c_{t}=\phi(1-l_{t})w_{t}$' , eq='\#9']
    (1-ph)*c=ph*(1-l)*w;

    
    [name='Marginal depreciation and energy cost equal marginal return, $d\delta$/d u_{t} k_{t}+p_{t}da_{t}/du_{t}=r_{t}k_{t}$' , eq='\#10']
    ome_0*u^(ome_1-1)*k(-1)+exp(p)*v0*u^(v1-1)*k(-1)=r*k(-1);

    [name='Euler s equation, $\partial U/\partial C_{t}=\beta \partial U/\partial C_{t+1}\Big(r_{t+1}u_{t+1}+1-\delta(u_{t+1})-p_{t+1}a(u_{t+1})\Big)$' , eq='\#11']
    c^(ph*(1-alph)-1)*(1-l)^((1-ph)*(1-alph))=bet*
    c(+1)^(ph*(1-alph)-1)*(1-l(+1))^((1-ph)*(1-alph))*
    (r(+1)*u(+1)+1-delt(+1)-p(+1)*a(+1));

    [name='Law of motion for the technology stochastic process, $z_t=\rho_z z_{t-1}+\epsilon_{t}^z$' , eq='\#12']
    z=rho_z*z(-1)+eps_z;

    [name='Law of motion for the energy price stochastic process, $p_t=\rho_p p_{t-1}+\epsilon_{t}^p$' , eq='\#13']
    p=rho_p*p(-1)+eps_p;
	end;
	
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Computation of the model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

    initval;
    y = 1;
    c = .774;
    i = .183;
    e = .043;
    p = 1;
    z = 1.5462;
    l = .3;
    k = 7.3217;
    u = .82;
    delt = .025;
    a = .0059;
	end;
    resid(1);
	steady;


    shocks;
    var eps_z=sig_z^2;
    var eps_p=sig_p^2;
    end;
    steady;
    stoch_simul(periods=1500,irf=50); 

    write_latex_original_model;
    write_latex_static_model;
    write_latex_dynamic_model(write_equation_tags);
    write_latex_definitions;
    write_latex_steady_state_model;
    