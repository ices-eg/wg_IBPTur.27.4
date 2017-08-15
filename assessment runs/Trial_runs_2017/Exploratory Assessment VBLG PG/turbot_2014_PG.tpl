// ADMB code for simple catch at age model with CAA matrix and three indices
// Turbot assessment

DATA_SECTION
  init_int    nyrs
  init_int    nages
  init_int    qplat
  init_int    minFbar
  init_int    maxFbar
  init_int    F_age_knots
  init_int    F_time_knots
  init_number time_surv1
  init_number time_surv2
  init_number time_surv3
  init_matrix obs_catch_at_age(1,nyrs,1,nages)
  init_matrix catch_weights(1,nyrs,1,nages)
  init_matrix stock_weights(1,nyrs,1,nages)
  init_matrix obs_surv1(1,nyrs,1,nages)
  init_matrix obs_surv2(1,nyrs,1,nages)
  init_matrix obs_surv3(1,nyrs,1,nages)
  init_vector landings(1,nyrs)
  init_number M
  init_vector maturity(1,nages)
  init_matrix bs1(1,F_age_knots,1,qplat)
  init_matrix bs2(1,F_time_knots,1,nyrs)
  vector Fvec(1,11806);
  number recTAC;                       // recruitment at start of tac year
  number ssbTACyr;                     // SSB at start of TAC year
  number intC;                         // intermediate year catches
  vector TAC(1,11806);
  vector SSBstf(1,11806);
  matrix obs_catch_at_age_NA(1,nyrs,1,nages);
  matrix obs_surv1_NA(1,nyrs,1,nages);
  matrix obs_surv2_NA(1,nyrs,1,nages);
  matrix obs_surv3_NA(1,nyrs,1,nages);
  matrix catch_weights_NA(1,nyrs,1,nages);
  matrix stock_weights_NA(1,nyrs,1,nages);
  
PARAMETER_SECTION
  init_number logsigmaC(1)
  init_number logsigmaU1(1)
  init_number logsigmaU2(1)
  init_number logsigmaU3(1)
  init_vector logsigmaSWTS(1,nages)
  init_vector logsigmaCWTS(1,nages) 
  init_number K(1)
  init_number loga0(1)
  init_number logCWfact(1)
  init_vector sigma_offset_CU(2,nages,1)
  init_vector log_sel_coff1(1,F_age_knots,1)
  init_vector log_sel_coff2(1,2,1)
  init_vector log_sel_cofU1(1,F_age_knots,1)
  init_vector log_sel_cofU2(1,F_age_knots,1)
  init_vector log_sel_cofU3(1,F_age_knots,1)
  init_vector log_temp_coff(1,F_time_knots,2)
  init_vector log_temp_wts_Linf(1,F_time_knots,2)
  init_vector log_initpop(1,nyrs+nages-1,1)
  vector sigmaC(1,nages)
  vector sigmaU1(1,nages)
  vector sigmaU2(1,nages)
  vector sigmaU3(1,nages)
  vector effort_devs(1,nyrs)
  vector Linf(1,nyrs)
  vector log_self1(1,nages)
  vector log_self2(1,nages)
  vector log_selU1(1,nages)
  vector log_selU2(1,nages)
  vector log_selU3(1,nages)
  vector TSB(1,nyrs)
  sdreport_vector SSB(1,nyrs)
  vector VB(1,nyrs)
  matrix F(1,nyrs,1,nages)
  matrix S(1,nyrs,1,nages)
  matrix N(1,nyrs,1,nages)
  matrix U1(1,nyrs,1,nages)
  matrix U2(1,nyrs,1,nages)
  matrix U3(1,nyrs,1,nages)
  matrix C(1,nyrs,1,nages)
  matrix CWT(1,nyrs,1,nages)
  matrix SWT(1,nyrs,1,nages)
  number f_c
  number f_s1
  number f_s2
  number f_s3
  number f_cw
  number f_sw
  sdreport_vector Fbar(1,nyrs)
  number Fmax
  number YPRmax
  //number TAC
  objective_function_value f

PRELIMINARY_CALCS_SECTION
  // Create a sequence of Fbar values to evaluate YPR, in order to calculate Fmax:
  // Fvec <- c(0,10^-(9:5),seq(0.0001,1,0.0001),seq(1.001,2,0.001),seq(2.01,10,0.01))
  Fvec(1) = 0;
  Fvec(2) = 1e-9;
  for (int i=3; i<=7; i++)         Fvec(i) = Fvec(i-1) * 10;
  for (int i=8; i<=10006; i++)     Fvec(i) = Fvec(i-1) + 0.0001;
  for (int i=10007; i<=11006; i++) Fvec(i) = Fvec(i-1) + 0.001;
  for (int i=11007; i<=11806; i++) Fvec(i) = Fvec(i-1) + 0.01;

  //for all datasources create a copy, but fill with [0,1] for likelihood function
  for (int t=1; t<=nyrs; t++)
    for (int a=1; a<=nages; a++){
      obs_catch_at_age_NA(t,a) = (obs_catch_at_age(t,a) <0)?0:1;
      obs_catch_at_age(t,a)    = (obs_catch_at_age(t,a) <0)? -obs_catch_at_age(t,a) :obs_catch_at_age(t,a) ;    
      obs_surv1_NA(t,a) = (obs_surv1(t,a) <0)?0:1;
      obs_surv1(t,a)    = (obs_surv1(t,a) <0)? -obs_surv1(t,a) :obs_surv1(t,a) ;    
      obs_surv2_NA(t,a) = (obs_surv2(t,a) <0)?0:1;
      obs_surv2(t,a)    = (obs_surv2(t,a) <0)? -obs_surv2(t,a) :obs_surv2(t,a) ;    
      obs_surv3_NA(t,a) = (obs_surv3(t,a) <0)?0:1;
      obs_surv3(t,a)    = (obs_surv3(t,a) <0)? -obs_surv3(t,a) :obs_surv3(t,a) ;    
      catch_weights_NA(t,a) = (catch_weights(t,a) <0)?0:1;
      catch_weights(t,a)    = (catch_weights(t,a) <0)? -catch_weights(t,a) :catch_weights(t,a) ;    
      stock_weights_NA(t,a) = (stock_weights(t,a) <0)?0:1;
      stock_weights(t,a)    = (stock_weights(t,a) <0)? -stock_weights(t,a) :stock_weights(t,a) ;    
    }

PROCEDURE_SECTION
  get_sigmas();
  get_mortality_and_survival_rates();
  get_numbers_at_age();
  get_catch_at_age();
  get_surveys_at_age();
  get_est_wts();
  calculate_biomass();
  evaluate_the_objective_function();
    get_fmax();
  if (mceval_phase())
  {

    //get_tac();
    write_mcmc();
  }
 

REPORT_SECTION
  report << "Likelihoods f, f_c, f_s1, f_s2, f_s3, f_cw, f_sw" << endl;
  report << f << endl << f_c << endl << f_s1 << endl << f_s2 << endl << f_s3 << endl << f_cw << endl << f_sw << endl;
  report << "log_self1"         << endl << log_self1 << endl;
  report << "log_self2"         << endl << log_self2 << endl;
  report << "log_selU1"         << endl << log_selU1 << endl;
  report << "log_selU2"         << endl << log_selU2 << endl;
  report << "log_selU3"         << endl << log_selU3 << endl;
  report << "sigmaC"            << endl << sigmaC    << endl;
  report << "sigmaU1"           << endl << sigmaU1   << endl;
  report << "sigmaU2"           << endl << sigmaU2   << endl;
  report << "sigmaU3"           << endl << sigmaU3   << endl;
  report << "Estimated c@a"     << endl << C         << endl;
  report << "Estimated survey1" << endl << U1        << endl;
  report << "Estimated survey2" << endl << U2        << endl;
  report << "Estimated survey3" << endl << U3        << endl;
  report << "Estimated N"       << endl << N         << endl;
  report << "Estimated F"       << endl << F         << endl;
  report << "Estimated Fbar (" << minFbar << "-" << maxFbar << ")" << endl << Fbar << endl ;
  report << "Estimated SSB"     << endl << SSB       << endl;
  report << "Estimated TSB"     << endl << TSB       << endl;
  report << "Estimated CWT"     << endl << CWT       << endl;
  report << "Estimated SWT"     << endl << SWT       << endl;
  report << "K"                 << endl << mfexp(-K)  << endl;
  report << "Linf"              << endl << mfexp(Linf) << endl;
  report << "STF" << endl;
  report << "Int_yr_rec Int_yr_ssb Int_yr_landings" << endl;
  report << recTAC << " " << ssbTACyr << " " << intC << endl;
  report << "Fvec TAC resultant_ssb" << endl;
  report << Fvec << endl;
  report << TAC << endl;
  report << SSBstf << endl;
  report << "value cwt nyrs " << value(row(CWT,nyrs)) << endl;
  report <<  Fmax <<  " " << YPRmax << endl;

FUNCTION dvariable dnorm(const dvariable& x, const dvariable& mu, const dvariable& sd)
  return 0.5 * (log(2*M_PI*sd*sd) + square(x-mu)/(sd*sd));

FUNCTION get_sigmas
  // Fleet sigma 
  sigmaC(1) = mfexp(logsigmaC);
  // Survey sigma
  sigmaU1(1) = mfexp(logsigmaU1);
  sigmaU2(1) = mfexp(logsigmaU2);
  sigmaU3(1) = mfexp(logsigmaU3);
  // swt sigma
  for (int a=2; a<=nages; a++){
    sigmaC(a)  = mfexp(logsigmaC + sigma_offset_CU(a));
    sigmaU1(a) = mfexp(logsigmaU1 + sigma_offset_CU(a)); 
    sigmaU2(a) = mfexp(logsigmaU2 + sigma_offset_CU(a));
    sigmaU3(a) = mfexp(logsigmaU3 + sigma_offset_CU(a));
  }

FUNCTION get_est_wts
   Linf = mfexp(log_temp_wts_Linf) * bs2;
  for (int t=1; t<=nyrs; t++){
    for (int a=1; a<=nages; a++){
      CWT(t,a) = 0.00001508 *  pow(mfexp(Linf(t))*(1-exp(-mfexp(-K)*(a + mfexp(loga0)))), 3.090 ) ; //1986 cefas report on literature   
      SWT(t,a) = CWT(t,a) * mfexp(logCWfact);
    }
  }

FUNCTION get_mortality_and_survival_rates
  // Calculate selectivity from sel_coffs, where selectivity is same for last two ages
  log_self1(1,qplat) = elem_div(exp(log_sel_coff1*bs1), 1+exp(log_sel_coff1*bs1));
  for (int a=qplat+1; a<=nages; a++)
    log_self1(a) =  log_self1(qplat) ;
  //calc sels in final period
  for (int a=1; a<=2; a++)
    log_self2(a) =  log_self1(a)+  exp(log_sel_coff2(a)) ;
  for (int a=3; a<=nages; a++)
    log_self2(a) = log_self1(a);
  effort_devs = exp(log_temp_coff * bs2);

  // F = outer_prod(effort_devs,log_self);
  for (int t=1; t<=27; t++)
    for (int a=1; a<=nages; a++)
      F(t,a) = effort_devs(t) * log_self1(a);
  for (int t=28; t<=nyrs; t++)
    for (int a=1; a<=nages; a++)
      F(t,a) = effort_devs(t) * log_self2(a);
  for (int t=1; t<=nyrs; t++)
    Fbar(t) = mean(row(F,t)(minFbar,maxFbar));
  S = mfexp(-(F+M));

FUNCTION get_numbers_at_age
  for (int t=1; t<=nyrs; t++)
    N(t,1) = mfexp(log_initpop(t));
  for (int a=2; a<=nages; a++)
    N(1,a) = mfexp(log_initpop(nyrs+a-1));
  for (int t=1; t<nyrs; t++)
    for (int a=1; a<(nages-1); a++)
      N(t+1,a+1) = N(t,a) * S(t,a);
// plusgroup
  for (int t=1; t<nyrs; t++)
    N(t+1,nages) = N(t,nages) * S(t,nages) + N(t,nages-1) * S(t,nages-1);

FUNCTION get_catch_at_age
  C = elem_prod(elem_div(F,(F+M)), elem_prod(1-S,N));

FUNCTION get_surveys_at_age
  log_selU1(1,qplat) =  elem_div(exp(log_sel_cofU1*bs1), 1+exp(log_sel_cofU1*bs1));
  for (int a=qplat+1; a<=nages; a++)
    log_selU1(a) =  log_selU1(qplat) ;
  log_selU2(1,qplat) =  elem_div(exp(log_sel_cofU2*bs1), 1+exp(log_sel_cofU2*bs1));
  for (int a=qplat+1; a<=nages; a++)
    log_selU2(a) =  log_selU2(qplat) ;
  log_selU3(1,qplat) = elem_div(exp(log_sel_cofU3*bs1), 1+exp(log_sel_cofU3*bs1));
  for (int a=qplat+1; a<=nages; a++)
    log_selU3(a) =  log_selU3(qplat) ;

  for (int t=1; t<=nyrs; t++){
    for (int a=1; a<=nages; a++){
      U1(t,a) = log_selU1(a) * N(t,a) * mfexp(-time_surv1*(F(t,a)+M));
      U2(t,a) = log_selU2(a) * N(t,a) * mfexp(-time_surv2*(F(t,a)+M));
      U3(t,a) = log_selU3(a) * N(t,a) * mfexp(-time_surv3*(F(t,a)+M));
    }
  }

FUNCTION  calculate_biomass
  SSB = maturity * trans(elem_prod(N, SWT));
  TSB = rowsum(elem_prod(N, SWT));
  for (int t=1; t<=28; t++)
    VB(t) = log_self1 * elem_prod(N(t), CWT(t));  // biomass vulnerable to fleet1
  for (int t=29; t<=nyrs; t++)
    VB(t) = log_self2 * elem_prod(N(t), CWT(t));  // biomass vulnerable to fleet2

FUNCTION evaluate_the_objective_function
  f_c  = 0.0;
  f_s1 = 0.0;
  f_s2 = 0.0;
  f_s3 = 0.0;
  f_cw  = 0.0;
  f_sw  = 0.0;
 // Commercial catch at age
  for (int t=1; t<=nyrs; t++){
    for (int a=1; a<=nages; a++){
      f_c  += (obs_catch_at_age_NA(t,a) * dnorm(log(C(t,a)),  log(obs_catch_at_age(t,a)), sigmaC(a)));
      f_s1 += (obs_surv1_NA(t,a)        * dnorm(log(U1(t,a)), log(obs_surv1(t,a)),        sigmaU1(a))); // Survey 1
      f_s2 += (obs_surv2_NA(t,a)        * dnorm(log(U2(t,a)), log(obs_surv2(t,a)),        sigmaU2(a)));  // Survey 2
      f_s3 += (obs_surv3_NA(t,a)        * dnorm(log(U3(t,a)), log(obs_surv3(t,a)),        sigmaU3(a)));  // Survey 3
      f_cw += (catch_weights_NA(t,a)    * dnorm(CWT(t,a),     catch_weights(t,a),         mfexp(logsigmaCWTS(a)) ));  // landing wts
      f_sw += (stock_weights_NA(t,a)    * dnorm(SWT(t,a),     stock_weights(t,a),         mfexp(logsigmaSWTS(a))  ));   // stockwts
    }
  }
  // Add all components
  f = f_c + f_s1 + f_s2 + f_s3 + f_cw + f_sw;


FUNCTION get_fmax
  int i = 0;                            // element in Fvec being evaluated
  bool found = false;                   // whether Fmax is found
  dvector sel = value(log_self2);       // selectivity to use, not necessarily between 0 and 1
  dvector f(1,nages);                   // F at age when Fvec(i) is applied
  dvector z(1,nages);                   // Z = F+M
  dvector n(1,nages);                   // equilibrium population
  dvector c(1,nages);                   // equilibrium catches
  dvector w = value(row(CWT,nyrs));     // catch weights to use
  double ypr = -1;                      // highest YPR found so far
  double proposal;                      // YPR being evaluated
  n(1) = 1;
  while (!found)
  {
    i++;
    f = Fvec(i) * sel;
    z = f + M;
    for (int a=2; a<=nages; a++)
      n(a) = n(a-1) * exp(-(f(a-1)+M));
    for (int a=1; a<=nages; a++)
      c(a) = f(a)/z(a) * n(a) * (1-exp(-z(a)));
    proposal = sum(elem_prod(c, w));
    if (proposal > ypr)
    {
      ypr = proposal;
    }
    else
    {
      i--;  // move i back to optimum
      found = true;
    }
  }
  Fmax = mean(f(minFbar,maxFbar));
  YPRmax = ypr;

FUNCTION write_mcmc
  // Fbar
  mcmc_f << Fbar << endl;
  // Recruitment
  mcmc_rec << column(N,1) << endl;
  // Biomass
  mcmc_ssb << SSB << endl;
  // Reference points
  mcmc_ref << Fmax << endl;
  mcmc_ypr << YPRmax << endl;
  for (int a=1; a<=nages; a++){
    mcmc_cwtfinalyr << CWT(nyrs,a) << " ";
  }
   mcmc_cwtfinalyr << endl;
  for (int a=1; a<=nages; a++){
    mcmc_swtfinalyr << SWT(nyrs,a) << " " ;
  } 
  mcmc_swtfinalyr << endl;
  for (int a=1; a<=nages; a++){
    mcmc_self       <<  log_self2(a) << " " ;
  }
  mcmc_self       <<   endl;

RUNTIME_SECTION
  maximum_function_evaluations 2000, 2000, 2000

GLOBALS_SECTION
  #include "admodel.h" 
  #define TRACE(object) tracefile << #object << "\n" << object << "\n\n" << endl;
  ofstream tracefile("data.log");
  ofstream mcmc_f("f.mcmc");
  ofstream mcmc_rec("rec.mcmc");
  ofstream mcmc_ssb("ssb.mcmc");
  ofstream mcmc_ref("ref.mcmc");
  ofstream mcmc_ypr("ypr.mcmc");
  ofstream mcmc_cwtfinalyr("cwtfinalyr.mcmc"); 
  ofstream mcmc_swtfinalyr("swtfinalyr.mcmc"); 
  ofstream mcmc_self("self.mcmc");
