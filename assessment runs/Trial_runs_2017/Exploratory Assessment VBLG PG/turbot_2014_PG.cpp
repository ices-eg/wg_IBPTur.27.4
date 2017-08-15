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
#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <turbot_2014_PG.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  nyrs.allocate("nyrs");
  nages.allocate("nages");
  qplat.allocate("qplat");
  minFbar.allocate("minFbar");
  maxFbar.allocate("maxFbar");
  F_age_knots.allocate("F_age_knots");
  F_time_knots.allocate("F_time_knots");
  time_surv1.allocate("time_surv1");
  time_surv2.allocate("time_surv2");
  time_surv3.allocate("time_surv3");
  obs_catch_at_age.allocate(1,nyrs,1,nages,"obs_catch_at_age");
  catch_weights.allocate(1,nyrs,1,nages,"catch_weights");
  stock_weights.allocate(1,nyrs,1,nages,"stock_weights");
  obs_surv1.allocate(1,nyrs,1,nages,"obs_surv1");
  obs_surv2.allocate(1,nyrs,1,nages,"obs_surv2");
  obs_surv3.allocate(1,nyrs,1,nages,"obs_surv3");
  landings.allocate(1,nyrs,"landings");
  M.allocate("M");
  maturity.allocate(1,nages,"maturity");
  bs1.allocate(1,F_age_knots,1,qplat,"bs1");
  bs2.allocate(1,F_time_knots,1,nyrs,"bs2");
  Fvec.allocate(1,11806);
  TAC.allocate(1,11806);
  SSBstf.allocate(1,11806);
  obs_catch_at_age_NA.allocate(1,nyrs,1,nages);
  obs_surv1_NA.allocate(1,nyrs,1,nages);
  obs_surv2_NA.allocate(1,nyrs,1,nages);
  obs_surv3_NA.allocate(1,nyrs,1,nages);
  catch_weights_NA.allocate(1,nyrs,1,nages);
  stock_weights_NA.allocate(1,nyrs,1,nages);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  logsigmaC.allocate(1,"logsigmaC");
  logsigmaU1.allocate(1,"logsigmaU1");
  logsigmaU2.allocate(1,"logsigmaU2");
  logsigmaU3.allocate(1,"logsigmaU3");
  logsigmaSWTS.allocate(1,nages,"logsigmaSWTS");
  logsigmaCWTS.allocate(1,nages,"logsigmaCWTS");
  K.allocate(1,"K");
  loga0.allocate(1,"loga0");
  logCWfact.allocate(1,"logCWfact");
  sigma_offset_CU.allocate(2,nages,1,"sigma_offset_CU");
  log_sel_coff1.allocate(1,F_age_knots,1,"log_sel_coff1");
  log_sel_coff2.allocate(1,2,1,"log_sel_coff2");
  log_sel_cofU1.allocate(1,F_age_knots,1,"log_sel_cofU1");
  log_sel_cofU2.allocate(1,F_age_knots,1,"log_sel_cofU2");
  log_sel_cofU3.allocate(1,F_age_knots,1,"log_sel_cofU3");
  log_temp_coff.allocate(1,F_time_knots,2,"log_temp_coff");
  log_temp_wts_Linf.allocate(1,F_time_knots,2,"log_temp_wts_Linf");
  log_initpop.allocate(1,nyrs+nages-1,1,"log_initpop");
  sigmaC.allocate(1,nages,"sigmaC");
  #ifndef NO_AD_INITIALIZE
    sigmaC.initialize();
  #endif
  sigmaU1.allocate(1,nages,"sigmaU1");
  #ifndef NO_AD_INITIALIZE
    sigmaU1.initialize();
  #endif
  sigmaU2.allocate(1,nages,"sigmaU2");
  #ifndef NO_AD_INITIALIZE
    sigmaU2.initialize();
  #endif
  sigmaU3.allocate(1,nages,"sigmaU3");
  #ifndef NO_AD_INITIALIZE
    sigmaU3.initialize();
  #endif
  effort_devs.allocate(1,nyrs,"effort_devs");
  #ifndef NO_AD_INITIALIZE
    effort_devs.initialize();
  #endif
  Linf.allocate(1,nyrs,"Linf");
  #ifndef NO_AD_INITIALIZE
    Linf.initialize();
  #endif
  log_self1.allocate(1,nages,"log_self1");
  #ifndef NO_AD_INITIALIZE
    log_self1.initialize();
  #endif
  log_self2.allocate(1,nages,"log_self2");
  #ifndef NO_AD_INITIALIZE
    log_self2.initialize();
  #endif
  log_selU1.allocate(1,nages,"log_selU1");
  #ifndef NO_AD_INITIALIZE
    log_selU1.initialize();
  #endif
  log_selU2.allocate(1,nages,"log_selU2");
  #ifndef NO_AD_INITIALIZE
    log_selU2.initialize();
  #endif
  log_selU3.allocate(1,nages,"log_selU3");
  #ifndef NO_AD_INITIALIZE
    log_selU3.initialize();
  #endif
  TSB.allocate(1,nyrs,"TSB");
  #ifndef NO_AD_INITIALIZE
    TSB.initialize();
  #endif
  SSB.allocate(1,nyrs,"SSB");
  VB.allocate(1,nyrs,"VB");
  #ifndef NO_AD_INITIALIZE
    VB.initialize();
  #endif
  F.allocate(1,nyrs,1,nages,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  S.allocate(1,nyrs,1,nages,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  N.allocate(1,nyrs,1,nages,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  U1.allocate(1,nyrs,1,nages,"U1");
  #ifndef NO_AD_INITIALIZE
    U1.initialize();
  #endif
  U2.allocate(1,nyrs,1,nages,"U2");
  #ifndef NO_AD_INITIALIZE
    U2.initialize();
  #endif
  U3.allocate(1,nyrs,1,nages,"U3");
  #ifndef NO_AD_INITIALIZE
    U3.initialize();
  #endif
  C.allocate(1,nyrs,1,nages,"C");
  #ifndef NO_AD_INITIALIZE
    C.initialize();
  #endif
  CWT.allocate(1,nyrs,1,nages,"CWT");
  #ifndef NO_AD_INITIALIZE
    CWT.initialize();
  #endif
  SWT.allocate(1,nyrs,1,nages,"SWT");
  #ifndef NO_AD_INITIALIZE
    SWT.initialize();
  #endif
  f_c.allocate("f_c");
  #ifndef NO_AD_INITIALIZE
  f_c.initialize();
  #endif
  f_s1.allocate("f_s1");
  #ifndef NO_AD_INITIALIZE
  f_s1.initialize();
  #endif
  f_s2.allocate("f_s2");
  #ifndef NO_AD_INITIALIZE
  f_s2.initialize();
  #endif
  f_s3.allocate("f_s3");
  #ifndef NO_AD_INITIALIZE
  f_s3.initialize();
  #endif
  f_cw.allocate("f_cw");
  #ifndef NO_AD_INITIALIZE
  f_cw.initialize();
  #endif
  f_sw.allocate("f_sw");
  #ifndef NO_AD_INITIALIZE
  f_sw.initialize();
  #endif
  Fbar.allocate(1,nyrs,"Fbar");
  Fmax.allocate("Fmax");
  #ifndef NO_AD_INITIALIZE
  Fmax.initialize();
  #endif
  YPRmax.allocate("YPRmax");
  #ifndef NO_AD_INITIALIZE
  YPRmax.initialize();
  #endif
  f.allocate("f");
}

void model_parameters::preliminary_calculations(void)
{

  admaster_slave_variable_interface(*this);
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
}

void model_parameters::userfunction(void)
{
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
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
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
}

dvariable model_parameters::dnorm(const dvariable& x, const dvariable& mu, const dvariable& sd)
{
  return 0.5 * (log(2*M_PI*sd*sd) + square(x-mu)/(sd*sd));
}

void model_parameters::get_sigmas(void)
{
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
}

void model_parameters::get_est_wts(void)
{
   Linf = mfexp(log_temp_wts_Linf) * bs2;
  for (int t=1; t<=nyrs; t++){
    for (int a=1; a<=nages; a++){
      CWT(t,a) = 0.00001508 *  pow(mfexp(Linf(t))*(1-exp(-mfexp(-K)*(a + mfexp(loga0)))), 3.090 ) ; //1986 cefas report on literature   
      SWT(t,a) = CWT(t,a) * mfexp(logCWfact);
    }
  }
}

void model_parameters::get_mortality_and_survival_rates(void)
{
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
}

void model_parameters::get_numbers_at_age(void)
{
  for (int t=1; t<=nyrs; t++)
    N(t,1) = mfexp(log_initpop(t));
  for (int a=2; a<=nages; a++)
    N(1,a) = mfexp(log_initpop(nyrs+a-1));
  for (int t=1; t<nyrs; t++)
    for (int a=1; a<(nages-1); a++)
      N(t+1,a+1) = N(t,a) * S(t,a);
  for (int t=1; t<nyrs; t++)
    N(t+1,nages) = N(t,nages) * S(t,nages) + N(t,nages-1) * S(t,nages-1);
}

void model_parameters::get_catch_at_age(void)
{
  C = elem_prod(elem_div(F,(F+M)), elem_prod(1-S,N));
}

void model_parameters::get_surveys_at_age(void)
{
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
}

void model_parameters::calculate_biomass(void)
{
  SSB = maturity * trans(elem_prod(N, SWT));
  TSB = rowsum(elem_prod(N, SWT));
  for (int t=1; t<=28; t++)
    VB(t) = log_self1 * elem_prod(N(t), CWT(t));  // biomass vulnerable to fleet1
  for (int t=29; t<=nyrs; t++)
    VB(t) = log_self2 * elem_prod(N(t), CWT(t));  // biomass vulnerable to fleet2
}

void model_parameters::evaluate_the_objective_function(void)
{
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
}

void model_parameters::get_fmax(void)
{
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
}

void model_parameters::write_mcmc(void)
{
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
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{2000, 2000, 2000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
  #if defined(__GNUDOS__) || defined(DOS386) || defined(__DPMI32__)  || \
     defined(__MSVC32__)
      if (!arrmblsize) arrmblsize=150000;
  #else
      if (!arrmblsize) arrmblsize=25000;
  #endif
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
