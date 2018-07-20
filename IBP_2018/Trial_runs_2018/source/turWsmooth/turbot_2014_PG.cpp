  #include "admodel.h" 
  #define TRACE(object) tracefile << #object << "\n" << object << "\n\n" << endl;
  ofstream tracefile("data.log");
 
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <turbot_2014_PG.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  nyrs.allocate("nyrs");
  nages.allocate("nages");
  F_time_knots.allocate("F_time_knots");
  catch_weights.allocate(1,nyrs,1,nages,"catch_weights");
  stock_weights.allocate(1,nyrs,1,nages,"stock_weights");
  bs2.allocate(1,F_time_knots,1,nyrs,"bs2");
  catch_weights_NA.allocate(1,nyrs,1,nages);
  stock_weights_NA.allocate(1,nyrs,1,nages);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  logsigmaSWTS.allocate(1,nages,"logsigmaSWTS");
  logsigmaCWTS.allocate(1,nages,"logsigmaCWTS");
  K.allocate(1,"K");
  loga0.allocate(1,"loga0");
  logCWfact.allocate(1,"logCWfact");
  log_temp_wts_Linf.allocate(1,F_time_knots,2,"log_temp_wts_Linf");
  Linf.allocate(1,nyrs,"Linf");
  #ifndef NO_AD_INITIALIZE
    Linf.initialize();
  #endif
  CWT.allocate(1,nyrs,1,nages,"CWT");
  #ifndef NO_AD_INITIALIZE
    CWT.initialize();
  #endif
  SWT.allocate(1,nyrs,1,nages,"SWT");
  #ifndef NO_AD_INITIALIZE
    SWT.initialize();
  #endif
  f_cw.allocate("f_cw");
  #ifndef NO_AD_INITIALIZE
  f_cw.initialize();
  #endif
  f_sw.allocate("f_sw");
  #ifndef NO_AD_INITIALIZE
  f_sw.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
 
  //for all datasources create a copy, but fill with [0,1] for likelihood function
  for (int t=1; t<=nyrs; t++)
    for (int a=1; a<=nages; a++){
      catch_weights_NA(t,a) = (catch_weights(t,a) <0)?0:1;
      catch_weights(t,a)    = (catch_weights(t,a) <0)? -catch_weights(t,a) :catch_weights(t,a) ;    
      stock_weights_NA(t,a) = (stock_weights(t,a) <0)?0:1;
      stock_weights(t,a)    = (stock_weights(t,a) <0)? -stock_weights(t,a) :stock_weights(t,a) ;    
    }
}

void model_parameters::userfunction(void)
{
  f =0.0;
  get_est_wts();
  evaluate_the_objective_function();
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report << "Likelihoods f, f_c, f_s1, f_s2, f_s3, f_cw, f_sw" << endl;
  report << f  << endl << f_cw << endl << f_sw << endl;
  report << "Estimated CWT"     << endl << CWT       << endl;
  report << "Estimated SWT"     << endl << SWT       << endl;
  report << "K"                 << endl << mfexp(-K)  << endl;
  report << "Linf"              << endl << mfexp(Linf) << endl;
  report << "value cwt nyrs " << value(row(CWT,nyrs)) << endl;
}

dvariable model_parameters::dnorm(const dvariable& x, const dvariable& mu, const dvariable& sd)
{
  return 0.5 * (log(2*M_PI*sd*sd) + square(x-mu)/(sd*sd));
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

void model_parameters::evaluate_the_objective_function(void)
{
  f_cw  = 0.0;
  f_sw  = 0.0;
 // Commercial catch at age
  for (int t=1; t<=nyrs; t++){
    for (int a=1; a<=nages; a++){
      f_cw += (catch_weights_NA(t,a)    * dnorm(CWT(t,a),     catch_weights(t,a),         mfexp(logsigmaCWTS(a)) ));  // landing wts
      f_sw += (stock_weights_NA(t,a)    * dnorm(SWT(t,a),     stock_weights(t,a),         mfexp(logsigmaSWTS(a))  ));   // stockwts
    }
  }
  // Add all components
  f = f_cw + f_sw;
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
    if (!arrmblsize) arrmblsize=15000000;
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
