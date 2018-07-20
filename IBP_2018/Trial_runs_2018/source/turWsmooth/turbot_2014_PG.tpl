// ADMB code for simple catch at age model with CAA matrix and three indices
// Turbot assessment

DATA_SECTION
  init_int    nyrs
  init_int    nages
  init_int    F_time_knots
  init_matrix catch_weights(1,nyrs,1,nages)
  init_matrix stock_weights(1,nyrs,1,nages)
  init_matrix bs2(1,F_time_knots,1,nyrs)
  matrix catch_weights_NA(1,nyrs,1,nages);
  matrix stock_weights_NA(1,nyrs,1,nages);
  
PARAMETER_SECTION
  init_vector logsigmaSWTS(1,nages)
  init_vector logsigmaCWTS(1,nages) 
  init_number K(1)
  init_number loga0(1)
  init_number logCWfact(1)
  init_vector log_temp_wts_Linf(1,F_time_knots,2)
  vector Linf(1,nyrs)
  matrix CWT(1,nyrs,1,nages)
  matrix SWT(1,nyrs,1,nages)
  number f_cw
  number f_sw
  //number TAC
  objective_function_value f

PRELIMINARY_CALCS_SECTION
 
  //for all datasources create a copy, but fill with [0,1] for likelihood function
  for (int t=1; t<=nyrs; t++)
    for (int a=1; a<=nages; a++){
      catch_weights_NA(t,a) = (catch_weights(t,a) <0)?0:1;
      catch_weights(t,a)    = (catch_weights(t,a) <0)? -catch_weights(t,a) :catch_weights(t,a) ;    
      stock_weights_NA(t,a) = (stock_weights(t,a) <0)?0:1;
      stock_weights(t,a)    = (stock_weights(t,a) <0)? -stock_weights(t,a) :stock_weights(t,a) ;    
    }

PROCEDURE_SECTION
  get_est_wts();
  evaluate_the_objective_function();
 
REPORT_SECTION
  report << "Likelihoods f, f_c, f_s1, f_s2, f_s3, f_cw, f_sw" << endl;
  report << f  << endl << f_cw << endl << f_sw << endl;
  report << "Estimated CWT"     << endl << CWT       << endl;
  report << "Estimated SWT"     << endl << SWT       << endl;
  report << "K"                 << endl << mfexp(-K)  << endl;
  report << "Linf"              << endl << mfexp(Linf) << endl;
  report << "value cwt nyrs " << value(row(CWT,nyrs)) << endl;

FUNCTION dvariable dnorm(const dvariable& x, const dvariable& mu, const dvariable& sd)
  return 0.5 * (log(2*M_PI*sd*sd) + square(x-mu)/(sd*sd));

FUNCTION get_est_wts
   Linf = mfexp(log_temp_wts_Linf) * bs2;
  for (int t=1; t<=nyrs; t++){
    for (int a=1; a<=nages; a++){
      CWT(t,a) = 0.00001508 *  pow(mfexp(Linf(t))*(1-exp(-mfexp(-K)*(a + mfexp(loga0)))), 3.090 ) ; //1986 cefas report on literature   
      SWT(t,a) = CWT(t,a) * mfexp(logCWfact);
    }
  }

FUNCTION evaluate_the_objective_function
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

RUNTIME_SECTION
  maximum_function_evaluations 2000, 2000, 2000

GLOBALS_SECTION
  #include "admodel.h" 
  #define TRACE(object) tracefile << #object << "\n" << object << "\n\n" << endl;
  ofstream tracefile("data.log");
 
