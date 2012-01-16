#include <iostream>
#include <fstream>
#include <ctime>

#include <cfit/parameter.hh>
#include <cfit/coef.hh>
#include <cfit/amplitude.hh>
#include <cfit/phasespace.hh>
#include <cfit/matrix.hh>

#include <cfit/models/decay3body.hh>
#include <cfit/models/relbreitwigner.hh>
#include <cfit/models/gounarissakurai.hh>
#include <cfit/models/glass.hh>

#ifdef MPI_ON
#include <mpi.h>
#endif



const std::vector< double > poles()
{
  std::vector< double > m0;
  m0.push_back( 0.65100 );
  m0.push_back( 1.20360 );
  m0.push_back( 1.55817 );
  m0.push_back( 1.21000 );
  m0.push_back( 1.82206 );

  return m0;
}



const Matrix< double > baseResidueFunctions()
{
  Matrix< double > g0( 5 );

  g0( 0, 0 ) =  0.22889;
  g0( 1, 0 ) =  0.94128;
  g0( 2, 0 ) =  0.36856;
  g0( 3, 0 ) =  0.33650;
  g0( 4, 0 ) =  0.18171;

  g0( 0, 1 ) = -0.55377;
  g0( 1, 1 ) =  0.55095;
  g0( 2, 1 ) =  0.23888;
  g0( 3, 1 ) =  0.40907;
  g0( 4, 1 ) = -0.17558;

  g0( 0, 2 ) =  0.00000;
  g0( 1, 2 ) =  0.00000;
  g0( 2, 2 ) =  0.55639;
  g0( 3, 2 ) =  0.85679;
  g0( 4, 2 ) = -0.79658;

  g0( 0, 3 ) = -0.39899;
  g0( 1, 3 ) =  0.39065;
  g0( 2, 3 ) =  0.18340;
  g0( 3, 3 ) =  0.19906;
  g0( 4, 3 ) = -0.00355;

  g0( 0, 4 ) = -0.34639;
  g0( 1, 4 ) =  0.31503;
  g0( 2, 4 ) =  0.18681;
  g0( 3, 4 ) = -0.00984;
  g0( 4, 4 ) =  0.22358;

  return g0;
}


const Matrix< double > backgroundK()
{
  Matrix< double > fSc( 5 );

  fSc( 0, 0 ) =  0.23399;
  fSc( 0, 1 ) =  0.15044;
  fSc( 0, 2 ) = -0.20545;
  fSc( 0, 3 ) =  0.32825;
  fSc( 0, 4 ) =  0.35412;

  for ( int r = 1; r < 5; ++r )
    fSc( r, 0 ) = fSc( 0, r );

  return fSc;
}


const std::vector< Coef > backgroundP()
{
  Parameter refPr_pipi( "refPr_pipi", -11.424400, 0.111700 );
  Parameter refPr_kk  ( "refPr_kk"  , - 6.602050, 0.367000 );
  Parameter refPr_4pi ( "refPr_4pi" , - 3.795820, 0.731100 );
  Parameter refPr_ee  ( "refPr_ee"  ,   0.000000, 0.000000 );
  Parameter refPr_eep ( "refPr_eep" ,   0.000000, 0.000000 );
  Parameter imfPr_pipi( "imfPr_pipi",   0.056421, 0.104400 );
  Parameter imfPr_kk  ( "imfPr_kk"  ,  14.002300, 0.407000 );
  Parameter imfPr_4pi ( "imfPr_4pi" , - 5.828540, 0.750500 );
  Parameter imfPr_ee  ( "imfPr_ee"  ,   0.000000, 0.000000 );
  Parameter imfPr_eep ( "imfPr_eep" ,   0.000000, 0.000000 );

  Coef fPr_pipi( refPr_pipi, imfPr_pipi );
  Coef fPr_kk  ( refPr_kk  , imfPr_kk   );
  Coef fPr_4pi ( refPr_4pi , imfPr_4pi  );
  Coef fPr_ee  ( refPr_ee  , imfPr_ee   );
  Coef fPr_eep ( refPr_eep , imfPr_eep  );

  std::vector< Coef > fPr;

  fPr.push_back( fPr_pipi );
  fPr.push_back( fPr_kk   );
  fPr.push_back( fPr_4pi  );
  fPr.push_back( fPr_ee   );
  fPr.push_back( fPr_eep  );

  return fPr;
}


const std::vector< Coef > betas()
{
  // Parameters of the P vector pole coefficients.
  Parameter reBeta0( "reBeta0", - 5.532630, 0.059820 );
  Parameter reBeta1( "reBeta1",  15.634400, 0.057220 );
  Parameter reBeta2( "reBeta2",  40.865600, 1.118000 );
  Parameter reBeta3( "reBeta3",   6.139390, 0.249900 );
  Parameter reBeta4( "reBeta4",   0.0     , 0.1      );
  Parameter imBeta0( "imBeta0",   0.298584, 0.040560 );
  Parameter imBeta1( "imBeta1",   0.262657, 0.070300 );
  Parameter imBeta2( "imBeta2", -17.780900, 0.754000 );
  Parameter imBeta3( "imBeta3", - 6.929000, 0.178100 );
  Parameter imBeta4( "imBeta4",   0.0     , 0.1      );

  Coef beta0( reBeta0, imBeta0 );
  Coef beta1( reBeta1, imBeta1 );
  Coef beta2( reBeta2, imBeta2 );
  Coef beta3( reBeta3, imBeta3 );
  Coef beta4( reBeta4, imBeta4 );

  std::vector< Coef > beta;

  beta.push_back( beta0 );
  beta.push_back( beta1 );
  beta.push_back( beta2 );
  beta.push_back( beta3 );
  beta.push_back( beta4 );

  return beta;
}



int main( int argc, char** argv )
{
//   Parameter reCoef1( "re1", 1., .5 );
//   Parameter imCoef1( "im1", 1., .5 );
//   Coef      coef1( reCoef1, imCoef1 );

//   Parameter reCoef2( "re2", 1., .5 );
//   Parameter imCoef2( "im2", 1., .5 );
//   Coef      coef2( reCoef2, imCoef2 );

//   Parameter       mass1 ( "m1", .893, .3   );
//   Parameter       width1( "w1", .047, .003 );
//   Parameter       r1    ( "r1", 1.5 , .5   );
//   RelBreitWigner  reso1 ( 1, 2, mass1, width1, r1, 1 );

//   Parameter       mass2 ( "m2", .893, .3   );
//   Parameter       width2( "w2", .047, .003 );
//   Parameter       r2    ( "r2", 1.5 , .5   );
//   GounarisSakurai reso2 ( 1, 3, mass2, width2, r2, 1 );

  Parameter reCoef_Kstm      ( "reCoef_Kstm"      , -1.196090, 0.005755 );
  Parameter reCoef_Kstp      ( "reCoef_Kstp"      ,  0.118156, 0.003275 );
  Parameter reCoef_rho       ( "reCoef_rho"       ,  1.0     , 0.1      );
  Parameter reCoef_omega     ( "reCoef_omega"     , -0.019235, 0.000670 );
  Parameter reCoef_f2_1270   ( "reCoef_f2_1270"   ,  0.384383, 0.014580 );
  Parameter reCoef_K0stm_1430( "reCoef_K0stm_1430", -0.393116, 0.037100 );
  Parameter reCoef_K0stp_1430( "reCoef_K0stp_1430",  0.059364, 0.029840 );
  Parameter reCoef_K2stm_1430( "reCoef_K2stm_1430",  1.041910, 0.012840 );
  Parameter reCoef_K2stp_1430( "reCoef_K2stp_1430",  0.103058, 0.012930 );
  Parameter reCoef_Kstm_1680 ( "reCoef_Kstm_1680" , -0.888847, 0.033120 );

  Parameter imCoef_Kstm      ( "imCoef_Kstm"      ,  1.256890, 0.006278 );
  Parameter imCoef_Kstp      ( "imCoef_Kstp"      , -0.114066, 0.002994 );
  Parameter imCoef_rho       ( "imCoef_rho"       ,  0.0     , 0.1      );
  Parameter imCoef_omega     ( "imCoef_omega"     ,  0.037376, 0.000529 );
  Parameter imCoef_f2_1270   ( "imCoef_f2_1270"   ,  0.023297, 0.014100 );
  Parameter imCoef_K0stm_1430( "imCoef_K0stm_1430", -5.284940, 0.030800 );
  Parameter imCoef_K0stp_1430( "imCoef_K0stp_1430", -0.284656, 0.028240 );
  Parameter imCoef_K2stm_1430( "imCoef_K2stm_1430", -0.782001, 0.016350 );
  Parameter imCoef_K2stp_1430( "imCoef_K2stp_1430", -0.050667, 0.012820 );
  Parameter imCoef_Kstm_1680 ( "imCoef_Kstm_1680" , -0.150273, 0.040110 );

  Coef coef_Kstm      ( reCoef_Kstm      , imCoef_Kstm       );
  Coef coef_Kstp      ( reCoef_Kstp      , imCoef_Kstp       );
  Coef coef_rho       ( reCoef_rho       , imCoef_rho        );
  Coef coef_omega     ( reCoef_omega     , imCoef_omega      );
  Coef coef_f2_1270   ( reCoef_f2_1270   , imCoef_f2_1270    );
  Coef coef_K0stm_1430( reCoef_K0stm_1430, imCoef_K0stm_1430 );
  Coef coef_K0stp_1430( reCoef_K0stp_1430, imCoef_K0stp_1430 );
  Coef coef_K2stm_1430( reCoef_K2stm_1430, imCoef_K2stm_1430 );
  Coef coef_K2stp_1430( reCoef_K2stp_1430, imCoef_K2stp_1430 );
  Coef coef_Kstm_1680 ( reCoef_Kstm_1680 , imCoef_Kstm_1680  );

  Parameter mKstm      ( "mKstm"      , 0.8937010, 0.1 );
  Parameter mKstp      ( "mKstm"      , 0.8929300, 0.1 );
  Parameter mRho       ( "mRho"       , 0.7758000, 0.1 );
  Parameter momega     ( "momega"     , 0.7825900, 0.1 );
  Parameter mf2_1270   ( "mf2_1270"   , 1.2754000, 0.1 );
  Parameter mK0stm_1430( "mK0stm_1430", 1.4215500, 0.1 );
  Parameter mK0stp_1430( "mK0stm_1430", 1.4470000, 0.1 );
  Parameter mK2stm_1430( "mK2stm_1430", 1.4256000, 0.1 );
  Parameter mK2stp_1430( "mK2stm_1430", 1.4256000, 0.1 );
  Parameter mKstm_1680 ( "mKstm_1680" , 1.6770000, 0.1 );

  Parameter wKstm      ( "wKstm"      , 0.0467401, 0.1 );
  Parameter wKstp      ( "wKstm"      , 0.0466400, 0.1 );
  Parameter wRho       ( "wRho"       , 0.1464000, 0.1 );
  Parameter womega     ( "womega"     , 0.0084900, 0.1 );
  Parameter wf2_1270   ( "wf2_1270"   , 0.1851000, 0.1 );
  Parameter wK0stm_1430( "wK0stm_1430", 0.2467310, 0.1 );
  Parameter wK0stp_1430( "wK0stm_1430", 0.2720000, 0.1 );
  Parameter wK2stm_1430( "wK2stm_1430", 0.0985000, 0.1 );
  Parameter wK2stp_1430( "wK2stm_1430", 0.0985000, 0.1 );
  Parameter wKstm_1680 ( "wKstm_1680" , 0.2050000, 0.1 );

  Parameter rBW        ( "rBW"        , 1.5      , 0.5 ); // Blatt-Weisskopf radius.

  // Generalized Lass parameters.
  Parameter gLassR     ( "lassR"   ,   1.000000, 0.1 );
  Parameter gLassB     ( "lassB"   ,   0.617734, 0.1 );
  Parameter gLassPhiR  ( "lassPhiR",   1.104390, 0.1 );
  Parameter gLassPhiB  ( "lassPhiB", - 0.099519, 0.1 );
  Parameter gLassr     ( "lassr"   , -15.010300, 0.1 );
  Parameter gLassa     ( "lassa"   ,   0.224004, 0.1 );

  RelBreitWigner  prop_Kstm      ( 1, 3, mKstm      , wKstm      , rBW, 1 );
  RelBreitWigner  prop_Kstp      ( 1, 2, mKstp      , wKstp      , rBW, 1 );
  GounarisSakurai prop_rho       ( 2, 3, mRho       , wRho       , rBW, 1 );
  RelBreitWigner  prop_omega     ( 2, 3, momega     , womega     , rBW, 1 );
  RelBreitWigner  prop_f2_1270   ( 2, 3, mf2_1270   , wf2_1270   , rBW, 2 );
  GLass           prop_K0stm_1430( 1, 3, mK0stm_1430, wK0stm_1430, rBW, gLassR, gLassB, gLassPhiR, gLassPhiB, gLassr, gLassa, 0 );
  GLass           prop_K0stp_1430( 1, 2, mK0stp_1430, wK0stp_1430, rBW, gLassR, gLassB, gLassPhiR, gLassPhiB, gLassr, gLassa, 0 );
  RelBreitWigner  prop_K2stm_1430( 1, 3, mK2stm_1430, wK2stm_1430, rBW, 2 );
  RelBreitWigner  prop_K2stp_1430( 1, 2, mK2stp_1430, wK2stp_1430, rBW, 2 );
  RelBreitWigner  prop_Kstm_1680 ( 1, 3, mKstm_1680 , wKstm_1680 , rBW, 1 );

  // Parameters of the F-vector component.
  std::vector< double > m0   = poles();
  Matrix< double >      g0   = baseResidueFunctions();
  Matrix< double >      fSc  = backgroundK();
  std::vector< Coef >   fPr  = backgroundP();
  std::vector< Coef >   beta = betas();

  Parameter             s0pr( "s0pr", -3.92637, 0.1 );

  Fvector fvec( 2, 3, m0, g0, fSc, -3.92637, -0.15, 1.0, beta, fPr, s0pr );

  Amplitude amp;

  amp += coef_Kstm       * prop_Kstm;
  amp += coef_Kstp       * prop_Kstp;
  amp += coef_rho        * prop_rho;
  amp += coef_omega      * prop_omega;
  amp += coef_f2_1270    * prop_f2_1270;

  amp += coef_K0stm_1430 * prop_K0stm_1430;
  amp += coef_K0stp_1430 * prop_K0stp_1430;

  amp += coef_K2stm_1430 * prop_K2stm_1430;
  amp += coef_K2stp_1430 * prop_K2stp_1430;
  amp += coef_Kstm_1680  * prop_Kstm_1680;
  amp += fvec;

  const double& mD0 = 1.8645;
  const double& mKs = 0.49767;
  const double& mPi = 0.139570;

  PhaseSpace ps( mD0, mKs, mPi, mPi );

  Variable mSq12( "mSq12", 1.  );
  Variable mSq13( "mSq13", 1.5 );
  Variable mSq23( "mSq23",  .5 );

  Decay3Body decayModel( mSq12, mSq13, mSq23, amp, ps );

  std::cout << "DEBUG: "                                                        << std::endl;
  std::cout << amp.evaluate( ps, 1.0, 1.5, .5 ) << " " << decayModel.evaluate() << std::endl;
  decayModel.setVar( "mSq23", .6 );
  std::cout << amp.evaluate( ps, 1.0, 1.5, .6 ) << " " << decayModel.evaluate() << std::endl;
  decayModel.setVar( "mSq23", .7 );
  std::cout << amp.evaluate( ps, 1.0, 1.5, .7 ) << " " << decayModel.evaluate() << std::endl;
  decayModel.setVar( "mSq23", .8 );
  std::cout << amp.evaluate( ps, 1.0, 1.5, .8 ) << " " << decayModel.evaluate() << std::endl;

  return 0;
}
