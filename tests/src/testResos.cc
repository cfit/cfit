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

#ifdef MPI_ON
#include <mpi.h>
#endif



std::vector< double > poles()
{
  std::vector< double > m0;
  m0.push_back( 0.65100 );
  m0.push_back( 1.20360 );
  m0.push_back( 1.55817 );
  m0.push_back( 1.21000 );
  m0.push_back( 1.82206 );

  return m0;
}



Matrix< double > baseResidueFunctions()
{
  Matrix< double > g0( 5 );

  g0( 0, 0 ) =  0.22889;
  g0( 0, 1 ) = -0.55377;
  g0( 0, 2 ) =  0.00000;
  g0( 0, 3 ) = -0.39899;
  g0( 0, 4 ) = -0.34639;

  g0( 1, 0 ) =  0.94128;
  g0( 1, 1 ) =  0.55095;
  g0( 1, 2 ) =  0.00000;
  g0( 1, 3 ) =  0.39065;
  g0( 1, 4 ) =  0.31503;

  g0( 2, 0 ) =  0.36856;
  g0( 2, 1 ) =  0.23888;
  g0( 2, 2 ) =  0.55639;
  g0( 2, 3 ) =  0.18340;
  g0( 2, 4 ) =  0.18681;

  g0( 3, 0 ) =  0.33650;
  g0( 3, 1 ) =  0.40907;
  g0( 3, 2 ) =  0.85679;
  g0( 3, 3 ) =  0.19906;
  g0( 3, 4 ) = -0.00984;

  g0( 4, 0 ) =  0.18171;
  g0( 4, 1 ) = -0.17558;
  g0( 4, 2 ) = -0.79658;
  g0( 4, 3 ) = -0.00355;
  g0( 4, 4 ) =  0.22358;

  return g0;
}


Matrix< double > background()
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


int main( int argc, char** argv )
{
  Parameter reCoef1( "re1", 1., .5 );
  Parameter imCoef1( "im1", 1., .5 );
  Coef      coef1( reCoef1, imCoef1 );

  Parameter reCoef2( "re2", 1., .5 );
  Parameter imCoef2( "im2", 1., .5 );
  Coef      coef2( reCoef2, imCoef2 );

  Parameter       mass1 ( "m1", .893, .3   );
  Parameter       width1( "w1", .047, .003 );
  Parameter       r1    ( "r1", 1.5 , .5   );
  RelBreitWigner  reso1 ( 1, 2, mass1, width1, r1, 1 );

  Parameter       mass2 ( "m2", .893, .3   );
  Parameter       width2( "w2", .047, .003 );
  Parameter       r2    ( "r2", 1.5 , .5   );
  GounarisSakurai reso2 ( 1, 3, mass2, width2, r2, 1 );


  std::vector< double > m0  = poles();
  Matrix< double >      g0  = baseResidueFunctions();
  Matrix< double >      fSc = background();

//   std::cout << g0 .dump() << std::endl;
//   std::cout << fSc.dump() << std::endl;

  std::vector< Coef >      beta;
  std::vector< Parameter > fPr;
  Parameter                s0pr;

  beta.push_back( coef1 );
  beta.push_back( coef1 );
  beta.push_back( coef1 );
  beta.push_back( coef2 );
  beta.push_back( coef2 );

  fPr.push_back( mass1 );
  fPr.push_back( mass1 );
  fPr.push_back( mass1 );
  fPr.push_back( mass2 );
  fPr.push_back( mass2 );

  s0pr = width1;

  Fvector fvec( 2, 3, m0, g0, fSc, -3.92637, -0.15, 1.0, beta, fPr, s0pr );

  Amplitude amp;
//   amp  = coef1 * reso1 + coef2 * reso2 + fvec;
  amp = fvec;

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
