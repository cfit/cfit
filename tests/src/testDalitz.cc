#include <iostream>
#include <fstream>
#include <ctime>

#include <cfit/parameter.hh>
#include <cfit/variable.hh>
#include <cfit/dataset.hh>
#include <cfit/chi2.hh>
#include <cfit/nll.hh>
#include <cfit/pdf.hh>

#include <cfit/models/d0mix.hh>

#include <Minuit/MnPrint.h>

#ifdef MPI_ON
#include <mpi.h>
#endif


void readData ( std::string fileName, Dataset& data )
{
  double m2AB;
  double m2AC;
  double m2BC;
  double time;

  std::ifstream file( fileName.c_str() );
  while ( file >> m2AB >> m2AC >> m2BC >> time )
    {
      data.push( "m2AB", m2AB );
      data.push( "m2AC", m2AC );
      data.push( "m2BC", m2BC );
      data.push( "t"   , time );
    }
  file.close();

  return;
}



int main( int argc, char** argv )
{
  // If working with MPI, initialize it.
#ifdef MPI_ON
  MPI::Init( argc, argv );
  const int rank = MPI::COMM_WORLD.Get_rank();
#endif

  clock_t initClock = clock();

  // Data container.
  Dataset data;

  // Read the data. If working with MPI, only the root
  //    process reads the data, then it gets scattered.
#ifdef MPI_ON
  if ( rank == 0 )
    readData( "data/dalitz.dat", data );
  data.scatter();
#else
  readData( "data/dalitz.dat", data );
#endif

  // Variables the model depends on.
  Variable m2AB( "m2AB" );
  Variable m2AC( "m2AC" );
  Variable m2BC( "m2BC" );
  Variable time( "t"    );

  // Parameters of the model.
  Parameter tau           ( "tau"           ,   .41     , .1 );
  Parameter x             ( "x"             ,   .0      , .1 );
  Parameter y             ( "y"             ,   .0      , .1 );

  // Masses and widths.
  Parameter mKStarc       ( "mKStarc"       , 0.893445  , .1 );
  Parameter wKStarc       ( "wKStarc"       , 0.0466057 , .1 );
  Parameter mRho          ( "mRho"          , 0.7758    , .1 );
  Parameter wRho          ( "wRho"          , 0.1464    , .1 );
  Parameter mOmega        ( "mOmega"        , 0.78259   , .1 );
  Parameter wOmega        ( "wOmega"        , 0.00849   , .1 );
  Parameter mF0_980       ( "mF0_980"       , 0.975     , .1 );
  Parameter wF0_980       ( "wF0_980"       , 0.044     , .1 );
  Parameter mF0_1370      ( "mF0_1370"      , 1.434     , .1 );
  Parameter wF0_1370      ( "wF0_1370"      , 0.173     , .1 );
  Parameter mF2_1270      ( "mF2_1270"      , 1.2754    , .1 );
  Parameter wF2_1270      ( "wF2_1270"      , 0.1851    , .1 );
  Parameter mK0Starc_1430 ( "mK0Starc_1430" , 1.459     , .1 );
  Parameter wK0Starc_1430 ( "wK0Starc_1430" , 0.175     , .1 );
  Parameter mK2Starc_1430 ( "mK2Starc_1430" , 1.4256    , .1 );
  Parameter wK2Starc_1430 ( "wK2Starc_1430" , 0.0985    , .1 );
  Parameter mSigma        ( "mSigma"        , 0.518288  , .1 );
  Parameter wSigma        ( "wSigma"        , 0.525990  , .1 );
  Parameter mSigma2       ( "mSigma2"       , 1.02784   , .1 );
  Parameter wSigma2       ( "wSigma2"       , 0.0952316 , .1 );
  Parameter mKStarm_1680  ( "mKStarm_1680"  , 1.677     , .1 );
  Parameter wKStarm_1680  ( "wKStarm_1680"  , 0.205     , .1 );

  // Real and imaginary parts of the isobar coefficients.
  Parameter rKstarm       ( "rKstarm"      , -1.16721   , .1 );
  Parameter iKstarm       ( "iKstarm"      ,  1.20719   , .1 );
  Parameter rKstarp       ( "rKstarp"      ,  0.108452  , .1 );
  Parameter iKstarp       ( "iKstarp"      , -0.115997  , .1 );
  Parameter rRho          ( "rRho"         ,  1.0       , .1 );
  Parameter iRho          ( "iRho"         ,  0.0       , .1 );
  Parameter rOmega        ( "rOmega"       , -0.0246859 , .1 );
  Parameter iOmega        ( "iOmega"       ,  0.0391864 , .1 );
  Parameter rF0_980       ( "rF0_980"      , -0.427196  , .1 );
  Parameter iF0_980       ( "iF0_980"      , -0.244929  , .1 );
  Parameter rF0_1370      ( "rF0_1370"     , -2.6895    , .1 );
  Parameter iF0_1370      ( "iF0_1370"     ,  3.82566   , .1 );
  Parameter rF2_1270      ( "rF2_1270"     ,  0.237268  , .1 );
  Parameter iF2_1270      ( "iF2_1270"     , -0.154790  , .1 );
  Parameter rK0starm_1430 ( "rK0starm_1430",  1.62243   , .1 );
  Parameter iK0starm_1430 ( "iK0starm_1430",  1.09156   , .1 );
  Parameter rK0starp_1430 ( "rK0starp_1430",  0.166794  , .1 );
  Parameter iK0starp_1430 ( "iK0starp_1430",  0.0854312 , .1 );
  Parameter rK2starm_1430 ( "rK2starm_1430",  1.16792   , .1 );
  Parameter iK2starm_1430 ( "iK2starm_1430", -0.760669  , .1 );
  Parameter rK2starp_1430 ( "rK2starp_1430",  0.128534  , .1 );
  Parameter iK2starp_1430 ( "iK2starp_1430", -0.169843  , .1 );
  Parameter rSigma        ( "rSigma"       , -1.64117   , .1 );
  Parameter iSigma        ( "iSigma"       , -0.833011  , .1 );
  Parameter rSigma2       ( "rSigma2"      , -0.275866  , .1 );
  Parameter iSigma2       ( "iSigma2"      ,  0.00733001, .1 );
  Parameter rKstarm_1680  ( "rKstarm_1680" , -1.74556   , .1 );
  Parameter iKstarm_1680  ( "iKstarm_1680" ,  0.0505433 , .1 );
  Parameter rNonReson     ( "rNonReson"    ,  1.47537   , .1 );
  Parameter iNonReson     ( "iNonReson"    ,  0.904678  , .1 );

  //tau        .fix();
  //x          .fix();
  //y          .fix();

  mKStarc      .fix();
  wKStarc      .fix();
  mRho         .fix();
  wRho         .fix();
  mOmega       .fix();
  wOmega       .fix();
  mF0_980      .fix();
  wF0_980      .fix();
  mF0_1370     .fix();
  wF0_1370     .fix();
  mF2_1270     .fix();
  wF2_1270     .fix();
  mK0Starc_1430.fix();
  wK0Starc_1430.fix();
  mK2Starc_1430.fix();
  wK2Starc_1430.fix();
  mSigma       .fix();
  wSigma       .fix();
  mSigma2      .fix();
  wSigma2      .fix();
  mKStarm_1680 .fix();
  wKStarm_1680 .fix();

  tau          .setLimits(  0.39, 0.43 );
  x            .fix();
  y            .fix();
//   x            .setLimits( -0.01, 0.01 );
//   y            .setLimits( -0.01, 0.01 );

  rKstarm      .setLimits( -5., 5. );
  iKstarm      .setLimits( -5., 5. );
  rKstarp      .setLimits( -5., 5. );
  iKstarp      .setLimits( -5., 5. );
  rRho         .fix();
  iRho         .fix();
  rOmega       .setLimits( -5., 5. );
  iOmega       .setLimits( -5., 5. );
  rF0_980      .setLimits( -5., 5. );
  iF0_980      .setLimits( -5., 5. );
  rF0_1370     .setLimits( -5., 5. );
  iF0_1370     .setLimits( -5., 5. );
  rF2_1270     .setLimits( -5., 5. );
  iF2_1270     .setLimits( -5., 5. );
  rK0starm_1430.setLimits( -5., 5. );
  iK0starm_1430.setLimits( -5., 5. );
  rK0starp_1430.setLimits( -5., 5. );
  iK0starp_1430.setLimits( -5., 5. );
  rK2starm_1430.setLimits( -5., 5. );
  iK2starm_1430.setLimits( -5., 5. );
  rK2starp_1430.setLimits( -5., 5. );
  iK2starp_1430.setLimits( -5., 5. );
  rSigma       .setLimits( -5., 5. );
  iSigma       .setLimits( -5., 5. );
  rSigma2      .setLimits( -5., 5. );
  iSigma2      .setLimits( -5., 5. );
  rKstarm_1680 .setLimits( -5., 5. );
  iKstarm_1680 .setLimits( -5., 5. );
  rNonReson    .setLimits( -5., 5. );
  iNonReson    .setLimits( -5., 5. );

  std::vector< Variable  > vars;
  std::vector< Parameter > pars;

  vars.push_back( m2AB );
  vars.push_back( m2AC );
  vars.push_back( m2BC );
  vars.push_back( time );

  pars.push_back( tau         );
  pars.push_back( x           );
  pars.push_back( y           );

  // Masses and widths.
  pars.push_back( mKStarc       );
  pars.push_back( wKStarc       );
  pars.push_back( mRho          );
  pars.push_back( wRho          );
  pars.push_back( mOmega        );
  pars.push_back( wOmega        );
  pars.push_back( mF0_980       );
  pars.push_back( wF0_980       );
  pars.push_back( mF0_1370      );
  pars.push_back( wF0_1370      );
  pars.push_back( mF2_1270      );
  pars.push_back( wF2_1270      );
  pars.push_back( mK0Starc_1430 );
  pars.push_back( wK0Starc_1430 );
  pars.push_back( mK2Starc_1430 );
  pars.push_back( wK2Starc_1430 );
  pars.push_back( mSigma        );
  pars.push_back( wSigma        );
  pars.push_back( mSigma2       );
  pars.push_back( wSigma2       );
  pars.push_back( mKStarm_1680  );
  pars.push_back( wKStarm_1680  );

  // Real and imaginary part of the isobar coefficients.
  pars.push_back( rKstarm       );
  pars.push_back( iKstarm       );
  pars.push_back( rKstarp       );
  pars.push_back( iKstarp       );
  pars.push_back( rRho          );
  pars.push_back( iRho          );
  pars.push_back( rOmega        );
  pars.push_back( iOmega        );
  pars.push_back( rF0_980       );
  pars.push_back( iF0_980       );
  pars.push_back( rF0_1370      );
  pars.push_back( iF0_1370      );
  pars.push_back( rF2_1270      );
  pars.push_back( iF2_1270      );
  pars.push_back( rK0starm_1430 );
  pars.push_back( iK0starm_1430 );
  pars.push_back( rK0starp_1430 );
  pars.push_back( iK0starp_1430 );
  pars.push_back( rK2starm_1430 );
  pars.push_back( iK2starm_1430 );
  pars.push_back( rK2starp_1430 );
  pars.push_back( iK2starp_1430 );
  pars.push_back( rSigma        );
  pars.push_back( iSigma        );
  pars.push_back( rSigma2       );
  pars.push_back( iSigma2       );
  pars.push_back( rKstarm_1680  );
  pars.push_back( iKstarm_1680  );
  pars.push_back( rNonReson     );
  pars.push_back( iNonReson     );

  // Definition of the pdf.
  DalitzD0mix dalitz( vars, pars );

  // Definition of the minimizer from the pdf.
  Nll nll( dalitz, data );

  // Compute the minimum.
  FunctionMinimum minimum = nll.minimize();

#ifdef MPI_ON
  if ( rank == 0 )
    {
      std::cout << minimum << std::endl;
      std::cerr << "Processing time: " << double( clock() - initClock ) / double( CLOCKS_PER_SEC ) << " s" << std::endl;
    }

  MPI::Finalize();
#else
  std::cout << minimum << std::endl;
  std::cerr << "Processing time: " << double( clock() - initClock ) / double( CLOCKS_PER_SEC ) << " s" << std::endl;
#endif

  return 0;
}
