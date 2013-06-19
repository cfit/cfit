
#include <cmath>
#include <limits>
#include <cfit/math.hh>


/*************************************************************************
Copyright (c) Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
>>> END OF LICENSE >>>
*************************************************************************/


/*************************************************************************
Incomplete gamma integral

The function is defined by

                          x
                           -
                  1       | |  -t  a-1
 igam(a,x)  =   -----     |   e   t   dt.
                 -      | |
                | (a)    -
                          0


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30       200000       3.6e-14     2.9e-15
   IEEE      0,100      300000       9.9e-14     1.5e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
const double Math::gamma_p( const double& a, const double& x )
{
  double igammaepsilon;
  double ans;
  double ax;
  double c;
  double r;

  igammaepsilon = 0.000000000000001;
  if( ( x <= 0 ) || ( a <= 0 ) )
    return 0.0;

  if( ( x > 1 ) && ( x > a ) )
    return 1.0 - gamma_q( a, x );

  ax = a * std::log( x ) - x - lngamma( a );
  if( ax < -709.78271289338399 )
    return 0.0;

  ax = std::exp( ax );
  r = a;
  c = 1;
  ans = 1;
  do
  {
    r = r+1;
    c = c*x/r;
    ans = ans+c;
  }
  while( c/ans > igammaepsilon );

  return ans*ax/a;
}




/*************************************************************************
Complemented incomplete gamma integral

The function is defined by


 igamc(a,x)   =   1 - igam(a,x)

                           inf.
                             -
                    1       | |  -t  a-1
              =   -----     |   e   t   dt.
                   -      | |
                  | (a)    -
                            x


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

Tested at random a, x.
               a         x                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
   IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
const double Math::gamma_q( const double& a, const double& x )
{
  double igammaepsilon;
  double igammabignumber;
  double igammabignumberinv;
  double ans;
  double ax;
  double c;
  double yc;
  double r;
  double t;
  double y;
  double z;
  double pk;
  double pkm1;
  double pkm2;
  double qk;
  double qkm1;
  double qkm2;


  igammaepsilon = 0.000000000000001;
  igammabignumber = 4503599627370496.0;
  igammabignumberinv = 2.22044604925031308085*0.0000000000000001;
  if( (x <= 0.0 ) || ( a <= 0 ) )
    return 1.0;

  if( ( x < 1.0 ) || ( x < a ) )
    return 1.0 - gamma_p( a, x );

  ax = a * std::log( x ) - x - lngamma( a );
  if( ax < -709.78271289338399 )
    return 0.0;

  ax = std::exp( ax );
  y = 1-a;
  z = x+y+1;
  c = 0;
  pkm2 = 1;
  qkm2 = x;
  pkm1 = x+1;
  qkm1 = z*x;
  ans = pkm1/qkm1;
  do
  {
    c = c+1;
    y = y+1;
    z = z+2;
    yc = y*c;
    pk = pkm1*z-pkm2*yc;
    qk = qkm1*z-qkm2*yc;
    if( qk != 0 )
    {
      r = pk/qk;
      t = std::fabs( (ans-r)/r );
      ans = r;
    }
    else
    {
      t = 1;
    }
    pkm2 = pkm1;
    pkm1 = pk;
    qkm2 = qkm1;
    qkm1 = qk;
    if( std::fabs( pk ) > igammabignumber )
    {
      pkm2 = pkm2*igammabignumberinv;
      pkm1 = pkm1*igammabignumberinv;
      qkm2 = qkm2*igammabignumberinv;
      qkm1 = qkm1*igammabignumberinv;
    }
  }
  while( t > igammaepsilon );

  return ans*ax;
}




const double Math::gammastirf( const double& x )
{
  double y;
  double w;
  double v;
  double stir;

  w = 1/x;
  stir = 7.87311395793093628397E-4;
  stir = -2.29549961613378126380E-4+w*stir;
  stir = -2.68132617805781232825E-3+w*stir;
  stir = 3.47222221605458667310E-3+w*stir;
  stir = 8.33333333333482257126E-2+w*stir;
  w = 1+w*stir;
  y = std::exp( x );
  if( x > 143.01608 )
  {
    v = std::pow( x, 0.5*x-0.25 );
    y = v*(v/y);
  }
  else
  {
    y = std::pow( x, x-0.5 )/y;
  }
  return 2.50662827463100050242*y*w;
}


/*************************************************************************
Gamma function

Input parameters:
    X   -   argument

Domain:
    0 < X < 171.6
    -170 < X < 0, X is not an integer.

Relative error:
 arithmetic   domain     # trials      peak         rms
    IEEE    -170,-33      20000       2.3e-15     3.3e-16
    IEEE     -33,  33     20000       9.4e-16     2.2e-16
    IEEE      33, 171.6   20000       2.3e-15     3.2e-16

Cephes Math Library Release 2.8:  June, 2000
Original copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
Translated to AlgoPascal by Bochkanov Sergey (2005, 2006, 2007).
*************************************************************************/
const double Math::gamma( double x )
{
  double p;
  double pp;
  double q;
  double qq;
  double z;
  unsigned i;
  double sgngam;

  sgngam = 1.0;
  q = std::fabs( x );
  if( q > 33.0 )
  {
    if( x < 0.0 )
    {
      p = std::floor( q );
      i = (unsigned) std::floor( p + 0.5 );
      if( i%2==0 )
      {
        sgngam = -1.0;
      }
      z = q-p;
      if( z > 0.5 )
      {
        p = p+1;
        z = q-p;
      }
      z = q * std::sin( M_PI * z );
      z = std::fabs( z );
      z = M_PI / ( z * gammastirf( q ) );
    }
    else
    {
      z = gammastirf( x );
    }
    return sgngam*z;
  }
  z = 1;
  while( x >= 3.0 )
  {
    x = x-1;
    z = z*x;
  }
  while( x < 0 )
  {
    if( x > -0.000000001 )
    {
      return z/((1+0.5772156649015329*x)*x);
    }
    z = z/x;
    x = x+1;
  }
  while( x < 2.0 )
  {
    if( x < 0.000000001 )
    {
      return z/((1+0.5772156649015329*x)*x);
    }
    z = z/x;
    x = x+1.0;
  }
  if( x == 2.0 )
  {
    return z;
  }
  x = x-2.0;
  pp = 1.60119522476751861407E-4;
  pp = 1.19135147006586384913E-3+x*pp;
  pp = 1.04213797561761569935E-2+x*pp;
  pp = 4.76367800457137231464E-2+x*pp;
  pp = 2.07448227648435975150E-1+x*pp;
  pp = 4.94214826801497100753E-1+x*pp;
  pp = 9.99999999999999996796E-1+x*pp;
  qq = -2.31581873324120129819E-5;
  qq = 5.39605580493303397842E-4+x*qq;
  qq = -4.45641913851797240494E-3+x*qq;
  qq = 1.18139785222060435552E-2+x*qq;
  qq = 3.58236398605498653373E-2+x*qq;
  qq = -2.34591795718243348568E-1+x*qq;
  qq = 7.14304917030273074085E-2+x*qq;
  qq = 1.00000000000000000320+x*qq;
  return z*pp/qq;
}






/*************************************************************************
Natural logarithm of gamma function

Input parameters:
    X       -   argument

Result:
    logarithm of the absolute value of the Gamma(X).

Output parameters:
    SgnGam  -   sign(Gamma(X))

Domain:
    0 < X < 2.55e305
    -2.55e305 < X < 0, X is not an integer.

ACCURACY:
arithmetic      domain        # trials     peak         rms
   IEEE    0, 3                 28000     5.4e-16     1.1e-16
   IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
The error criterion was relative when the function magnitude
was greater than one but absolute when it was less than one.

The following test used the relative error criterion, though
at certain points the relative error could be much higher than
indicated.
   IEEE    -200, -4             10000     4.8e-16     1.3e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
Translated to AlgoPascal by Bochkanov Sergey (2005, 2006, 2007).
*************************************************************************/
const double Math::lngamma( double x )
{
  double a;
  double b;
  double c;
  double p;
  double q;
  double u;
  double w;
  double z;
  double logpi;
  double ls2pi;

  logpi = 1.14472988584940017414;
  ls2pi = 0.91893853320467274178;
  if( x < -34.0 )
  {
    q = -x;
    w = lngamma( q );
    p = std::floor( q );
    z = q-p;
    if( z > 0.5 )
    {
      p = p+1;
      z = p-q;
    }
    z = q * std::sin( M_PI * z );
    return logpi - std::log( z ) - w;
  }
  if( x < 13.0 )
  {
    z = 1;
    p = 0;
    u = x;
    while( u >= 3.0 )
    {
      p = p-1;
      u = x+p;
      z = z*u;
    }
    while( u < 2.0 )
    {
      z = z/u;
      p = p+1;
      u = x+p;
    }
    if( z < 0.0 )
    {
      z = -z;
    }
    if( u == 2.0 )
    {
      return std::log( z );
    }
    p = p-2;
    x = x+p;
    b = -1378.25152569120859100;
    b = -38801.6315134637840924+x*b;
    b = -331612.992738871184744+x*b;
    b = -1162370.97492762307383+x*b;
    b = -1721737.00820839662146+x*b;
    b = -853555.664245765465627+x*b;
    c = 1;
    c = -351.815701436523470549+x*c;
    c = -17064.2106651881159223+x*c;
    c = -220528.590553854454839+x*c;
    c = -1139334.44367982507207+x*c;
    c = -2532523.07177582951285+x*c;
    c = -2018891.41433532773231+x*c;
    p = x*b/c;
    return std::log( z ) + p;
  }
  q = (x-0.5) * std::log( x ) - x + ls2pi;
  if( x > 100000000.0 )
  {
    return q;
  }
  p = 1/(x*x);
  if( x >= 1000.0 )
  {
    q = q+((7.9365079365079365079365*0.0001*p-2.7777777777777777777778*0.001)*p+0.0833333333333333333333)/x;
  }
  else
  {
    a = 8.11614167470508450300*0.0001;
    a = -5.95061904284301438324*0.0001+p*a;
    a = 7.93650340457716943945*0.0001+p*a;
    a = -2.77777777730099687205*0.001+p*a;
    a = 8.33333333333331927722*0.01+p*a;
    q = q+a/x;
  }
  return q;
}


/*************************************************************************
Error function

The integral is

                          x
                           -
                2         | |          2
  erf(x)  =  --------     |    exp( - t  ) dt.
             sqrt(pi)   | |
                         -
                          0

For 0 <= |x| < 1, erf(x) = x * P4(x**2)/Q5(x**2); otherwise
erf(x) = 1 - erfc(x).


ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,1         30000       3.7e-16     1.0e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
*************************************************************************/
const double Math::erf( const double& x )
{
  double absx;
  double xsq;
  double s;
  double p;
  double q;

  absx = std::fabs( x );
  s = ( x != 0 ) ? x / absx : 0.0; // Sign of x.
  if( absx < 0.5 )
  {
    xsq = x*x;
    p = 0.007547728033418631287834;
    p = 0.288805137207594084924010+xsq*p;
    p = 14.3383842191748205576712+xsq*p;
    p = 38.0140318123903008244444+xsq*p;
    p = 3017.82788536507577809226+xsq*p;
    p = 7404.07142710151470082064+xsq*p;
    p = 80437.3630960840172832162+xsq*p;
    q = 0.0;
    q = 1.00000000000000000000000+xsq*q;
    q = 38.0190713951939403753468+xsq*q;
    q = 658.070155459240506326937+xsq*q;
    q = 6379.60017324428279487120+xsq*q;
    q = 34216.5257924628539769006+xsq*q;
    q = 80437.3630960840172826266+xsq*q;
    return s*1.1283791670955125738961589031 * absx * p / q;
  }
  if( absx >= 10.0 )
  {
    return s;
  }
  return s * ( 1.0 - erfc( absx ) );
}


/*************************************************************************
Complementary error function

 1 - erf(x) =

                          inf.
                            -
                 2         | |          2
  erfc(x)  =  --------     |    exp( - t  ) dt
              sqrt(pi)   | |
                          -
                           x


For small x, erfc(x) = 1 - erf(x); otherwise rational
approximations are computed.


ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,26.6417   30000       5.7e-14     1.5e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
*************************************************************************/
const double Math::erfc( const double& x )
{
  double p;
  double q;

  if( x < 0 )
  {
    return 2.0 - erfc( -x );
  }
  if( x < 0.5 )
  {
    return 1.0 - erf( x );
  }
  if( x >= 10.0 )
  {
    return 0.0;
  }
  p = 0.0;
  p = 0.5641877825507397413087057563+x*p;
  p = 9.675807882987265400604202961+x*p;
  p = 77.08161730368428609781633646+x*p;
  p = 368.5196154710010637133875746+x*p;
  p = 1143.262070703886173606073338+x*p;
  p = 2320.439590251635247384768711+x*p;
  p = 2898.0293292167655611275846+x*p;
  p = 1826.3348842295112592168999+x*p;
  q = 1.0;
  q = 17.14980943627607849376131193+x*q;
  q = 137.1255960500622202878443578+x*q;
  q = 661.7361207107653469211984771+x*q;
  q = 2094.384367789539593790281779+x*q;
  q = 4429.612803883682726711528526+x*q;
  q = 6089.5424232724435504633068+x*q;
  q = 4958.82756472114071495438422+x*q;
  q = 1826.3348842295112595576438+x*q;
  return std::exp( -std::pow( x, 2 ) )*p/q;
}




/*************************************************************************
Inverse of the error function

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
*************************************************************************/
const double Math::inverf( const double& x )
{
  return invnormal( 0.5 * ( x + 1.0 ) ) / std::sqrt( 2.0 );
}




/*************************************************************************
Inverse of Normal distribution function

Returns the argument, x, for which the area under the
Gaussian probability density function (integrated from
minus infinity to x) is equal to y.


For small arguments 0 < y < exp(-2), the program computes
z = sqrt( -2.0 * log(y) );  then the approximation is
x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).
There are two rational functions P/Q, one for 0 < y < exp(-32)
and the other for y up to exp(-2).  For larger arguments,
w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).

ACCURACY:

                     Relative error:
arithmetic   domain        # trials      peak         rms
   IEEE     0.125, 1        20000       7.2e-16     1.3e-16
   IEEE     3e-308, 0.135   50000       4.6e-16     9.8e-17

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
*************************************************************************/
const double Math::invnormal( const double& y0 )
{
  double expm2;
  double s2pi;
  double x;
  double y;
  double z;
  double y2;
  double x0;
  double x1;
  bool code;
  double p0;
  double q0;
  double p1;
  double q1;
  double p2;
  double q2;

  expm2 = 0.13533528323661269189;
  s2pi = 2.50662827463100050242;

  if ( y0 <= 0.0 )
    return -std::numeric_limits< double >::infinity();

  if( y0 >= 1.0 )
    return std::numeric_limits< double >::infinity();

  code = true;
  y = y0;

  if( y > 1.0 - expm2 )
  {
    y = 1.0 - y;
    code = false;
  }

  if( y > expm2 )
  {
    y = y-0.5;
    y2 = y*y;
    p0 = -59.9633501014107895267;
    p0 = 98.0010754185999661536+y2*p0;
    p0 = -56.6762857469070293439+y2*p0;
    p0 = 13.9312609387279679503+y2*p0;
    p0 = -1.23916583867381258016+y2*p0;
    q0 = 1;
    q0 = 1.95448858338141759834+y2*q0;
    q0 = 4.67627912898881538453+y2*q0;
    q0 = 86.3602421390890590575+y2*q0;
    q0 = -225.462687854119370527+y2*q0;
    q0 = 200.260212380060660359+y2*q0;
    q0 = -82.0372256168333339912+y2*q0;
    q0 = 15.9056225126211695515+y2*q0;
    q0 = -1.18331621121330003142+y2*q0;
    x = y+y*y2*p0/q0;
    x = x*s2pi;

    return x;
  }

  x = std::sqrt( -2.0 * std::log( y ) );
  x0 = x - std::log( x ) / x;
  z = 1.0/x;
  if( x < 8.0 )
  {
    p1 = 4.05544892305962419923;
    p1 = 31.5251094599893866154+z*p1;
    p1 = 57.1628192246421288162+z*p1;
    p1 = 44.0805073893200834700+z*p1;
    p1 = 14.6849561928858024014+z*p1;
    p1 = 2.18663306850790267539+z*p1;
    p1 = -1.40256079171354495875*0.1+z*p1;
    p1 = -3.50424626827848203418*0.01+z*p1;
    p1 = -8.57456785154685413611*0.0001+z*p1;
    q1 = 1;
    q1 = 15.7799883256466749731+z*q1;
    q1 = 45.3907635128879210584+z*q1;
    q1 = 41.3172038254672030440+z*q1;
    q1 = 15.0425385692907503408+z*q1;
    q1 = 2.50464946208309415979+z*q1;
    q1 = -1.42182922854787788574*0.1+z*q1;
    q1 = -3.80806407691578277194*0.01+z*q1;
    q1 = -9.33259480895457427372*0.0001+z*q1;
    x1 = z*p1/q1;
  }
  else
  {
    p2 = 3.23774891776946035970;
    p2 = 6.91522889068984211695+z*p2;
    p2 = 3.93881025292474443415+z*p2;
    p2 = 1.33303460815807542389+z*p2;
    p2 = 2.01485389549179081538*0.1+z*p2;
    p2 = 1.23716634817820021358*0.01+z*p2;
    p2 = 3.01581553508235416007*0.0001+z*p2;
    p2 = 2.65806974686737550832*0.000001+z*p2;
    p2 = 6.23974539184983293730*0.000000001+z*p2;
    q2 = 1;
    q2 = 6.02427039364742014255+z*q2;
    q2 = 3.67983563856160859403+z*q2;
    q2 = 1.37702099489081330271+z*q2;
    q2 = 2.16236993594496635890*0.1+z*q2;
    q2 = 1.34204006088543189037*0.01+z*q2;
    q2 = 3.28014464682127739104*0.0001+z*q2;
    q2 = 2.89247864745380683936*0.000001+z*q2;
    q2 = 6.79019408009981274425*0.000000001+z*q2;
    x1 = z*p2/q2;
  }
  x = x0-x1;

  if( code )
    return -x;

  return x;
}




/*************************************************************************
Inverse of complemented imcomplete gamma integral

Given p, the function finds x such that

 igamc( a, x ) = p.

Starting with the approximate value

        3
 x = a t

 where

 t = 1 - d - ndtri(p) sqrt(d)

and

 d = 1/9a,

the routine performs up to 10 Newton iterations to find the
root of igamc(a,x) - p = 0.

ACCURACY:

Tested at random a, p in the intervals indicated.

               a        p                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,0.5       100000       1.0e-14     1.7e-15
   IEEE     0.01,0.5  0,0.5       100000       9.0e-14     3.4e-15
   IEEE    0.5,10000  0,0.5        20000       2.3e-13     3.8e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
const double Math::invgamma_q( const double& a, const double& y0 )
{
  double igammaepsilon;
  double iinvgammabignumber;
  double x0;
  double x1;
  double x;
  double yl;
  double yh;
  double y;
  double d;
  double lgm;
  double dithresh;
  int i;
  int dir;

  igammaepsilon = 0.000000000000001;
  iinvgammabignumber = 4503599627370496.0;
  x0 = iinvgammabignumber;
  yl = 0;
  x1 = 0;
  yh = 1;
  dithresh = 5.0 * igammaepsilon;
  d = 1.0 / ( 9.0 * a );
  y = 1.0 - d - invnormal( y0 ) * std::sqrt( d );
  x = a*y*y*y;
  lgm = lngamma( a );
  i = 0;
  while( i < 10 )
  {
    if( ( x > x0 ) || ( x < x1 ) )
    {
      d = 0.0625;
      break;
    }
    y = gamma_q( a, x );
    if( ( x < yl ) || ( y > yh ) )
    {
      d = 0.0625;
      break;
    }
    if( y < y0 )
    {
      x0 = x;
      yl = y;
    }
    else
    {
      x1 = x;
      yh = y;
    }
    d = ( a - 1.0 ) * std::log( x ) - x - lgm;
    if( d < -709.78271289338399 )
    {
      d = 0.0625;
      break;
    }
    d = -std::exp( d );
    d = (y-y0)/d;

    if( std::fabs( d/x ) < igammaepsilon )
      return x;

    x = x-d;
    i = i+1;
  }
  if( x0 == iinvgammabignumber )
  {
    if( x <= 0.0 )
      x = 1;

    while( x0 == iinvgammabignumber )
    {
      x = ( 1.0 + d ) * x;
      y = gamma_q( a, x );
      if( y < y0 )
      {
        x0 = x;
        yl = y;
        break;
      }
      d = d + d;
    }
  }
  d = 0.5;
  dir = 0;
  i = 0;
  while( i < 400 )
  {
    x = x1+d*(x0-x1);
    y = gamma_q( a, x );
    lgm = (x0-x1)/(x1+x0);
    if( std::fabs( lgm ) < dithresh )
      break;

    lgm = (y-y0)/y0;
    if( std::fabs( lgm ) < dithresh )
      break;

    if( x <= 0.0 )
      break;

    if( y >= y0 )
    {
      x1 = x;
      yh = y;
      if( dir < 0 )
      {
        dir = 0;
        d = 0.5;
      }
      else
      {
        if( dir > 1 )
        {
          d = 0.5*d+0.5;
        }
        else
        {
          d = (y0-yl)/(yh-yl);
        }
      }
      dir = dir+1;
    }
    else
    {
      x0 = x;
      yl = y;
      if( dir > 0 )
      {
        dir = 0.0;
        d = 0.5;
      }
      else
      {
        if( dir < -1 )
          d = 0.5*d;
        else
          d = (y0-yl)/(yh-yl);
      }
      dir = dir-1;
    }
    i = i+1;
  }

  return x;
}

