
#include <cmath>
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
  {
    return 0.0;
  }
  if( ( x > 1 ) && ( x > a ) )
  {
    return 1.0 - gamma_q( a, x );
  }
  ax = a * std::log( x ) - x - lngamma( a );
  if( ax < -709.78271289338399 )
  {
    return 0.0;
  }
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
  {
    return 1.0;
  }
  if( ( x < 1.0 ) || ( x < a ) )
  {
    return 1.0 - gamma_p( a, x );
  }
  ax = a * std::log( x ) - x - lngamma( a );
  if( ax < -709.78271289338399 )
  {
    return 0.0;
  }
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
  unsigned i;
  double logpi;
  double ls2pi;

  logpi = 1.14472988584940017414;
  ls2pi = 0.91893853320467274178;
  if( x < -34.0 )
  {
    q = -x;
    w = lngamma( q );
    p = std::floor( q );
    i = (unsigned) std::floor( p + 0.5 );
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
  if( x < 0.5 )
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
  if( x >= 10.0 )
  {
    return s;
  }
  return s * ( 1.0 - erfc( x ) );
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

