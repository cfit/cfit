
#include <iostream> // TEMPORAL
#include <iterator> // TEMPORAL
// FALTA CONSTRUIR EL PARVECT QUAN S'OPERA SOBRE UN OBJECTE PDF.

#include <sstream>
#include <vector>
#include <stack>
#include <algorithm>

#include <cfit/pdfmodel.hh>
#include <cfit/pdf.hh>

#include <cfit/pdfops.hh>


std::string int2str( int value )
{
  std::ostringstream str;
  str << value;
  return str.str();
}

/*
class SameVariable
{
  typedef std::map< std::string, Variable > varMap;
  bool operator() ( const varMap& x, const varMap& y ) const
  {
    
    return x.first == y.first;
  }
};
*/



Pdf& operator+( const PdfModel& left, const PdfModel& right )
{
  Pdf* pdf = new Pdf();

  // Add the operands in the expression.
  pdf->_expression += "mm+";

  // Build the variables map.
  pdf->_vars.insert( left ._vars.begin(), left ._vars.end() );
  pdf->_vars.insert( right._vars.begin(), right._vars.end() );

  // Build the parameters map.
  pdf->_pars.insert( left ._pars.begin(), left ._pars.end() );
  pdf->_pars.insert( right._pars.begin(), right._pars.end() );

  // Build the pdfs vector.
  pdf->_pdfVect.push_back( const_cast<PdfModel*>( &left  ) );
  pdf->_pdfVect.push_back( const_cast<PdfModel*>( &right ) );

  return *pdf;
}


Pdf& operator*( const PdfModel& left, const PdfModel& right )
{
  Pdf* pdf = new Pdf();

  // Add the operands in the expression.
  pdf->_expression += "mm*";

  // Build the variables map.
  pdf->_vars.insert( left ._vars.begin(), left ._vars.end() );
  pdf->_vars.insert( right._vars.begin(), right._vars.end() );

  // Build the parameters map.
  pdf->_pars.insert( left ._pars.begin(), left ._pars.end() );
  pdf->_pars.insert( right._pars.begin(), right._pars.end() );

  // Build the pdfs vector.
  pdf->_pdfVect.push_back( const_cast<PdfModel*>( &left  ) );
  pdf->_pdfVect.push_back( const_cast<PdfModel*>( &right ) );

  return *pdf;
}


Pdf& operator+( const Pdf& left, const PdfModel& right )
{
  Pdf* pdf = new Pdf();

  // Add the operands in the expression.
  pdf->_expression += left._expression + "m+";

  // Build the variables map.
  pdf->_vars.insert( left ._vars.begin(), left ._vars.end() );
  pdf->_vars.insert( right._vars.begin(), right._vars.end() );

  // Build the parameters map.
  pdf->_pars.insert( left ._pars.begin(), left ._pars.end() );
  pdf->_pars.insert( right._pars.begin(), right._pars.end() );

  // Build the parameters vector.
  pdf->_parVect.insert( pdf->_parVect.end(), left._parVect.begin(), left._parVect.end() );

  // Build the pdfs vector.
  pdf->_pdfVect.insert( pdf->_pdfVect.end(), left._pdfVect.begin(), left._pdfVect.end() );
  pdf->_pdfVect.push_back( const_cast<PdfModel*>( &right ) );

  return *pdf;
}


Pdf& operator*( const Pdf& left, const PdfModel& right )
{
  Pdf* pdf = new Pdf();

  // Add the operands in the expression.
  pdf->_expression += left._expression + "m*";

  // Build the variables map.
  pdf->_vars.insert( left ._vars.begin(), left ._vars.end() );
  pdf->_vars.insert( right._vars.begin(), right._vars.end() );

  // Build the parameters map.
  pdf->_pars.insert( left ._pars.begin(), left ._pars.end() );
  pdf->_pars.insert( right._pars.begin(), right._pars.end() );

  // Build the parameters vector.
  pdf->_parVect.insert( pdf->_parVect.end(), left._parVect.begin(), left._parVect.end() );

  // Build the pdfs vector.
  pdf->_pdfVect.insert( pdf->_pdfVect.end(), left._pdfVect.begin(), left._pdfVect.end() );
  pdf->_pdfVect.push_back( const_cast<PdfModel*>( &right ) );

  return *pdf;
}


Pdf& operator+( const Pdf& left, const Pdf& right )
{
  Pdf* pdf = new Pdf();

  // Add the operands in the expression.
  pdf->_expression = left._expression + right._expression + "+";

  // Build the variables map.
  pdf->_vars.insert( left ._vars.begin(), left ._vars.end() );
  pdf->_vars.insert( right._vars.begin(), right._vars.end() );

  // Build the parameters map.
  pdf->_pars.insert( left ._pars.begin(), left ._pars.end() );
  pdf->_pars.insert( right._pars.begin(), right._pars.end() );

  // Build the parameters vector.
  pdf->_parVect.insert( pdf->_parVect.end(), left ._parVect.begin(), left ._parVect.end() );
  pdf->_parVect.insert( pdf->_parVect.end(), right._parVect.begin(), right._parVect.end() );

  // Build the pdfs vector.
  pdf->_pdfVect.insert( pdf->_pdfVect.end(), left ._pdfVect.begin(), left ._pdfVect.end() );
  pdf->_pdfVect.insert( pdf->_pdfVect.end(), right._pdfVect.begin(), right._pdfVect.end() );

  return *pdf;
}


Pdf& operator*( const Pdf& left, const Pdf& right )
{
  Pdf* pdf = new Pdf();

  // Add the operands in the expression.
  pdf->_expression = left._expression + right._expression + "*";

  // Build the variables map.
  pdf->_vars.insert( left ._vars.begin(), left ._vars.end() );
  pdf->_vars.insert( right._vars.begin(), right._vars.end() );

  // Build the parameters map.
  pdf->_pars.insert( left ._pars.begin(), left ._pars.end() );
  pdf->_pars.insert( right._pars.begin(), right._pars.end() );

  // Build the parameters vector.
  pdf->_parVect.insert( pdf->_parVect.end(), left ._parVect.begin(), left ._parVect.end() );
  pdf->_parVect.insert( pdf->_parVect.end(), right._parVect.begin(), right._parVect.end() );

  // Build the pdfs vector.
  pdf->_pdfVect.insert( pdf->_pdfVect.end(), left ._pdfVect.begin(), left ._pdfVect.end() );
  pdf->_pdfVect.insert( pdf->_pdfVect.end(), right._pdfVect.begin(), right._pdfVect.end() );

  return *pdf;
}


Pdf& operator+( const Parameter& left, const PdfModel& right )
{
  Pdf* pdf = new Pdf();

  // Add the operands in the expression.
  pdf->_expression += "pm+";

  // Build the variables map.
  pdf->_vars.insert( right._vars.begin(), right._vars.end() );

  // Build the parameters map.
  pdf->_pars[ left.name() ] = left;
  pdf->_pars.insert( right._pars.begin(), right._pars.end() );

  // Build the parameters vector.
  pdf->_parVect.push_back( left.name() );

  // Build the pdfs vector.
  pdf->_pdfVect.push_back( const_cast<PdfModel*>( &right ) );

  return *pdf;
}


Pdf& operator*( const Parameter& left, const PdfModel& right )
{
  Pdf* pdf = new Pdf();

  // Add the operands in the expression.
  pdf->_expression += "pm*";

  // Build the variables map.
  pdf->_vars.insert( right._vars.begin(), right._vars.end() );

  // Build the parameters map.
  pdf->_pars[ left.name() ] = left;
  pdf->_pars.insert( right._pars.begin(), right._pars.end() );

  // Build the parameters vector.
  pdf->_parVect.push_back( left.name() );

  // Build the pdfs vector.
  pdf->_pdfVect.push_back( const_cast<PdfModel*>( &right ) );

  return *pdf;
}


Pdf& operator^( const Pdf& left, const Pdf& right ) throw( PdfException )
{
  // Find the sets of variables that are common in all the
  //    models of the left and the right.
  std::vector< std::string > commonLeft  = left .commonVars();
  std::vector< std::string > commonRight = right.commonVars();

  // Make sure that one of the previous sets contains the other.
  std::vector< std::string > common;
  std::set_intersection( commonLeft.begin() , commonLeft.end() ,
			 commonRight.begin(), commonRight.end(),
			 std::back_inserter( common ) );

  const Pdf* resolution;
  const Pdf* pdf;

  if ( common == commonLeft )
    {
      resolution = &left;
      pdf        = &right;
    }
  else if ( common == commonRight )
    {
      resolution = &right;
      pdf        = &left;
    }
  else
    throw PdfException( "Error computing convolution: left and right operands do not share enough variables." );

  // ALERTA: MANCA FER TOTA LA IMPLEMENTACIÓ DE LA CONVOLUCIÓ.

  return *( new Pdf() );
}
