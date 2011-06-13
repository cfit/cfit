
#include <iostream> // TEMPORAL
#include <iterator> // TEMPORAL
// FALTA CONSTRUIR EL PARVECT QUAN S'OPERA SOBRE UN OBJECTE PDF.

#include <sstream>
#include <vector>
#include <stack>
#include <algorithm>

#include <cfit/pdfmodel.hh>
#include <cfit/pdf.hh>

#include <cfit/parameter.hh>
#include <cfit/parameterexpr.hh>
#include <cfit/operation.hh>


// Operations with two parameters.
const ParameterExpr operator+( const Parameter& left, const Parameter& right )
{
  return ParameterExpr( left, right, Operation::plus );
}

const ParameterExpr operator-( const Parameter& left, const Parameter& right )
{
  return ParameterExpr( left, right, Operation::minus );
}

const ParameterExpr operator*( const Parameter& left, const Parameter& right )
{
  return ParameterExpr( left, right, Operation::mult );
}

const ParameterExpr operator/( const Parameter& left, const Parameter& right )
{
  return ParameterExpr( left, right, Operation::div );
}


// Operations with a constant and a parameter.
const ParameterExpr operator+( const double& left, const Parameter& right )
{
  return ParameterExpr( left, right, Operation::plus );
}

const ParameterExpr operator-( const double& left, const Parameter& right )
{
  return ParameterExpr( left, right, Operation::minus );
}

const ParameterExpr operator*( const double& left, const Parameter& right )
{
  return ParameterExpr( left, right, Operation::mult );
}

const ParameterExpr operator/( const double& left, const Parameter& right )
{
  return ParameterExpr( left, right, Operation::div );
}


// Operations with a parameter and a constant.
const ParameterExpr operator+( const Parameter& left, const double& right )
{
  return ParameterExpr( left, right, Operation::plus );
}

const ParameterExpr operator-( const Parameter& left, const double& right )
{
  return ParameterExpr( left, right, Operation::minus );
}

const ParameterExpr operator*( const Parameter& left, const double& right )
{
  return ParameterExpr( left, right, Operation::mult );
}

const ParameterExpr operator/( const Parameter& left, const double& right )
{
  return ParameterExpr( left, right, Operation::div );
}


// Operations with two parameter expressions.
const ParameterExpr operator+( const ParameterExpr& left, const ParameterExpr& right )
{
  return ParameterExpr( left, right, Operation::plus );
}

const ParameterExpr operator-( const ParameterExpr& left, const ParameterExpr& right )
{
  return ParameterExpr( left, right, Operation::minus );
}

const ParameterExpr operator*( const ParameterExpr& left, const ParameterExpr& right )
{
  return ParameterExpr( left, right, Operation::mult );
}

const ParameterExpr operator/( const ParameterExpr& left, const ParameterExpr& right )
{
  return ParameterExpr( left, right, Operation::div );
}


// Operations with a constant and a parameter expression.
const ParameterExpr operator+( const double& left, const ParameterExpr& right )
{
  return ParameterExpr( left, right, Operation::plus );
}

const ParameterExpr operator-( const double& left, const ParameterExpr& right )
{
  return ParameterExpr( left, right, Operation::minus );
}

const ParameterExpr operator*( const double& left, const ParameterExpr& right )
{
  return ParameterExpr( left, right, Operation::mult );
}

const ParameterExpr operator/( const double& left, const ParameterExpr& right )
{
  return ParameterExpr( left, right, Operation::div );
}


// Operations with a parameter expression and a constant.
const ParameterExpr operator+( const ParameterExpr& left, const double& right )
{
  return ParameterExpr( left, right, Operation::plus );
}

const ParameterExpr operator-( const ParameterExpr& left, const double& right )
{
  return ParameterExpr( left, right, Operation::minus );
}

const ParameterExpr operator*( const ParameterExpr& left, const double& right )
{
  return ParameterExpr( left, right, Operation::mult );
}

const ParameterExpr operator/( const ParameterExpr& left, const double& right )
{
  return ParameterExpr( left, right, Operation::div );
}


// Operations with a parameter and a parameter expression.
const ParameterExpr operator+( const Parameter& left, const ParameterExpr& right )
{
  return ParameterExpr( left, right, Operation::plus );
}

const ParameterExpr operator-( const Parameter& left, const ParameterExpr& right )
{
  return ParameterExpr( left, right, Operation::minus );
}

const ParameterExpr operator*( const Parameter& left, const ParameterExpr& right )
{
  return ParameterExpr( left, right, Operation::mult );
}

const ParameterExpr operator/( const Parameter& left, const ParameterExpr& right )
{
  return ParameterExpr( left, right, Operation::div );
}


// Operations with a parameter expression and a parameter.
const ParameterExpr operator+( const ParameterExpr& left, const Parameter& right )
{
  return ParameterExpr( left, right, Operation::plus );
}

const ParameterExpr operator-( const ParameterExpr& left, const Parameter& right )
{
  return ParameterExpr( left, right, Operation::minus );
}

const ParameterExpr operator*( const ParameterExpr& left, const Parameter& right )
{
  return ParameterExpr( left, right, Operation::mult );
}

const ParameterExpr operator/( const ParameterExpr& left, const Parameter& right )
{
  return ParameterExpr( left, right, Operation::div );
}



// Unary minus.
const ParameterExpr operator-( const Parameter& par )
{
  return ParameterExpr( par, Operation::minus );
}

const ParameterExpr operator-( const ParameterExpr& par )
{
  return ParameterExpr( par, Operation::minus );
}


// Binary operations with a parameter.
const ParameterExpr pow( const Parameter& left, const double& right )
{
  return ParameterExpr( left, right, Operation::pow );
}

const ParameterExpr pow( const ParameterExpr& left, const double& right )
{
  return ParameterExpr( left, right, Operation::pow );
}


// Unary operations with a parameter.
const ParameterExpr exp( const Parameter& par )
{
  return ParameterExpr( par, Operation::exp );
}

const ParameterExpr log( const Parameter& par )
{
  return ParameterExpr( par, Operation::log );
}

const ParameterExpr sin( const Parameter& par )
{
  return ParameterExpr( par, Operation::sin );
}

const ParameterExpr cos( const Parameter& par )
{
  return ParameterExpr( par, Operation::cos );
}

const ParameterExpr tan( const Parameter& par )
{
  return ParameterExpr( par, Operation::tan );
}


const ParameterExpr exp( const ParameterExpr& par )
{
  return ParameterExpr( par, Operation::exp );
}

const ParameterExpr log( const ParameterExpr& par )
{
  return ParameterExpr( par, Operation::log );
}

const ParameterExpr sin( const ParameterExpr& par )
{
  return ParameterExpr( par, Operation::sin );
}

const ParameterExpr cos( const ParameterExpr& par )
{
  return ParameterExpr( par, Operation::cos );
}

const ParameterExpr tan( const ParameterExpr& par )
{
  return ParameterExpr( par, Operation::tan );
}




// Operations with two pdf models.
const Pdf operator+( const PdfModel& left, const PdfModel& right )
{
  // Cannot add two pdfs that do not depend on exactly the same variables.
  if ( left.varNames() != right.varNames() )
    throw PdfException( "Cannot add two pdfs that do not depend on the same variables." );

  return Pdf( left, right, Operation::plus );
}

const Pdf operator*( const PdfModel& left, const PdfModel& right )
{
  // Cannot multiply two pdfs that share some variable.
  std::vector< std::string > lVars = left. varNames();
  std::vector< std::string > rVars = right.varNames();
  std::vector< std::string > intersect;
  std::set_intersection( lVars.begin(), lVars.end(),
			 rVars.begin(), rVars.end(),
			 std::back_inserter( intersect ) );

  if ( ! intersect.empty() )
    throw PdfException( "Cannot multiply two pdfs that depend on some common variable." );

  return Pdf( left, right, Operation::mult );
}

// Operations with a pdf model and a pdf expression.
const Pdf operator+( const PdfModel& left, const Pdf& right )
{
  // Cannot add two pdfs that do not depend on exactly the same variables.
  if ( left.varNames() != right.varNames() )
    throw PdfException( "Cannot add two pdfs that do not depend on the same variables." );

  return Pdf( left, right, Operation::plus );
}

const Pdf operator*( const PdfModel& left, const Pdf& right )
{
  // Cannot multiply two pdfs that share some variable.
  std::vector< std::string > lVars = left. varNames();
  std::vector< std::string > rVars = right.varNames();
  std::vector< std::string > intersect;
  std::set_intersection( lVars.begin(), lVars.end(),
			 rVars.begin(), rVars.end(),
			 std::back_inserter( intersect ) );

  if ( ! intersect.empty() )
    throw PdfException( "Cannot multiply two pdfs that depend on some common variable." );

  return Pdf( left, right, Operation::mult );
}

// Operations with a pdf expression and a pdf model.
const Pdf operator+( const Pdf& left, const PdfModel& right )
{
  // Cannot add two pdfs that do not depend on exactly the same variables.
  if ( left.varNames() != right.varNames() )
    throw PdfException( "Cannot add two pdfs that do not depend on the same variables." );

  return Pdf( left, right, Operation::plus );
}

const Pdf operator*( const Pdf& left, const PdfModel& right )
{
  // Cannot multiply two pdfs that share some variable.
  std::vector< std::string > lVars = left. varNames();
  std::vector< std::string > rVars = right.varNames();
  std::vector< std::string > intersect;
  std::set_intersection( lVars.begin(), lVars.end(),
			 rVars.begin(), rVars.end(),
			 std::back_inserter( intersect ) );

  if ( ! intersect.empty() )
    throw PdfException( "Cannot multiply two pdfs that depend on some common variable." );

  return Pdf( left, right, Operation::mult );
}

// Operations with two pdf expressions.
const Pdf operator+( const Pdf& left, const Pdf& right )
{
  // Cannot add two pdfs that do not depend on exactly the same variables.
  if ( left.varNames() != right.varNames() )
    throw PdfException( "Cannot add two pdfs that do not depend on the same variables." );

  return Pdf( left, right, Operation::plus );
}

const Pdf operator*( const Pdf& left, const Pdf& right )
{
  // Cannot multiply two pdfs that share some variable.
  std::vector< std::string > lVars = left. varNames();
  std::vector< std::string > rVars = right.varNames();
  std::vector< std::string > intersect;
  std::set_intersection( lVars.begin(), lVars.end(),
			 rVars.begin(), rVars.end(),
			 std::back_inserter( intersect ) );

  if ( ! intersect.empty() )
    throw PdfException( "Cannot multiply two pdfs that depend on some common variable." );

  return Pdf( left, right, Operation::mult );
}


// Operations with a parameter and a model.
const Pdf operator*( const Parameter& left, const PdfModel& right )
{
  return Pdf( left, right, Operation::mult );
}

const Pdf operator*( const PdfModel& left, const Parameter& right )
{
  return Pdf( left, right, Operation::mult );
}

const Pdf operator/( const PdfModel& left, const Parameter& right )
{
  return Pdf( left, right, Operation::div );
}


// Operations with a parameter expression and a model.
const Pdf operator*( const ParameterExpr& left, const PdfModel& right )
{
  return Pdf( left, right, Operation::mult );
}

const Pdf operator*( const PdfModel& left, const ParameterExpr& right )
{
  return Pdf( left, right, Operation::mult );
}

const Pdf operator/( const PdfModel& left, const ParameterExpr& right )
{
  return Pdf( left, right, Operation::div );
}


// Operations with a constant and a model.
const Pdf operator*( const double& left, const PdfModel& right )
{
  return Pdf( left, right, Operation::mult );
}

const Pdf operator*( const PdfModel& left, const double& right )
{
  return Pdf( left, right, Operation::mult );
}

const Pdf operator/( const PdfModel& left, const double& right )
{
  return Pdf( left, right, Operation::div );
}


// Operations with a parameter and a pdf expression.
const Pdf operator*( const Parameter& left, const Pdf& right )
{
  return Pdf( left, right, Operation::mult );
}

const Pdf operator*( const Pdf& left, const Parameter& right )
{
  return Pdf( left, right, Operation::mult );
}

const Pdf operator/( const Pdf& left, const Parameter& right )
{
  return Pdf( left, right, Operation::div );
}


// Operations with a parameter expression and a pdf expression.
const Pdf operator*( const ParameterExpr& left, const Pdf& right )
{
  return Pdf( left, right, Operation::mult );
}

const Pdf operator*( const Pdf& left, const ParameterExpr& right )
{
  return Pdf( left, right, Operation::mult );
}

const Pdf operator/( const Pdf& left, const ParameterExpr& right )
{
  return Pdf( left, right, Operation::div );
}


// Operations with a constant and a pdf expression.
const Pdf operator*( const double& left, const Pdf& right )
{
  return Pdf( left, right, Operation::mult );
}

const Pdf operator*( const Pdf& left, const double& right )
{
  return Pdf( left, right, Operation::mult );
}

const Pdf operator/( const Pdf& left, const double& right )
{
  return Pdf( left, right, Operation::div );
}


const Pdf operator^( const Pdf& left, const Pdf& right ) throw( PdfException )
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

