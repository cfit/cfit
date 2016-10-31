
#include <sstream>
#include <vector>
#include <stack>
#include <algorithm>

#include <cfit/pdfmodel.hh>
#include <cfit/pdfexpr.hh>

#include <cfit/parameter.hh>
#include <cfit/parameterexpr.hh>
#include <cfit/coef.hh>
#include <cfit/coefexpr.hh>
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

const ParameterExpr tanh( const Parameter& par )
{
  return ParameterExpr( par, Operation::tanh );
}

const ParameterExpr atanh( const Parameter& par )
{
  return ParameterExpr( par, Operation::atanh );
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

const ParameterExpr tanh( const ParameterExpr& par )
{
  return ParameterExpr( par, Operation::tanh );
}

const ParameterExpr atanh( const ParameterExpr& par )
{
  return ParameterExpr( par, Operation::atanh );
}




// Operations with two coefficients.
const CoefExpr operator+( const Coef& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const Coef& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const Coef& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const Coef& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::div );
}


// Operations with a constant and a coefficient.
const CoefExpr operator+( const double& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const double& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const double& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const double& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::div );
}


// Operations with a coefficient and a constant.
const CoefExpr operator+( const Coef& left, const double& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const Coef& left, const double& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const Coef& left, const double& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const Coef& left, const double& right )
{
  return CoefExpr( left, right, Operation::div );
}




// Operations with a parameter and a coefficient.
const CoefExpr operator+( const Parameter& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const Parameter& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const Parameter& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const Parameter& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::div );
}


// Operations with a coefficient and a parameter.
const CoefExpr operator+( const Coef& left, const Parameter& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const Coef& left, const Parameter& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const Coef& left, const Parameter& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const Coef& left, const Parameter& right )
{
  return CoefExpr( left, right, Operation::div );
}




// Operations with a parameter expression and a coefficient.
const CoefExpr operator+( const ParameterExpr& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const ParameterExpr& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const ParameterExpr& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const ParameterExpr& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::div );
}


// Operations with a coefficient and a parameter expression.
const CoefExpr operator+( const Coef& left, const ParameterExpr& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const Coef& left, const ParameterExpr& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const Coef& left, const ParameterExpr& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const Coef& left, const ParameterExpr& right )
{
  return CoefExpr( left, right, Operation::div );
}



// Operations with a parameter and a coefficient.
const CoefExpr operator+( const Parameter& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const Parameter& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const Parameter& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const Parameter& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::div );
}


// Operations with a coefficient and a parameter.
const CoefExpr operator+( const CoefExpr& left, const Parameter& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const CoefExpr& left, const Parameter& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const CoefExpr& left, const Parameter& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const CoefExpr& left, const Parameter& right )
{
  return CoefExpr( left, right, Operation::div );
}



// Operations with a parameter expression and a coefficient.
const CoefExpr operator+( const ParameterExpr& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const ParameterExpr& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const ParameterExpr& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const ParameterExpr& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::div );
}


// Operations with a coefficient and a parameter expression.
const CoefExpr operator+( const CoefExpr& left, const ParameterExpr& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const CoefExpr& left, const ParameterExpr& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const CoefExpr& left, const ParameterExpr& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const CoefExpr& left, const ParameterExpr& right )
{
  return CoefExpr( left, right, Operation::div );
}



// Operations with a complex constant and a coefficient.
const CoefExpr operator+( const std::complex< double >& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const std::complex< double >& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const std::complex< double >& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const std::complex< double >& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::div );
}


// Operations with a coefficient and a complex constant.
const CoefExpr operator+( const Coef& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const Coef& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const Coef& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const Coef& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::div );
}



// Operations with a complex constant and a coefficient.
const CoefExpr operator+( const std::complex< double >& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const std::complex< double >& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const std::complex< double >& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const std::complex< double >& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::div );
}


// Operations with a coefficient and a constant.
const CoefExpr operator+( const CoefExpr& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const CoefExpr& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const CoefExpr& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const CoefExpr& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::div );
}



// Operations with a complex constant and a coefficient.
const CoefExpr operator+( const std::complex< double >& left, const Parameter& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const std::complex< double >& left, const Parameter& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const std::complex< double >& left, const Parameter& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const std::complex< double >& left, const Parameter& right )
{
  return CoefExpr( left, right, Operation::div );
}


// Operations with a coefficient and a complex constant.
const CoefExpr operator+( const Parameter& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const Parameter& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const Parameter& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const Parameter& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::div );
}



// Operations with a complex constant and a parameter expression.
const CoefExpr operator+( const std::complex< double >& left, const ParameterExpr& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const std::complex< double >& left, const ParameterExpr& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const std::complex< double >& left, const ParameterExpr& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const std::complex< double >& left, const ParameterExpr& right )
{
  return CoefExpr( left, right, Operation::div );
}


// Operations with a parameter expression and a constant.
const CoefExpr operator+( const ParameterExpr& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const ParameterExpr& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const ParameterExpr& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const ParameterExpr& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::div );
}



// Operations with two coefficient expressions.
const CoefExpr operator+( const CoefExpr& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const CoefExpr& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const CoefExpr& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const CoefExpr& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::div );
}


// Operations with a constant and a coefficient expression.
const CoefExpr operator+( const double& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const double& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const double& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const double& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::div );
}


// Operations with a coefficient expression and a constant.
const CoefExpr operator+( const CoefExpr& left, const double& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const CoefExpr& left, const double& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const CoefExpr& left, const double& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const CoefExpr& left, const double& right )
{
  return CoefExpr( left, right, Operation::div );
}


// Operations with a coefficient and a coefficient expression.
const CoefExpr operator+( const Coef& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const Coef& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const Coef& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const Coef& left, const CoefExpr& right )
{
  return CoefExpr( left, right, Operation::div );
}


// Operations with a coefficient expression and a coefficient.
const CoefExpr operator+( const CoefExpr& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::plus );
}

const CoefExpr operator-( const CoefExpr& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::minus );
}

const CoefExpr operator*( const CoefExpr& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::mult );
}

const CoefExpr operator/( const CoefExpr& left, const Coef& right )
{
  return CoefExpr( left, right, Operation::div );
}



// Unary minus.
const CoefExpr operator-( const Coef& coef )
{
  return CoefExpr( coef, Operation::minus );
}

const CoefExpr operator-( const CoefExpr& coef )
{
  return CoefExpr( coef, Operation::minus );
}


// Binary operations with a coefficient.
const CoefExpr pow( const Coef& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::pow );
}

const CoefExpr pow( const CoefExpr& left, const std::complex< double >& right )
{
  return CoefExpr( left, right, Operation::pow );
}


// Unary operations with a coefficient.
const CoefExpr exp( const Coef& coef )
{
  return CoefExpr( coef, Operation::exp );
}

const CoefExpr log( const Coef& coef )
{
  return CoefExpr( coef, Operation::log );
}

const CoefExpr sin( const Coef& coef )
{
  return CoefExpr( coef, Operation::sin );
}

const CoefExpr cos( const Coef& coef )
{
  return CoefExpr( coef, Operation::cos );
}

const CoefExpr tan( const Coef& coef )
{
  return CoefExpr( coef, Operation::tan );
}

const CoefExpr tanh( const Coef& coef )
{
  return CoefExpr( coef, Operation::tanh );
}

const CoefExpr atanh( const Coef& coef )
{
  return CoefExpr( coef, Operation::atanh );
}


const CoefExpr exp( const CoefExpr& coef )
{
  return CoefExpr( coef, Operation::exp );
}

const CoefExpr log( const CoefExpr& coef )
{
  return CoefExpr( coef, Operation::log );
}

const CoefExpr sin( const CoefExpr& coef )
{
  return CoefExpr( coef, Operation::sin );
}

const CoefExpr cos( const CoefExpr& coef )
{
  return CoefExpr( coef, Operation::cos );
}

const CoefExpr tan( const CoefExpr& coef )
{
  return CoefExpr( coef, Operation::tan );
}

const CoefExpr tanh( const CoefExpr& coef )
{
  return CoefExpr( coef, Operation::tanh );
}

const CoefExpr atanh( const CoefExpr& coef )
{
  return CoefExpr( coef, Operation::atanh );
}




// Operations with two pdf models.
const PdfExpr operator+( const PdfModel& left, const PdfModel& right )
{
  // Cannot add two pdfs that do not depend on exactly the same variables.
  if ( left.varNames() != right.varNames() )
    throw PdfException( "Cannot add two pdfs that do not depend on the same variables." );

  return PdfExpr( left, right, Operation::plus );
}

const PdfExpr operator*( const PdfModel& left, const PdfModel& right )
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

  return PdfExpr( left, right, Operation::mult );
}

// Operations with a pdf model and a pdf expression.
const PdfExpr operator+( const PdfModel& left, const PdfExpr& right )
{
  // Cannot add two pdfs that do not depend on exactly the same variables.
  if ( left.varNames() != right.varNames() )
    throw PdfException( "Cannot add two pdfs that do not depend on the same variables." );

  return PdfExpr( left, right, Operation::plus );
}

const PdfExpr operator*( const PdfModel& left, const PdfExpr& right )
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

  return PdfExpr( left, right, Operation::mult );
}

// Operations with a pdf expression and a pdf model.
const PdfExpr operator+( const PdfExpr& left, const PdfModel& right )
{
  // Cannot add two pdfs that do not depend on exactly the same variables.
  if ( left.varNames() != right.varNames() )
    throw PdfException( "Cannot add two pdfs that do not depend on the same variables." );

  return PdfExpr( left, right, Operation::plus );
}

const PdfExpr operator*( const PdfExpr& left, const PdfModel& right )
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

  return PdfExpr( left, right, Operation::mult );
}

// Operations with two pdf expressions.
const PdfExpr operator+( const PdfExpr& left, const PdfExpr& right )
{
  // Cannot add two pdfs that do not depend on exactly the same variables.
  if ( left.varNames() != right.varNames() )
    throw PdfException( "Cannot add two pdfs that do not depend on the same variables." );

  return PdfExpr( left, right, Operation::plus );
}

const PdfExpr operator*( const PdfExpr& left, const PdfExpr& right )
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

  return PdfExpr( left, right, Operation::mult );
}


// Operations with a parameter and a model.
const PdfExpr operator*( const Parameter& left, const PdfModel& right )
{
  return PdfExpr( left, right, Operation::mult );
}

const PdfExpr operator*( const PdfModel& left, const Parameter& right )
{
  return PdfExpr( left, right, Operation::mult );
}

const PdfExpr operator/( const PdfModel& left, const Parameter& right )
{
  return PdfExpr( left, right, Operation::div );
}


// Operations with a parameter expression and a model.
const PdfExpr operator*( const ParameterExpr& left, const PdfModel& right )
{
  return PdfExpr( left, right, Operation::mult );
}

const PdfExpr operator*( const PdfModel& left, const ParameterExpr& right )
{
  return PdfExpr( left, right, Operation::mult );
}

const PdfExpr operator/( const PdfModel& left, const ParameterExpr& right )
{
  return PdfExpr( left, right, Operation::div );
}


// Operations with a constant and a model.
const PdfExpr operator*( const double& left, const PdfModel& right )
{
  return PdfExpr( left, right, Operation::mult );
}

const PdfExpr operator*( const PdfModel& left, const double& right )
{
  return PdfExpr( left, right, Operation::mult );
}

const PdfExpr operator/( const PdfModel& left, const double& right )
{
  return PdfExpr( left, right, Operation::div );
}


// Operations with a parameter and a pdf expression.
const PdfExpr operator*( const Parameter& left, const PdfExpr& right )
{
  return PdfExpr( left, right, Operation::mult );
}

const PdfExpr operator*( const PdfExpr& left, const Parameter& right )
{
  return PdfExpr( left, right, Operation::mult );
}

const PdfExpr operator/( const PdfExpr& left, const Parameter& right )
{
  return PdfExpr( left, right, Operation::div );
}


// Operations with a parameter expression and a pdf expression.
const PdfExpr operator*( const ParameterExpr& left, const PdfExpr& right )
{
  return PdfExpr( left, right, Operation::mult );
}

const PdfExpr operator*( const PdfExpr& left, const ParameterExpr& right )
{
  return PdfExpr( left, right, Operation::mult );
}

const PdfExpr operator/( const PdfExpr& left, const ParameterExpr& right )
{
  return PdfExpr( left, right, Operation::div );
}


// Operations with a constant and a pdf expression.
const PdfExpr operator*( const double& left, const PdfExpr& right )
{
  return PdfExpr( left, right, Operation::mult );
}

const PdfExpr operator*( const PdfExpr& left, const double& right )
{
  return PdfExpr( left, right, Operation::mult );
}

const PdfExpr operator/( const PdfExpr& left, const double& right )
{
  return PdfExpr( left, right, Operation::div );
}

