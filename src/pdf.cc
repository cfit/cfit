
#include <iostream> // TEMPORAL
#include <iterator> // TEMPORAL

#include <sstream>
#include <vector>
#include <stack>

#include <cfit/parameter.hh>
#include <cfit/pdfmodel.hh>
#include <cfit/pdf.hh>



void Pdf::setVars( const std::vector< double >& vars ) throw( PdfException )
{
  if ( _vars.size() != vars.size() )
    throw PdfException( "Number of arguments passed does not match number of required arguments." );

  // Set the local values of the variables.
  typedef std::map< std::string, Variable >::iterator vIter;
  int index = 0;
  for ( vIter var = _vars.begin(); var != _vars.end(); ++var )
    var->second.setValue( vars[ index++ ] );

  // Propagate the values to the list of pdfs.
  typedef std::vector< PdfModel* >::const_iterator pdfIter;
  for ( pdfIter pdf = _pdfVect.begin(); pdf != _pdfVect.end(); ++pdf )
    {
      std::map< std::string, Variable >& pdfVars = (*pdf)->_vars;
      for ( vIter var = pdfVars.begin(); var != pdfVars.end(); var++ )
	var->second.setValue( _vars[ var->second.name() ].value() );
    }
}


void Pdf::setPars( const std::vector< double >& pars ) throw( PdfException )
{
  if ( _pars.size() != pars.size() )
    throw PdfException( "Number of arguments passed does not match number of required arguments." );

  // Set the local values of the parameters.
  typedef std::map< std::string, Parameter >::iterator pIter;
  int index = 0;
  for ( pIter par = _pars.begin(); par != _pars.end(); par++ )
    par->second.setValue( pars[ index++ ] );

  // Propagate the values to the list of pdfs.
  typedef std::vector< PdfModel* >::const_iterator pdfIter;
  for ( pdfIter pdf = _pdfVect.begin(); pdf != _pdfVect.end(); pdf++ )
    {
      std::map< std::string, Parameter >& pdfPars = (*pdf)->_pars;
      for ( pIter par = pdfPars.begin(); par != pdfPars.end(); par++ )
	par->second.setValue( _pars[ par->second.name() ].value() );
    }
}


void Pdf::cache()
{
  typedef std::vector< PdfModel* >::const_iterator pIter;
  for ( pIter pdf = _pdfVect.begin(); pdf != _pdfVect.end(); pdf++ )
    (*pdf)->cache();

  return;
}


// Before running this function, the Pdf::setVars( vars ) function must be called.
//    To avoid the risk of forgetting it, run Pdf::evaluate( vars ).
double Pdf::evaluate() const throw( PdfException )
{
  std::stack< double > values;

  double x;
  double y;
  std::vector< PdfModel*   >::const_iterator pdf = _pdfVect.begin();
  std::vector< std::string >::const_iterator par = _parVect.begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ch++ )
    if ( *ch == 'm' )
      values.push( (*pdf++)->evaluate() );
    else if ( *ch == 'p' )
      values.push( _pars.find( *par++ )->second.value() );
    else
      if ( values.size() < 2 )
	throw PdfException( "Parse error: not enough values in the stack." );
      else
	{
	  x = values.top();
	  values.pop();
	  y = values.top();
	  values.pop();
	  if ( *ch == '+' )
	    values.push( x + y );
	  else if ( *ch == '*' )
	    values.push( x * y );
	  else
	    throw PdfException( std::string( "Parse error: unknown operation " ) + *ch + "." );
	}

  if ( values.size() == 1 )
    return values.top();
  else
    throw PdfException( "Parse error: too many values have been supplied." );

  return 0.;
}


double Pdf::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  if ( _vars.size() != vars.size() )
    throw PdfException( "Number of arguments passed does not match number of required arguments." );

  // Dictionary of the variable names with the values passed.
  std::map< std::string, double > localVars;
  std::vector< double > modelVars;

  // Set the local values of the variables.
  int index = 0;
  typedef std::map< std::string, Variable >::const_iterator vIter;
  for ( vIter var = _vars.begin(); var != _vars.end(); var++ )
    localVars[ var->first ] = vars[ index++ ];

  std::stack< double > values;

  double x;
  double y;
  std::vector< PdfModel*   >::const_iterator pdf = _pdfVect.begin();
  std::vector< std::string >::const_iterator par = _parVect.begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ch++ )
    if ( *ch == 'm' )
      {
	// Determine the variables that the pdf depends on.
	modelVars.clear();
	std::map< std::string, Variable >& pdfVars = (*pdf)->_vars;
	for ( vIter var = pdfVars.begin(); var != pdfVars.end(); ++var )
	  modelVars.push_back( localVars.find( var->second.name() )->second );

	// Evaluate the function at the given point.
	values.push( (*pdf++)->evaluate( modelVars ) );
      }
    else if ( *ch == 'p' )
      values.push( _pars.find( *par++ )->second.value() );
    else
      if ( values.size() < 2 )
	throw PdfException( "Parse error: not enough values in the stack." );
      else
	{
	  x = values.top();
	  values.pop();
	  y = values.top();
	  values.pop();
	  if ( *ch == '+' )
	    values.push( x + y );
	  else if ( *ch == '*' )
	    values.push( x * y );
	  else
	    throw PdfException( std::string( "Parse error: unknown operation " ) + *ch + "." );
	}

  if ( values.size() == 1 )
    return values.top();
  else
    throw PdfException( "Parse error: too many values have been supplied." );

  return 0.;
}



// Get a vector with the names of the variables common in all the
//    products of pdfs. This is the list of variables that can be
//    integrated over in a convolution operation.
std::vector< std::string > Pdf::commonVars() const throw( PdfException )
{
  // Vectors that will contain the variables each model depends on.
  std::vector< std::string > x;
  std::vector< std::string > y;
  std::vector< std::string > z;

  // Empty vector.
  std::vector< std::string > voidVec;

  // Stack for partial calculations.
  std::stack< std::vector< std::string > > calcs;

  // Iterator over the models that the pdf depends on.
  std::vector< PdfModel* >::const_iterator models = _pdfVect.begin();

  typedef std::string::const_iterator eIter;
  std::vector< std::string > varNames;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ch++ )
    if ( *ch == 'm' )
      {
	std::map< std::string, Variable >& vars = (*models++)->_vars;
	std::transform( vars.begin(), vars.end(), std::back_inserter( varNames ), Select1st() );
	calcs.push( varNames );
      }
    else if ( *ch == 'p' )
      calcs.push( voidVec );
    else
      if ( calcs.size() < 2 )
	throw PdfException( "Parse error computing convolution: not enough values in the stack." );
      else
	{
	  x = calcs.top();
	  calcs.pop();
	  y = calcs.top();
	  calcs.pop();
	  z.clear();

	  if ( *ch == '+' )
	    std::set_intersection( x.begin(), x.end(), y.begin(), y.end(), std::back_inserter( z ) );
	  else if ( *ch == '*' )
	    std::set_union       ( x.begin(), x.end(), y.begin(), y.end(), std::back_inserter( z ) );

	  calcs.push( z );
	}

  if ( calcs.size() == 1 )
    return calcs.top();
  else
    throw PdfException( "Parse error computing convolution: too many values have been supplied." );

  return voidVec;
}

