
#include <sstream>
#include <vector>
#include <stack>
#include <cmath>
#include <random>

#include <Minuit/FunctionMinimum.h>
#include <Minuit/MnUserParameters.h>
#include <Minuit/MinuitParameter.h>

#include <cfit/parameter.hh>
#include <cfit/parameterexpr.hh>
#include <cfit/pdfmodel.hh>
#include <cfit/pdfexpr.hh>
#include <cfit/operation.hh>

#include <cfit/random.hh>



PdfExpr::PdfExpr( const PdfExpr& right )
{
  _varMap = right._varMap;
  _parMap = right._parMap;

  _expression  = right._expression;
  _opers       = right._opers;
  _ctnts       = right._ctnts;
  _parms       = right._parms;

  std::transform( right._pdfs.begin(), right._pdfs.end(),
                  std::back_inserter( _pdfs ), std::mem_fun( &PdfModel::copy ) );

  _limits = right._limits;
  _scale  = right._scale;
}


PdfExpr::~PdfExpr()
{
  // Delete all the pointers to pdf models, since they have been allocated
  //    by each copy() function of each pdf model.
  typedef std::vector< PdfModel* >::iterator mIter;
  for ( mIter model = _pdfs.begin(); model != _pdfs.end(); ++model )
    delete *model;
}


void PdfExpr::clear()
{
  _varMap.clear();
  _parMap.clear();
  _opers .clear();
  _ctnts .clear();
  _parms .clear();

  // Delete all the pointers to pdf models, since they have been allocated
  //    by each copy() function of each pdf model.
  typedef std::vector< PdfModel* >::iterator mIter;
  for ( mIter model = _pdfs.begin(); model != _pdfs.end(); ++model )
    delete *model;

  _pdfs.clear();

  _expression.clear();
}


// Append a model.
void PdfExpr::append( const PdfModel& model )
{
  _varMap.insert( model._varMap.begin(), model._varMap.end() );
  _parMap.insert( model._parMap.begin(), model._parMap.end() );
  _pdfs.push_back( model.copy() );

  _expression += "m"; // m = model.
}


// Append a pdf expression.
void PdfExpr::append( const PdfExpr& pdf )
{
  _varMap.insert(               pdf._varMap.begin(), pdf._varMap .end() );
  _parMap.insert(               pdf._parMap.begin(), pdf._parMap .end() );
  _opers .insert( _opers.end(), pdf._opers .begin(), pdf._opers  .end() );
  _ctnts .insert( _ctnts.end(), pdf._ctnts .begin(), pdf._ctnts  .end() );
  _parms .insert( _parms.end(), pdf._parms .begin(), pdf._parms  .end() );

  std::transform( pdf._pdfs.begin(), pdf._pdfs.end(),
                  std::back_inserter( _pdfs ), std::mem_fun( &PdfModel::copy ) );

  _expression += pdf._expression;
}

// Append a parameter.
void PdfExpr::append( const Parameter& par )
{
  _parMap[ par.name() ] = par;
  _parms.push_back( par );

  _expression += "p"; // p = parameter.
}

// Append a parameter expression.
void PdfExpr::append( const ParameterExpr& expr )
{
  for ( std::vector< Parameter >::const_iterator par = expr._parms.begin(); par != expr._parms.end(); ++par )
    _parMap[ par->name() ] = *par;
  _opers.insert( _opers.end(), expr._opers.begin(), expr._opers.end() );
  _ctnts.insert( _ctnts.end(), expr._ctnts.begin(), expr._ctnts.end() );
  _parms.insert( _parms.end(), expr._parms.begin(), expr._parms.end() );

  _expression += expr._expression;
}

// Append a constant.
void PdfExpr::append( const double& ctnt )
{
  _ctnts.push_back( ctnt );

  _expression += "c"; // c = constant.
}

// Append a binary operation. No unary operation should ever be appended.
void PdfExpr::append( const Operation::Op& oper )
{
  _opers.push_back( oper );

  _expression += "b";
}


// Assignment operations.
const PdfExpr& PdfExpr::operator=( const PdfModel& right )
{
  clear();

  append( right );
  return *this;
}


const PdfExpr& PdfExpr::operator=( const PdfExpr& right )
{
  clear();

  append( right );
  return *this;
}


// Assignment operations with a pdf model.
const PdfExpr& PdfExpr::operator+=( const PdfModel& right ) throw( PdfException )
{
  // Cannot add two pdfs that do not depend on exactly the same variables.
  if ( ! this->_expression.empty() )
    if ( this->varNames() != right.varNames() )
      throw PdfException( "Cannot add two pdfs that do not depend on the same variables." );

  // If this is the first element to be appended to the expression,
  //    it does not make sense to operate with previous elements.
  if ( this->_expression.empty() )
  {
    append( right );
    return *this;
  }

  append( right           );
  append( Operation::plus );
  return *this;
}

const PdfExpr& PdfExpr::operator*=( const PdfModel& right ) throw( PdfException )
{
  // Cannot multiply two pdfs that share some variable.
  std::vector< std::string > theseVars = this->varNames();
  std::vector< std::string > rightVars = right.varNames();
  std::vector< std::string > intersect;
  std::set_intersection( theseVars.begin(), theseVars.end(),
			 rightVars.begin(), rightVars.end(),
			 std::back_inserter( intersect ) );

  if ( ! intersect.empty() )
    throw PdfException( "Cannot multiply two pdfs that depend on some common variable." );

  append( right           );
  append( Operation::mult );
  return *this;
}

// Assignment operations with a pdf expression.
const PdfExpr& PdfExpr::operator+=( const PdfExpr& right ) throw( PdfException )
{
  // Cannot add two pdfs that do not depend on exactly the same variables.
  if ( ! _expression.empty() )
    if ( this->varNames() != right.varNames() )
      throw PdfException( "Cannot add two pdfs that do not depend on the same variables." );

  // If this is the first element to be appended to the expression,
  //    it does not make sense to operate with previous elements.
  if ( this->_expression.empty() )
  {
    append( right );
    return *this;
  }

  append( right           );
  append( Operation::plus );
  return *this;
}

const PdfExpr& PdfExpr::operator*=( const PdfExpr& right ) throw( PdfException )
{
  // Cannot multiply two pdfs that share some variable.
  std::vector< std::string > theseVars = this->varNames();
  std::vector< std::string > rightVars = right.varNames();
  std::vector< std::string > intersect;
  std::set_intersection( theseVars.begin(), theseVars.end(),
			 rightVars.begin(), rightVars.end(),
			 std::back_inserter( intersect ) );

  if ( ! intersect.empty() )
    throw PdfException( "Cannot multiply two pdfs that depend on some common variable." );

  append( right           );
  append( Operation::mult );
  return *this;
}

// Assignment operations with a parameter.
const PdfExpr& PdfExpr::operator*=( const Parameter& right )
{
  append( right           );
  append( Operation::mult );
  return *this;
}

const PdfExpr& PdfExpr::operator/=( const Parameter& right )
{
  append( right          );
  append( Operation::div );
  return *this;
}

// Assignment operations with a parameter expression.
const PdfExpr& PdfExpr::operator*=( const ParameterExpr& right )
{
  append( right           );
  append( Operation::mult );
  return *this;
}

const PdfExpr& PdfExpr::operator/=( const ParameterExpr& right )
{
  append( right          );
  append( Operation::div );
  return *this;
}

// Assignment operations with a constant.
const PdfExpr& PdfExpr::operator*=( const double& right )
{
  append( right           );
  append( Operation::mult );
  return *this;
}

const PdfExpr& PdfExpr::operator/=( const double& right )
{
  append( right          );
  append( Operation::div );
  return *this;
}


// Setter for individual parameter.
void PdfExpr::setPar( const std::string& name, const double& val, const double& err ) throw( PdfException )
{
  if ( ! _parMap.count( name ) )
    throw PdfException( "Cannot set unexisting parameter " + name + "." );

  _parMap[ name ].set( val, err );

  // Propagate the values to the list of pdfs.
  typedef std::vector< PdfModel* >::const_iterator pdfIter;
  for ( pdfIter pdf = _pdfs.begin(); pdf != _pdfs.end(); ++pdf )
    if ( (*pdf)->_parMap.count( name ) )
      (*pdf)->_parMap[ name ].set( val, err );
}


void PdfExpr::setPars( const std::vector< double >& pars ) throw( PdfException )
{
  if ( _parMap.size() != pars.size() )
    throw PdfException( "PdfExpr::setPars( vector ): Number of arguments passed does not match number of required arguments." );

  // Set the local values of the parameters.
  typedef std::map< std::string, Parameter >::iterator pIter;
  int index = 0;
  for ( pIter par = _parMap.begin(); par != _parMap.end(); ++par )
    par->second.setValue( pars[ index++ ] );

  // Propagate the values to the list of pdfs.
  typedef std::vector< PdfModel* >::const_iterator pdfIter;
  for ( pdfIter pdf = _pdfs.begin(); pdf != _pdfs.end(); ++pdf )
    (*pdf)->setPars( _parMap );
}


void PdfExpr::setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException )
{
  // Set the local values of the parameters.
  typedef std::map< std::string, Parameter >::const_iterator pIter;
  for ( pIter par = pars.begin(); par != pars.end(); ++par )
    _parMap[ par->first ].set( par->second.value(), par->second.error() );

  // Propagate the values to the list of pdfs.
  typedef std::vector< PdfModel* >::const_iterator pdfIter;
  for ( pdfIter pdf = _pdfs.begin(); pdf != _pdfs.end(); ++pdf )
    (*pdf)->setPars( _parMap );
}


void PdfExpr::setPars( const FunctionMinimum& min ) throw( PdfException )
{
  const MnUserParameters& pars = min.userParameters();
  const std::vector< MinuitParameter >& parVec = pars.parameters();

  if ( _parMap.size() != parVec.size() )
    throw PdfException( "PdfExpr::setPars( minimum ): Number of arguments passed does not match number of required arguments." );

  // Set the local values of the parameters.
  typedef std::vector< MinuitParameter >::const_iterator pIter;
  for ( pIter par = parVec.begin(); par != parVec.end(); ++par )
    _parMap[ par->name() ].set( par->value(), par->error() );

  // Propagate the values to the list of pdfs.
  typedef std::vector< PdfModel* >::const_iterator pdfIter;
  for ( pdfIter pdf = _pdfs.begin(); pdf != _pdfs.end(); ++pdf )
    (*pdf)->setPars( _parMap );
}


void PdfExpr::cache()
{
  typedef std::vector< PdfModel* >::const_iterator pIter;
  for ( pIter pdf = _pdfs.begin(); pdf != _pdfs.end(); ++pdf )
    (*pdf)->cache();

  return;
}



const std::map< std::string, double > PdfExpr::generate() const throw( PdfException )
{
  std::map< std::string, double > genVals;
  double val = 0.0;
  double min = 0.0;
  double max = 0.0;

  typedef std::map< std::string, std::pair< double, double > >::const_iterator lIter;

  // Generate full events until one of them is valid.
  bool valid = false;
  while ( ! valid )
  {
    valid = true;
    genVals = generateFull();
    for ( lIter lim = _limits.begin(); valid && ( lim != _limits.end() ); ++lim )
    {
      val = genVals[ lim->first ];
      min = lim->second.first;
      max = lim->second.second;

      if ( ( val < min ) || ( val > max ) )
        valid = false;
    }
  }

  return genVals;
}


// Generate over the full range of definition of all the models. It may happen that
//    some variable gets generated outside of its specified range.
const std::map< std::string, double > PdfExpr::generateFull() const throw( PdfException )
{
  std::stack< double >                                values;
  std::stack< std::pair< std::size_t, std::size_t > > ranges;

  std::pair< std::size_t, std::size_t > xrange;
  std::pair< std::size_t, std::size_t > yrange;

  double x;
  double y;
  std::vector< Parameter     >::const_iterator par = _parms.begin();
  std::vector< double        >::const_iterator ctt = _ctnts.begin();
  std::vector< Operation::Op >::const_iterator ops = _opers.begin();

  std::vector< bool > keep( _pdfs.size(), true );

  // Range of the current block of pdfs.
  std::size_t first = 0;
  std::size_t last  = 0;

  double rnd = 0.0;

  Operation::Op op;

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
  {
    if ( *ch == 'm' )
    {
      values.push( 1.0 );
      ranges.push( std::make_pair( first++, ++last ) );
    }
    else if ( *ch == 'p' )
    {
      values.push( _parMap.find( par++->name() )->second.value() );
      ranges.push( std::make_pair( first, last ) );
    }
    else if ( *ch == 'c' )
    {
      values.push( *ctt++ );
      ranges.push( std::make_pair( first, last ) );
    }
    else
    {
      if ( *ch == 'u' )
      {
        if ( values.empty() )
          throw PdfException( "Parse error: not enough values in the stack." );
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, *ops++ ) );
        xrange = ranges.top();
        if ( xrange.first != xrange.second )
          throw PdfException( "Generation error: trying to apply unary operation to a pdf expression." );
      }
      else if ( *ch == 'b' )
      {
        if ( values.size() < 2 )
          throw PdfException( "Parse error: not enough values in the stack." );

        op = *ops++;

        y = values.top();
        values.pop();
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, y, op ) );

        yrange = ranges.top();
        ranges.pop();
        xrange = ranges.top();
        ranges.pop();
        if ( ( op == Operation::plus ) || ( op == Operation::minus ) )
        {
          if ( xrange == yrange )
            ranges.push( xrange );
          else
          {
            if ( xrange.first == xrange.second )
              throw PdfException( "Generation error: this should have never happened." );

            if ( op == Operation::plus )
              if ( ( x < 0 ) || ( y < 0 ) )
                throw PdfException( "Generation error: negative coefficient multiplying a pdf." );

            if ( op == Operation::minus )
              if ( ( x < 0 ) || ( y > 0 ) )
                throw PdfException( "Generation error: negative coefficient multiplying a pdf." );

            rnd = Random::uniform( 0, x + std::fabs( y ) );

            // Choose one of the two distributions.
            if ( rnd < x )
            {
              ranges.push( xrange );
              std::fill( keep.begin() + yrange.first, keep.begin() + yrange.second, false );
            }
            else
            {
              ranges.push( yrange );
              std::fill( keep.begin() + xrange.first, keep.begin() + xrange.second, false );
            }
          }
        }
        else // Operations not plus or minus.
          ranges.push( std::make_pair( xrange.first, yrange.second ) );
      }
      else
        throw PdfException( std::string( "Parse error: unknown operation " ) + *ch + "." );
    }
  }

  if ( values.size() != 1 )
    throw PdfException( "PdfExpr parse error: too many values have been supplied." );

  std::map< std::string, double > entry;
  std::map< std::string, double > partial;
  for ( unsigned pdf = 0; pdf < _pdfs.size(); ++pdf )
    if ( keep[ pdf ] )
    {
      partial = _pdfs[ pdf ]->generate();
      entry.insert( partial.begin(), partial.end() );
    }

  return entry;
}


const std::map< unsigned, std::vector< double > > PdfExpr::cacheReal( const Dataset& data )
{
  std::map< unsigned, std::vector< double > > cache;

  for ( std::vector< PdfModel* >::iterator pdf = _pdfs.begin(); pdf != _pdfs.end(); ++pdf )
    insert( cache, (*pdf)->cacheReal( data ) );

  return cache;
}



const std::map< unsigned, std::vector< std::complex< double > > > PdfExpr::cacheComplex( const Dataset& data )
{
  std::map< unsigned, std::vector< std::complex< double > > > cache;

  for ( std::vector< PdfModel* >::iterator pdf = _pdfs.begin(); pdf != _pdfs.end(); ++pdf )
    insert( cache, (*pdf)->cacheComplex( data ) );

  return cache;
}



const double PdfExpr::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  if ( _varMap.size() != vars.size() )
    throw PdfException( "PdfExpr::evaluate: Number of arguments passed does not match number of required arguments." );

  // Dictionary of the variable names with the values passed.
  std::map< std::string, double > localVars;
  std::vector< double > modelVars;

  // Set the local values of the variables.
  int index = 0;
  typedef std::map< std::string, Variable >::const_iterator vIter;
  for ( vIter var = _varMap.begin(); var != _varMap.end(); ++var )
    localVars[ var->first ] = vars[ index++ ];

  std::stack< double > values;

  double x;
  double y;
  std::vector< PdfModel*     >::const_iterator pdf = _pdfs .begin();
  std::vector< Parameter     >::const_iterator par = _parms.begin();
  std::vector< double        >::const_iterator ctt = _ctnts.begin();
  std::vector< Operation::Op >::const_iterator ops = _opers.begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
    if ( *ch == 'm' )
    {
      // Determine the variables that the pdf depends on.
      modelVars.clear();
      std::map< std::string, Variable >& pdfVars = (*pdf)->_varMap;
      for ( vIter var = pdfVars.begin(); var != pdfVars.end(); ++var )
        modelVars.push_back( localVars.find( var->second.name() )->second );

      // Evaluate the function at the given point.
      values.push( (*pdf++)->evaluate( modelVars ) );
    }
    else if ( *ch == 'p' )
      values.push( _parMap.find( par++->name() )->second.value() );
    else if ( *ch == 'c' )
      values.push( *ctt++ );
    else
    {
      if ( *ch == 'b' )
      {
        if ( values.size() < 2 )
          throw PdfException( "Parse error: not enough values in the stack." );
        y = values.top();
        values.pop();
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, y, *ops++ ) );
      }
      else if ( *ch == 'u' )
      {
        if ( values.empty() )
          throw PdfException( "Parse error: not enough values in the stack." );
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, *ops++ ) );
      }
      else
        throw PdfException( std::string( "Parse error: unknown operation " ) + *ch + "." );
    }

  if ( values.size() != 1 )
    throw PdfException( "PdfExpr parse error: too many values have been supplied." );

  return values.top();
}



const double PdfExpr::evaluate( const std::vector< double                 >& vars  ,
                                const std::vector< double                 >& cacheR,
                                const std::vector< std::complex< double > >& cacheC  ) const throw( PdfException )
{
  if ( _varMap.size() != vars.size() )
    throw PdfException( "PdfExpr::evaluate: Number of arguments passed does not match number of required arguments." );

  // Dictionary of the variable names with the values passed.
  std::map< std::string, double > localVars;
  std::vector< double > modelVars;

  // Set the local values of the variables.
  int index = 0;
  typedef std::map< std::string, Variable >::const_iterator vIter;
  for ( vIter var = _varMap.begin(); var != _varMap.end(); ++var )
    localVars[ var->first ] = vars[ index++ ];

  std::stack< double > values;

  double x;
  double y;
  std::vector< PdfModel*     >::const_iterator pdf = _pdfs .begin();
  std::vector< Parameter     >::const_iterator par = _parms.begin();
  std::vector< double        >::const_iterator ctt = _ctnts.begin();
  std::vector< Operation::Op >::const_iterator ops = _opers.begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
    if ( *ch == 'm' )
    {
      // Determine the variables that the pdf depends on.
      modelVars.clear();
      std::map< std::string, Variable >& pdfVars = (*pdf)->_varMap;
      for ( vIter var = pdfVars.begin(); var != pdfVars.end(); ++var )
        modelVars.push_back( localVars.find( var->second.name() )->second );

      // Evaluate the function at the given point.
      values.push( (*pdf++)->evaluate( modelVars, cacheR, cacheC ) );
    }
    else if ( *ch == 'p' )
      values.push( _parMap.find( par++->name() )->second.value() );
    else if ( *ch == 'c' )
      values.push( *ctt++ );
    else
    {
      if ( *ch == 'b' )
      {
        if ( values.size() < 2 )
          throw PdfException( "Parse error: not enough values in the stack." );
        y = values.top();
        values.pop();
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, y, *ops++ ) );
      }
      else if ( *ch == 'u' )
      {
        if ( values.empty() )
          throw PdfException( "Parse error: not enough values in the stack." );
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, *ops++ ) );
      }
      else
        throw PdfException( std::string( "Parse error: unknown operation " ) + *ch + "." );
    }

  if ( values.size() != 1 )
    throw PdfException( "PdfExpr parse error: too many values have been supplied." );

  return values.top();
}


void PdfExpr::setLimits( const Variable& var, const double& min, const double& max )
{
  const double& totalArea = this->area();
  const double& limitArea = this->area( var.name(), min, max );

  _scale *= totalArea / limitArea;
  _limits[ var.name() ] = std::make_pair( min, max );
}


// Get a vector with the names of the variables common in all the
//    products of pdfs. This is the list of variables that can be
//    integrated over in a convolution operation.
std::vector< std::string > PdfExpr::commonVars() const throw( PdfException )
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
  std::vector< PdfModel* >::const_iterator models = _pdfs.begin();

  typedef std::string::const_iterator eIter;
  std::vector< std::string > varNames;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
  {
    if ( *ch == 'm' )
    {
      std::map< std::string, Variable >& vars = (*models++)->_varMap;
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
  }

  if ( calcs.size() != 1 )
    throw PdfException( "Parse error computing convolution: too many values have been supplied." );

  return calcs.top();
}



double PdfExpr::area( const std::string& var, const double& min, const double& max ) const throw( PdfException )
{
  std::stack< double > values;

  double x;
  double y;
  std::vector< PdfModel*     >::const_iterator pdf = _pdfs .begin();
  std::vector< Parameter     >::const_iterator par = _parms.begin();
  std::vector< double        >::const_iterator ctt = _ctnts.begin();
  std::vector< Operation::Op >::const_iterator ops = _opers.begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
    if ( *ch == 'm' )
    {
      if ( (*pdf)->dependsOn( var ) )
        values.push( (*pdf++)->area( min, max ) );
      else
      {
        values.push( 1.0 );
        pdf++;
      }
    }
    else if ( *ch == 'p' )
      values.push( _parMap.find( par++->name() )->second.value() );
    else if ( *ch == 'c' )
      values.push( *ctt++ );
    else
    {
      if ( *ch == 'b' )
      {
        if ( values.size() < 2 )
          throw PdfException( "Parse error: not enough values in the stack." );
        y = values.top();
        values.pop();
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, y, *ops++ ) );
      }
      else if ( *ch == 'u' )
      {
        if ( values.empty() )
          throw PdfException( "Parse error: not enough values in the stack." );
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, *ops++ ) );
      }
      else
        throw PdfException( std::string( "Parse error: unknown operation " ) + *ch + "." );
    }

  if ( values.size() != 1 )
    throw PdfException( "PdfExpr parse error: too many values have been supplied." );

  return values.top();
}



const double PdfExpr::yield() const
{
  std::stack< double > values;

  double x;
  double y;
  std::vector< Parameter     >::const_iterator par = _parms.begin();
  std::vector< double        >::const_iterator ctt = _ctnts.begin();
  std::vector< Operation::Op >::const_iterator ops = _opers.begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
    if ( *ch == 'm' )
      values.push( 1.0 );
    else if ( *ch == 'p' )
      values.push( _parMap.find( par++->name() )->second.value() );
    else if ( *ch == 'c' )
      values.push( *ctt++ );
    else
    {
      if ( *ch == 'b' )
      {
        if ( values.size() < 2 )
          throw PdfException( "Parse error: not enough values in the stack." );
        y = values.top();
        values.pop();
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, y, *ops++ ) );
      }
      else if ( *ch == 'u' )
      {
        if ( values.empty() )
          throw PdfException( "Parse error: not enough values in the stack." );
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, *ops++ ) );
      }
      else
        throw PdfException( std::string( "Parse error: unknown operation " ) + *ch + "." );
    }

  if ( values.size() != 1 )
    throw PdfException( "PdfExpr parse error: too many values have been supplied." );

  return values.top();
}



double PdfExpr::area() const throw( PdfException )
{
  std::stack< double > values;

  double x;
  double y;
  std::vector< Parameter     >::const_iterator par = _parms.begin();
  std::vector< double        >::const_iterator ctt = _ctnts.begin();
  std::vector< Operation::Op >::const_iterator ops = _opers.begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
    if ( *ch == 'm' )
      values.push( 1.0 );
    else if ( *ch == 'p' )
      values.push( _parMap.find( par++->name() )->second.value() );
    else if ( *ch == 'c' )
      values.push( *ctt++ );
    else
    {
      if ( *ch == 'b' )
      {
        if ( values.size() < 2 )
          throw PdfException( "Parse error: not enough values in the stack." );
        y = values.top();
        values.pop();
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, y, *ops++ ) );
      }
      else if ( *ch == 'u' )
      {
        if ( values.empty() )
          throw PdfException( "Parse error: not enough values in the stack." );
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, *ops++ ) );
      }
      else
        throw PdfException( std::string( "Parse error: unknown operation " ) + *ch + "." );
    }

  if ( values.size() != 1 )
    throw PdfException( "PdfExpr parse error: too many values have been supplied." );

  return values.top();
}



const double PdfExpr::project( const std::string& varName, const double& value ) const throw( PdfException )
{
  std::stack< double > values;

  double x;
  double y;
  std::vector< PdfModel*     >::const_iterator pdf = _pdfs .begin();
  std::vector< Parameter     >::const_iterator par = _parms.begin();
  std::vector< double        >::const_iterator ctt = _ctnts.begin();
  std::vector< Operation::Op >::const_iterator ops = _opers.begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
    if ( *ch == 'm' )
      // Project the function at the given point.
      values.push( (*pdf++)->project( varName, value ) );
    else if ( *ch == 'p' )
      values.push( _parMap.find( par++->name() )->second.value() );
    else if ( *ch == 'c' )
      values.push( *ctt++ );
    else
    {
      if ( *ch == 'b' )
      {
        if ( values.size() < 2 )
          throw PdfException( "Parse error: not enough values in the stack." );
        y = values.top();
        values.pop();
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, y, *ops++ ) );
      }
      else if ( *ch == 'u' )
      {
        if ( values.empty() )
          throw PdfException( "Parse error: not enough values in the stack." );
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, *ops++ ) );
      }
      else
        throw PdfException( std::string( "Parse error: unknown operation " ) + *ch + "." );
    }

  if ( values.size() != 1 )
    throw PdfException( "PdfExpr parse error: too many values have been supplied." );

  return values.top();
}



const double PdfExpr::project( const std::string& var1, const std::string& var2,
                               const double&      val1, const double&      val2  ) const throw( PdfException )
{
  std::stack< double > values;

  double x;
  double y;
  std::vector< PdfModel*     >::const_iterator pdf = _pdfs .begin();
  std::vector< Parameter     >::const_iterator par = _parms.begin();
  std::vector< double        >::const_iterator ctt = _ctnts.begin();
  std::vector< Operation::Op >::const_iterator ops = _opers.begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
    if ( *ch == 'm' )
      // Project the function at the given point.
      values.push( (*pdf++)->project( var1, var2, val1, val2 ) );
    else if ( *ch == 'p' )
      values.push( _parMap.find( par++->name() )->second.value() );
    else if ( *ch == 'c' )
      values.push( *ctt++ );
    else
    {
      if ( *ch == 'b' )
      {
        if ( values.size() < 2 )
          throw PdfException( "Parse error: not enough values in the stack." );
        y = values.top();
        values.pop();
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, y, *ops++ ) );
      }
      else if ( *ch == 'u' )
      {
        if ( values.empty() )
          throw PdfException( "Parse error: not enough values in the stack." );
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, *ops++ ) );
      }
      else
        throw PdfException( std::string( "Parse error: unknown operation " ) + *ch + "." );
    }

  if ( values.size() != 1 )
    throw PdfException( "PdfExpr parse error: too many values have been supplied." );

  return values.top();
}



const double PdfExpr::project( const std::string& varName,
                               const double&      value  ,
                               const Region&      region   ) const throw( PdfException )
{
  std::stack< double > values;

  double x;
  double y;
  std::vector< PdfModel*     >::const_iterator pdf = _pdfs .begin();
  std::vector< Parameter     >::const_iterator par = _parms.begin();
  std::vector< double        >::const_iterator ctt = _ctnts.begin();
  std::vector< Operation::Op >::const_iterator ops = _opers.begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
    if ( *ch == 'm' )
      // Project the function at the given point.
      values.push( (*pdf++)->project( varName, value, region ) );
    else if ( *ch == 'p' )
      values.push( _parMap.find( par++->name() )->second.value() );
    else if ( *ch == 'c' )
      values.push( *ctt++ );
    else
    {
      if ( *ch == 'b' )
      {
        if ( values.size() < 2 )
          throw PdfException( "Parse error: not enough values in the stack." );
        y = values.top();
        values.pop();
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, y, *ops++ ) );
      }
      else if ( *ch == 'u' )
      {
        if ( values.empty() )
          throw PdfException( "Parse error: not enough values in the stack." );
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, *ops++ ) );
      }
      else
        throw PdfException( std::string( "Parse error: unknown operation " ) + *ch + "." );
    }

  if ( values.size() != 1 )
    throw PdfException( "PdfExpr parse error: too many values have been supplied." );

  return values.top();
}



const double PdfExpr::project( const std::string& var1, const std::string& var2,
                               const double&      val1, const double&      val2,
                               const Region&      region                         ) const throw( PdfException )
{
  std::stack< double > values;

  double x;
  double y;
  std::vector< PdfModel*     >::const_iterator pdf = _pdfs .begin();
  std::vector< Parameter     >::const_iterator par = _parms.begin();
  std::vector< double        >::const_iterator ctt = _ctnts.begin();
  std::vector< Operation::Op >::const_iterator ops = _opers.begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
    if ( *ch == 'm' )
      // Project the function at the given point.
      values.push( (*pdf++)->project( var1, var2, val1, val2, region ) );
    else if ( *ch == 'p' )
      values.push( _parMap.find( par++->name() )->second.value() );
    else if ( *ch == 'c' )
      values.push( *ctt++ );
    else
    {
      if ( *ch == 'b' )
      {
        if ( values.size() < 2 )
          throw PdfException( "Parse error: not enough values in the stack." );
        y = values.top();
        values.pop();
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, y, *ops++ ) );
      }
      else if ( *ch == 'u' )
      {
        if ( values.empty() )
          throw PdfException( "Parse error: not enough values in the stack." );
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, *ops++ ) );
      }
      else
        throw PdfException( std::string( "Parse error: unknown operation " ) + *ch + "." );
    }

  if ( values.size() != 1 )
    throw PdfException( "PdfExpr parse error: too many values have been supplied." );

  return values.top();
}
