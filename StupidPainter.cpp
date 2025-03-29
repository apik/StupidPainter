/*-------------------------------------------------------------------

  Usage:   ./StupidPainter <NC>

  Syntax:

  [T(1)T(1)]=tr[T^a1_ij*T^a1_ij]

  f(1,2,3)f(1,2,3)=f^{abc}f^{abc}

  --------------------------------------------------------------------*/

#include <regex>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <list>
#include <stdexcept>
#include <set>
#include <map>
#include <tuple>
#include <iomanip>
// GNU readline for prompt editing
#include <readline/readline.h>
#include <readline/history.h>
// color prompt
#define KNRM  "\001\x1B[0m\002"
#define KRED  "\001\x1B[31m\002"

using Eigen::MatrixXcd;


std::list<size_t> int2base(size_t input, size_t base)
{
  std::list<size_t> result;
  while(input) {
    result.push_front(input%base);
    input /= base;
  }
  return result;
}

std::list<size_t> fixLenSeq(size_t input, size_t base, size_t len)
{
  std::list<size_t> res(len,0);

  std::list<size_t> ib = int2base(input,base);

  ib.insert (ib.begin(),len - ib.size(),0);
  return ib;
}

class Timer
{
  std::clock_t    start;
public:
  Timer():start(std::clock())
  {
  }

  double elapsed()
  {
    return (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
  }
};


bool float2rat(long double x, long int& num, long int& den,
               size_t maxit = 10, long double eps = 1e-9)
{
  std::vector<long int> p(maxit, 0);
  std::vector<long int> q(maxit, 0);
  std::vector<long int> a(maxit, 0);
  long int len;
  int i;
  //The first two convergents are 0/1 and 1/0
  p[0] = 0; q[0] = 1;
  p[1] = 1; q[1] = 0;
  //The rest of the convergents (and continued fraction)
  for(size_t i=2; i < maxit; ++i)
    {
      a[i] = lrint(floor(x));
      p[i] = a[i]*p[i-1] + p[i-2];
      q[i] = a[i]*q[i-1] + q[i-2];
      // std::cout << a[i] << ":   " << p[i] << "/" << q[i] << std::endl;
      len = i;
      if(fabs(x-a[i]) < eps)
        {
          num = p[i];
          den = q[i];
          return true;
        }
      x = 1.0/(x - a[i]);
    }
  return false;
}

typedef std::list<size_t> AdjVec;
typedef std::vector<size_t> abc;

typedef std::list<std::list<size_t> > TrVec;
typedef std::list<abc> fVec;

class Color
{
  const double TF = 0.5;
  size_t nc;                    // Fundamental representation
  size_t na;                    // Adjoint representation

  std::vector<MatrixXcd> vgelm; // Gell-Mann matrices

  std::map<std::tuple<size_t,size_t,size_t>, double> fabcMap; // Map of sorted f^abc

  TrVec traces;
  fVec structConsts;

  size_t nIndex;

  size_t seedCntr;
  std::map<size_t,size_t> amap;
public:
  Color(size_t NC): nc(NC), na(NC*NC-1)
  {
    // Constructing basis of Generalized Gell-Mann matrices
    // we use definition from the mathworld:
    // http://mathworld.wolfram.com/GeneralizedGell-MannMatrix.html

    std::cout << "Constructing T^a   >>> ";

    int stepT = std::max(1, (int)ceil(na/10.));
    size_t cntrT = 0;
    Timer  timeT;

    // Symmetric set 1
    for (size_t j = 0; j < nc; j++)
      for (size_t k = j + 1; k < nc; k++)
        {
          MatrixXcd m = MatrixXcd::Constant(nc, nc, 0);
          m(j,k) += 1;
          m(k,j) += 1;

          // vGelm = T
          vgelm.push_back(m/2.);
          cntrT++;
          if((cntrT % stepT) == 0) std::cout <<" ." << std::flush;
        }

    // Antisymmetric set 2
    for (size_t j = 0; j < nc; j++)
      for (size_t k = j + 1; k < nc; k++)
        {
          MatrixXcd m = MatrixXcd::Constant(nc,nc,0);
          m(j,k) -= std::complex<double>(0,1);
          m(k,j) += std::complex<double>(0,1);

          // vGelm = T
          vgelm.push_back(m/2.);
          cntrT++;
          if((cntrT % stepT) == 0) std::cout <<" ." << std::flush;
        }

    // Diagonal set 3
    for (size_t l = 1; l < nc; l++)
      {
        MatrixXcd m = MatrixXcd::Constant(nc,nc,0);
        double coef = sqrt(2./double(l)/double(l+1));

        for (size_t j = 1; j <= l; j++)
          m(j-1,j-1) += coef;

        m(l,l) -= coef * l;

        // vGelm = T
        vgelm.push_back(m/2.);
        cntrT++;
        if((cntrT % stepT) == 0) std::cout <<" ." << std::flush;
      }
    std::cout << " done in " << std::setw(12) << timeT.elapsed() << " ms" << std::endl;

    // Structure constants initialization

    std::cout << "Constructing f^abc >>> ";
    size_t cntrf = 1;
    Timer  timef;
    double nfabc = na*(na - 1)*(na - 2)/6.;
    int stepf = std::max(1, (int)floor(nfabc/10.));

    for (size_t c = 0; c < na; c++)
      for (size_t b = 0; b < c; b++)
        for (size_t a = 0; a < b; a++)
          {
            double fabc = fCalc(a,b,c);

            // We add only non-zero elements
            if(fabc != 0) fabcMap[std::make_tuple(a,b,c)] = fabc;

            if((cntrf % stepf) == 0) std::cout <<" ." << std::flush;
            cntrf++;
          }
    std::cout << " " << fabcMap.size() << "/" << nfabc << " done in " << std::setw(12) << timef.elapsed() << " ms" << std::endl;
  }


  bool nextSeed()
  {
    if(seedCntr < pow(na,amap.size()) - 1)
      {
        seedCntr++;

        std::list<size_t> result = fixLenSeq(seedCntr,na,amap.size());

        std::map<size_t,size_t>::iterator begin1 = amap.begin();
        std::list<size_t>::iterator begin2 = result.begin();
        std::map<size_t,size_t>::iterator end1 = amap.end();
        std::list<size_t>::iterator end2 = result.end();
        std::map<size_t,size_t>::iterator i1;
        std::list<size_t>::iterator i2;
        for (i1 = begin1, i2 = begin2; (i1 != end1) && (i2 != end2); ++i1, ++i2)
          amap[i1->first] = *i2;

        return true;
      }
    return false;
  }

  std::map<size_t,size_t> seed()
  {
    return amap;
  }

  bool parse(std::string str)
  {
    // Erase white spaces and all possible delimiters from the input
    str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
    str.erase(std::remove(str.begin(), str.end(), ';'), str.end());
    str.erase(std::remove(str.begin(), str.end(), ':'), str.end());
    str.erase(std::remove(str.begin(), str.end(), '*'), str.end());

    // We allow for input only
    // "[","]"
    // "(",")"
    // "f","T"
    // ","
    // integer numbers
    // spaces
    std::regex instrRegex("[[:digit:]()\\[\\]Tf,]+");
    if(!std::regex_match (str, instrRegex))
      return false;

    std::cout << std::endl;
    std::cout << "input: " << std::setw(72) <<  str << std::endl;

    // Clear state
    amap.clear();
    traces.clear();
    structConsts.clear();

    // Find traces
    std::size_t found_left = str.find_first_of("[");
    std::size_t found_right = found_left;

    std::regex traceRegex("(?:T\\(\\d+\\)){2,}");
    std::regex tRegex("T\\((\\d+)\\)");

    while (found_left!=std::string::npos)
      {
        found_right = str.find_first_of("]",found_left+1);

        if (found_right!=std::string::npos)
          {
            std::string traceStr = str.substr(found_left+1,found_right - found_left - 1);

            // Parse T^a_ij only if they combined in proper trace string
            if(std::regex_match (traceStr,traceRegex))
              {

                std::regex_token_iterator<std::string::iterator> rend;
                std::regex_token_iterator<std::string::iterator> a ( traceStr.begin(), traceStr.end(), tRegex,1 );
                AdjVec newTrace;
                while (a!=rend)
                  {
                    std::istringstream buffer(*a++);
                    size_t value;
                    buffer >> value; 
                    newTrace.push_back(value);
                  }

                traces.push_back(newTrace);
              }
            else
              throw std::runtime_error("Could not parse trace");
          }
        else
          throw std::runtime_error("Trace is not closed");

        found_left = found_right+1;
        found_left = str.find_first_of("[",found_right+1);
      }


    // Find f^{abc}
    found_left = str.find_first_of("f");
    found_right = found_left;

    std::regex traceFabc("\\((\\d+),(\\d+),(\\d+)\\)");

    while (found_left != std::string::npos)
      {
        found_right = str.find_first_of(")",found_left+1);

        if (found_right != std::string::npos)
          {
            std::string fStr = str.substr(found_left+1,found_right - found_left );

            // Parse f_a,b,c only if they combined in proper trace
            // string
            std::smatch sm;
            if(std::regex_match (fStr,sm,traceFabc))
              {
                abc newFabc;
                for (unsigned i=1; i<sm.size(); ++i)
                  {
                    std::istringstream buffer(sm[i]);
                    size_t value;
                    buffer >> value;
                    newFabc.push_back(value);
                  }

                structConsts.push_back(newFabc);
              }
            else
              throw std::runtime_error("Could not parse f(a,b,c)");
          }
        else
          throw std::runtime_error("f(a,b,c) is not closed");

        found_left = found_right+1;
        found_left = str.find_first_of("f",found_right+1);
      }

    std::multiset<size_t> adjIndices;
    // Add indices from traces
    for(TrVec::const_iterator it = traces.begin(); it != traces.end(); ++it)
      adjIndices.insert(it->begin(),it->end());

    // Add indices from structure constsnts
    for(fVec::const_iterator it = structConsts.begin(); it != structConsts.end(); ++it)
      adjIndices.insert(it->begin(),it->end());

    for(std::multiset<size_t>::const_iterator it = adjIndices.begin(); it != adjIndices.end(); ++it)
      if(adjIndices.count(*it) != 2)
        {
          std::stringstream ss;
          ss << "Index " << *it << " is not paired";
          throw std::runtime_error(ss.str());
        }
      else
        amap[*it] = 0;
    // Map filled with zeros and
    // initial value for counter set
    seedCntr = 0;

    nIndex = adjIndices.size()/2;

    std::cout << "[T^a...] traces:" << std::setw(9) << traces.size() << "  ";
    std::cout << "f^abc constants:" << std::setw(10) << structConsts.size() << "  ";
    std::cout << "gluon lines:"     << std::setw(13) << nIndex << std::endl;

    return true;
  } // parse()


  std::complex<double> contract()
  {
    std::complex<double> sumres(0,0);

    do
      {

        std::complex<double> term(1,0);

        // First multiply by f_abc
        for(fVec::const_iterator it = structConsts.begin(); it != structConsts.end(); ++it)
          term *= f(amap[it->at(0)], amap[it->at(1)], amap[it->at(2)]);

        // Second multiply by traces
        for(TrVec::const_iterator it = traces.begin(); it != traces.end(); ++it)
          {
            std::complex<double> c(trMap(*it));
            term *= c;

            // TODO if zero stop loop
            if(std::norm(c) == 0) break;
          }

        sumres += term;
      }
    while(nextSeed());

    return sumres;
  }

  inline size_t NC()
  {
    return nc;
  }
  inline size_t NA()
  {
    return na;
  }

  const MatrixXcd& t(size_t a) const
  {
    if(a >= na)                 // numbering starts from 0
      throw std::runtime_error("Adjoint index to big");
    return vgelm[a];
  }

  double f(size_t a, size_t b, size_t c)
  {
    if((a == b) || (b == c) || (c == a))
      return 0;

    std::map<std::tuple<size_t,size_t,size_t>, double>::const_iterator fabcit;
    int permSign = 1;

    if((a < b) && (b < c))
      {
        fabcit = fabcMap.find(std::make_tuple(a,b,c));
        return  (fabcit != fabcMap.end()) ? fabcit->second : 0; 
      }
    else if((a < c) && (c < b))
      {
        fabcit = fabcMap.find(std::make_tuple(a,c,b));
        return  (fabcit != fabcMap.end()) ? -fabcit->second : 0; 
      }
    else if((b < a) && (a < c))
      {
        fabcit = fabcMap.find(std::make_tuple(b,a,c));
        return  (fabcit != fabcMap.end()) ? -fabcit->second : 0; 
      }
    else if((b < c) && (c < a))
      {
        fabcit = fabcMap.find(std::make_tuple(b,c,a));
        return  (fabcit != fabcMap.end()) ? fabcit->second : 0; 
      }
    else if((c < a) && (a < b))
      {
        fabcit = fabcMap.find(std::make_tuple(c,a,b));
        return  (fabcit != fabcMap.end()) ? fabcit->second : 0; 
      }
    else if((c < b) && (b < a))
      {
        fabcit = fabcMap.find(std::make_tuple(c,b,a));
        return  (fabcit != fabcMap.end()) ? - fabcit->second : 0; 
      }
    else
      throw std::runtime_error("Not found f^abc value");
  }

  double fCalc(size_t a, size_t b, size_t c)
  {
    return std::real(((t(a)*t(b) - t(b)*t(a))*t(c)).trace()*std::complex<double>(0,-2));
  }

  std::complex<double> tr(AdjVec av)
  {
    MatrixXcd m = MatrixXcd::Identity(nc, nc);

    for (AdjVec::iterator it = av.begin(); it != av.end(); ++it)
      m *= t(*it);
    return m.trace();
  }

  std::complex<double> trMap(const AdjVec& av)
  {
    MatrixXcd m = MatrixXcd::Identity(nc, nc);
    for (AdjVec::const_iterator it = av.begin(); it != av.end(); ++it)
      m *= t(amap[*it]);
    return m.trace();
  }
};



int main(int argc, char* argv[])
{

  try
    {
      if (argc != 2)
        {
          std::cout << "Number of colors as argument required!" << std::endl;
          std::cout << std::endl;
          std::cout << "Usage ./StupidPainter <Nc>" << std::endl;
          std::cout << std::endl;
          std::cout << "Syntax:" << std::endl; 
          std::cout << std::endl;
          std::cout << "  [T(1)T(1)]=tr[T^a1_ij*T^a1_ij]  - fundamental rep trace" << std::endl;
          std::cout << "  f(1,2,3)f(1,2,3)=f^{abc}f^{abc} - adjoint rep matrices" <<std::endl;
          std::cout << std::endl;

          return 1;
        }

      std::cout << "Number of colors : " << argv[1] << std::endl;

      int Nc(atoi(argv[1]));

      if(Nc < 2)
        throw std::logic_error("Allowed N>=2 only for SU(N)");

      Color col(Nc);

      char* input;
      std::string* lineInput;
      std::string shell_prompt("StupidPainter>");
      for(;;)
        {
          input = readline(KRED "StupidPainter>" KNRM);
          // Exit if received EOLN or readline return
          // zero pointer
          if (input == 0)
            break;

          lineInput = new std::string(input);


          // Add non empty input to the history
          if(!lineInput->empty())
            add_history(input);
          try
            {
              if(col.parse(*lineInput))
                {
                  Timer t1;
                  std::complex<double> res = col.contract();
                  std::stringstream resStr;
                  if ( res.imag() == 0 )
                    resStr << res.real();
                  else if( res.real() == 0 )
                    resStr << "I*" << res.imag();
                  else
                    resStr << res.real() << " + I*" << res.imag();
                  
                  std::cout << std::setfill ('=') << std::setw (80) << "" << std::endl;
                  std::cout << std::setfill (' ');
                  std::cout << "result:" << std::setw(33) << resStr.str()
                            << "  Time:"  << std::setw(30) << t1.elapsed() << " ms"
                            << std::endl;

                  // Clear sstream
                  std::stringstream().swap(resStr); 
                  // Try to rationalize result
                  long int reNum,reDen,imNum,imDen;
                  if(float2rat(res.real(), reNum, reDen) && float2rat(res.imag(), imNum, imDen))
                    {
                      if ( res.imag() == 0 )
                        {
                          resStr << reNum;
                          if (reDen != 1) resStr << "/" << reDen;
                        }
                      else if( res.real() == 0 )
                        {
                          resStr << "I*" << imNum;
                          if (imDen != 1) resStr << "/" << imDen;
                        }
                      else
                        {
                          resStr << reNum;
                          if (reDen != 1) resStr << "/" << reDen;
                          resStr << "I*" << imNum;
                          if (imDen != 1) resStr << "/" << imDen;
                        }
                      std::cout << "rationalized:" << std::setw(27) << resStr.str() << std::endl;
                    }
                  std::cout << std::endl;
                }
            }  catch (std::exception &p) 
            {
              std::cerr << p.what() << std::endl;
            }
          
          delete(lineInput);
          free(input);

        }

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }

  std::cout << "bye!" << std::endl;
  return 0;
}
