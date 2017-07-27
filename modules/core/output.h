#include <streambuf>
#include <ostream>

class mstream : public std::streambuf {
public:
protected:
  virtual std::streamsize xsputn(const char *s, std::streamsize n)
  {
      mexPrintf("%.*s",n,s);
      return n;        
  }
  
  virtual int overflow(int c = EOF)
  {
    if (c != EOF) {
      mexPrintf("%.1s",&c);
    }
    return 1;
  }
}; 
