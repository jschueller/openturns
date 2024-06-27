// SWIG file OTexceptions.i

%{
#include <csignal>
#include "openturns/Exception.hxx"
%}

namespace std {
  class exception
  {
  public:
    exception() throw() { }
    virtual ~exception() throw();
    /** Returns a C-style character string describing the general cause
     *  of the current error.  */
    virtual const char* what() const throw();
  };
}

%include exception.i

%exception {
  std::cout <<"OT%exception begin" << std::endl;
  const OT::Bool ok = (std::signal(SIGINT, OT::SignalHandler) != SIG_ERR);
  std::cout <<"OT%exception signal ok="<<ok << std::endl;
  try
  {
    std::cout <<"OT%exception try..."<< std::endl;
    $action
    std::cout <<"OT%exception try ok"<< std::endl;
  }
  catch (const OT::InvalidArgumentException & ex) {
    std::cout <<"OT%exception catch InvalidArgumentException" << std::endl;
    SWIG_exception(SWIG_TypeError, ex.__repr__().c_str());
  }
  catch (const OT::OutOfBoundException & ex) {
    std::cout <<"OT%exception catch OutOfBoundException" << std::endl;
    SWIG_exception(SWIG_IndexError, ex.__repr__().c_str());
  }
  catch (const OT::InterruptionException &) {
    std::cout <<"OT%exception catch InterruptionException" << std::endl;
    SWIG_exception(SWIG_RuntimeError, "Interruption in $decl");
  }
  catch (const OT::Exception & ex) {
    std::cout <<"OT%exception catch Exception" << std::endl;
    SWIG_exception(SWIG_RuntimeError, ex.__repr__().c_str());
  }
  catch (const std::range_error & ex) {
    std::cout <<"OT%exception catch range_error" << std::endl;
    SWIG_exception(SWIG_IndexError, ex.what());
  }
  catch (const std::out_of_range & ex) {
    std::cout <<"OT%exception catch out_of_range" << std::endl;
    SWIG_exception(SWIG_IndexError, ex.what());
  }
  catch (const std::logic_error & ex) {
    std::cout <<"OT%exception catch logic_error" << std::endl;
    SWIG_exception(SWIG_IndexError, ex.what());
  }
  catch (const std::exception & ex) {
    std::cout <<"OT%exception catch logic_error" << std::endl;
    SWIG_exception(SWIG_RuntimeError, ex.what());
  }
  catch(...)
  {
    std::cout <<"OT%exception catch ..." << std::endl;
    SWIG_exception(SWIG_RuntimeError, "Unknown exception");
  }
  std::cout <<"OT%exception end" << std::endl;
}
