// SWIG file RandomGenerator.i

%{
#include "openturns/RandomGenerator.hxx"

namespace OT {

template <>
struct traitsPythonType<OT::RandomGenerator>
{
  typedef _PyObject_ Type;
};

template <>
inline
OT::RandomGenerator
convert<_PyObject_, OT::RandomGenerator>(PyObject * pyObj)
{
  void * ptr = 0;
  if (SWIG_IsOK(SWIG_ConvertPtr(pyObj, &ptr, SWIG_TypeQuery("OT::RandomGenerator *"), SWIG_POINTER_NO_NULL))) {
    OT::RandomGenerator * p_it = reinterpret_cast< OT::RandomGenerator * >(ptr);
    return *p_it;
  }
  else {
    throw OT::InvalidArgumentException(HERE) << "Object passed as argument is not convertible to a RandomGenerator";
  }
}

} // namespace OT
%}

%include RandomGenerator_doc.i

%copyctor OT::RandomGenerator;

%include openturns/RandomGenerator.hxx

%template(RandomGeneratorCollection) OT::Collection<OT::RandomGenerator>;

namespace OT {}
