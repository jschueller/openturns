// SWIG file DualBootstrapExperiment.i

// do not pass argument by reference, return it as tuple item
%typemap(in, numinputs=0) OT::Sample & sample2Out ($*ltype temp) %{ temp = OT::Sample(); $1 = &temp; %}
%typemap(argout) OT::Sample & sample2Out %{ $result = SWIG_Python_AppendOutput($result, SWIG_NewPointerObj(new OT::Sample(*$1), SWIGTYPE_p_OT__Sample, SWIG_POINTER_OWN |  0 )); %}

%{
#include "openturns/DualBootstrapExperiment.hxx"
%}

%include DualBootstrapExperiment_doc.i

%include openturns/DualBootstrapExperiment.hxx
namespace OT { %extend DualBootstrapExperiment { DualBootstrapExperiment(const DualBootstrapExperiment & other) { return new OT::DualBootstrapExperiment(other); } } }
