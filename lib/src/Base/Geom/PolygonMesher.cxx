//                                               -*- C++ -*-
/**
 *  @brief Meshing algorithm for intervals
 *
 *  Copyright 2005-2024 Airbus-EDF-IMACS-ONERA-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/PolygonMesher.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/Exception.hxx"
#include "openturns/Tuples.hxx"
#include "openturns/KPermutations.hxx"
#include "openturns/SpecFunc.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(PolygonMesher)
static const Factory<PolygonMesher> Factory_PolygonMesher;


/* Default constructor */
PolygonMesher::PolygonMesher()
  : PersistentObject()
{
  // Nothing to do
}

/* Virtual constructor */
PolygonMesher * PolygonMesher::clone() const
{
  return new PolygonMesher(*this);
}

/* String converter */
String PolygonMesher::__repr__() const
{
  OSS oss(true);
  oss << "class=" << PolygonMesher::GetClassName();
  return oss;
}

/* String converter */
String PolygonMesher::__str__(const String & ) const
{
  return __repr__();
}

static Scalar calculate_normal(const Sample & polygon)
{
  Scalar sum = 0.0;
  const UnsignedInteger size = polygon.getSize();
  for (UnsignedInteger i = 0; i < size; ++ i)
  {
    const Scalar x0 = polygon(i, 0);
    const Scalar x1 = polygon((i + 1) % size, 0);
    const Scalar y0 = polygon(i, 1);
    const Scalar y1 = polygon((i + 1) % size, 1);
    sum += (x1 - x0) * (y1 + y0);
  }
  if (sum > SpecFunc::Precision)
    return 1.0;
  else if (sum < -SpecFunc::Precision)
    return -1.0;
  else
    throw InvalidArgumentException(HERE) << "PolygonMesher: No normal found";
}


Mesh PolygonMesher::build(const Sample & polygonVerticesC) const
{
    Sample polygonVertices(polygonVerticesC);

  const UnsignedInteger dimension = polygonVertices.getDimension();
  if (dimension != 2) throw InvalidArgumentException(HERE) << "Only 2-d";
  Sample vertices;
  IndicesCollection simplices;
  
  
//   const Scalar normal = calculate_normal(polygonVertices);
  UnsignedInteger i = 0;
  while (polygonVertices.getSize() > 2)
  {
    if (i >= polygonVertices.getSize())
      throw InvalidArgumentException(HERE) << "PolygonMesher: Triangulation failed";
    
    
  }
  
  
  
  
  
  
  return Mesh(vertices, simplices);
}

END_NAMESPACE_OPENTURNS

