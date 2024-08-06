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
#ifndef OPENTURNS_POLYGONMESHER_HXX
#define OPENTURNS_POLYGONMESHER_HXX

#include "openturns/Interval.hxx"
#include "openturns/Mesh.hxx"
#include "openturns/ResourceMap.hxx"

BEGIN_NAMESPACE_OPENTURNS

/**
 * @class PolygonMesher
 */
class OT_API PolygonMesher
  : public PersistentObject
{
  CLASSNAME
public:

  /** Default constructor */
  PolygonMesher();

  /** Virtual constructor */
  PolygonMesher * clone() const override;

  /** String converter */
  String __repr__() const override;

  /** String converter */
  String __str__(const String & offset = "") const override;

  /* Here is the interface that all derived class must implement */
  virtual Mesh build(const Sample & vertices) const;

protected:

private:

}; /* class PolygonMesher */

END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_POLYGONMESHER_HXX */
