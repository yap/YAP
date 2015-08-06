/*  YAP - Yet another PWA toolkit
    Copyright 2015, Technische Universitaet Muenchen,
    Authors: Daniel Greenwald, Johannes Rauch

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef yap_SpinAmplitude_h
#define yap_SpinAmplitude_h

#include "DataAccessor.h"
#include "QuantumNumbers.h"
#include <array>

namespace yap {

class SpinAmplitude : public DataAccessor {
public:
  SpinAmplitude(const QuantumNumbers& initial, const QuantumNumbers& final1, const QuantumNumbers& final2)
  : InitialQuantumNumbers_(initial), FinalQuantumNumbers_({final1, final2}) {;}
  ~SpinAmplitude();

  virtual Amp amplitude(DataPoint& d);
  virtual bool checkConsistency() const;

private:
  QuantumNumbers InitialQuantumNumbers_;
  std::array<QuantumNumbers, 2> FinalQuantumNumbers_;

};

}

#endif
