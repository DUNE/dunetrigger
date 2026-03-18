// dunetrigger/DataProducts/MyProduct.h
#ifndef DUNE_DATAPRODUCTS_MYPRODUCT_H
#define DUNE_DATAPRODUCTS_MYPRODUCT_H

#include <vector>

namespace dunetrigger {

  /**
   * @brief Stores a custom reconstructed quantity.
   *
   * Always document clearly what physics concept this class represents.
   * Immutable after construction: all data must be passed to the constructor.
   */
  class MyProduct {
  public:

    // Default constructor required by ROOT
    MyProduct() = default;

    // Main constructor: fully defines the object
    MyProduct(double energy, int nhits)
      : fEnergy(energy), fNHits(nhits) {}

    // Accessors (no setters — immutable design)
    double Energy() const { return fEnergy; }
    int    NHits()  const { return fNHits;  }

  private:
    double fEnergy = 0.0;
    int    fNHits  = 0;
  };

} // namespace dunetrigger

#endif // DUNE_DATAPRODUCTS_MYPRODUCT_H