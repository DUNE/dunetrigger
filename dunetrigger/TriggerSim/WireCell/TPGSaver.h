/** wclsTPGSaver

    Collects the trigger primitives produced by all wclsTPGFilter nodes of a
    Wire-Cell graph and puts them into the art::Event as a single
    std::vector<dunedaq::trgdataformats::TriggerPrimitive>.

    This component is not a graph node. Each wclsTPGFilter refers to it by
    "type:name" and calls add() (possibly concurrently, from TBB threads)
    while the graph runs. List it in the WCLS "outputers" so that its visit()
    is called after the graph has finished for the event.

    Configuration:
      instance  art product instance name (default "")
      sort      sort TPs by (time_start, channel, ...) before saving, so the
                output does not depend on thread scheduling (default true)
*/

#ifndef DUNETRIGGER_TRIGGERSIM_WIRECELL_TPGSAVER_H
#define DUNETRIGGER_TRIGGERSIM_WIRECELL_TPGSAVER_H

#include "larwirecell/Interfaces/IArtEventVisitor.h"

#include "WireCellAux/Logger.h"
#include "WireCellIface/IConfigurable.h"

#include "detdataformats/trigger/TriggerPrimitive.hpp"

#include <mutex>
#include <string>
#include <vector>

namespace dunetrigger {

  class TPGSaver : public WireCell::Aux::Logger,
                   public wcls::IArtEventVisitor,
                   public WireCell::IConfigurable {
  public:
    using tp_t = dunedaq::trgdataformats::TriggerPrimitive;

    TPGSaver();
    virtual ~TPGSaver();

    // IArtEventVisitor
    void produces(art::ProducesCollector& collector) override;
    void visit(art::Event& event) override;

    // IConfigurable
    void configure(const WireCell::Configuration& config) override;
    WireCell::Configuration default_configuration() const override;

    // Thread safe. Called by wclsTPGFilter nodes while the graph runs.
    void add(std::vector<tp_t>&& tps);

  private:
    std::string m_instance{""};
    bool m_sort{true};

    std::mutex m_mutex;
    std::vector<tp_t> m_tps;
  };

}

#endif
