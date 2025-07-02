#include <iostream>


#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "dunecore/ArtSupport/ArtServiceHelper.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcorealg/Geometry/Exceptions.h"

// #include "TTree.h"
// #include "TFile.h"

int main() {
    std::cout << "Geometry Inspector" << std::endl;

    std::stringstream config;

    config << "#include \"geometry_dune.fcl\"" << std::endl;
    config << "services.Geometry:     @local::dunevd10kt_1x8x6_3view_30deg_geo" << std::endl;
    config << "services.WireReadout:      @local::dune_wire_readout" << std::endl;

    ArtServiceHelper::load_services(config);

    art::ServiceHandle<geo::Geometry> pgeo;
    auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    auto const& geo = art::ServiceHandle<geo::Geometry>();

    std::cout << geo->Info() << std::endl;
    std::cout << "\n\n---------------------------------------------------------------------\n\n" << std::endl;

    std::cout << "Detector Name : " << geo->DetectorName() << std::endl;

    auto const& cryogeo = geo->Cryostats().front();
    std::cout << "Cryostat ID : " << cryogeo.ID() << std::endl;
    std::cout << " - origin = " << cryogeo.GetCenter() << std::endl;
    std::cout << " - x-range = [" << cryogeo.MinX() << ", " << cryogeo.MaxX() << "]" << std::endl;
    std::cout << " - y-range = [" << cryogeo.MinY() << ", " << cryogeo.MaxY() << "]" << std::endl;
    std::cout << " - z-range = [" << cryogeo.MinZ() << ", " << cryogeo.MaxZ() << "]" << std::endl;

    std::cout<<"Total number of TPCs "<<pgeo->NTPC()<<std::endl;
    for (geo::TPCGeo const& tpc: pgeo->Iterate<geo::TPCGeo>(geo::CryostatID{0})) {
        size_t const t = tpc.ID().TPC;
        std::cout << "TPC ID : " << t << std::endl;
        std::cout << " - origin = " << tpc.GetCenter() << std::endl;
        std::cout << " - x-range = [" << tpc.MinX() << ", " << tpc.MaxX() << "]" << std::endl;
        std::cout << " - y-range = [" << tpc.MinY() << ", " << tpc.MaxY() << "]" << std::endl;
        std::cout << " - z-range = [" << tpc.MinZ() << ", " << tpc.MaxZ() << "]" << std::endl;


        std::vector<raw::ChannelID_t> channels;
        for (auto const& wire : wireReadout.Iterate<geo::WireID>(tpc.ID())) {
            raw::ChannelID_t ch = wireReadout.PlaneWireToChannel(wire);
            channels.push_back(ch);

            // auto& wiregeo = wireReadout.Wire(wire);
            // std::cout << " - " <<  ch << " " << wiregeo.Direction() << std::endl;
           
        }

        std::cout << " - N Channels = " << channels.size() << " (First channel =" << channels[0] << ")" << std::endl;
        
    }

    std::cout << "MaxTPCSet " << wireReadout.MaxTPCsets() << std::endl;

}