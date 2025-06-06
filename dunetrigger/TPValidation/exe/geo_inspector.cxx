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

#include "TTree.h"
#include "TFile.h"

int main() {
    std::cout << "Geometry dumper" << std::endl;

    std::stringstream config;

    config << "#include \"geometry_dune.fcl\"" << std::endl;
    config << "services.Geometry:     @local::dunevd10kt_1x8x6_3view_30deg_geo" << std::endl;
    config << "services.WireReadout:      @local::dune_wire_readout" << std::endl;

    ArtServiceHelper::load_services(config);

    art::ServiceHandle<geo::Geometry> pgeo;
    auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();


    TFile f("tree.root", "recreate");

    TTree* tree = new TTree("a_ttree", "test_ttree");
    std::map<raw::ChannelID_t, std::vector<short>>  channel_buffer;


    std::cout<<"Total number of TPCs "<<pgeo->NTPC()<<std::endl;
    for (geo::TPCGeo const& tpc: pgeo->Iterate<geo::TPCGeo>(geo::CryostatID{0})) {
        size_t const t = tpc.ID().TPC;
        std::cout<<"TPC ID:" << t << " origin: " << tpc.GetCenter() << std::endl;

        std::vector<raw::ChannelID_t> channels;
        for (auto const& wire : wireReadout.Iterate<geo::WireID>(tpc.ID())) {
            raw::ChannelID_t ch = wireReadout.PlaneWireToChannel(wire);
            channels.push_back(ch);
           
            auto& buffer = channel_buffer[ch];
            tree->Branch(std::to_string(ch).c_str(), &buffer);
        }

        std::cout << "N Channels: " << channels.size() << "(First channel " << channels[0] << ")" << std::endl;
        
    }

    std::vector<short> my_adcs = {1,2,3};
    for( auto& [ch, adcs] : channel_buffer ) {
        adcs = my_adcs;
    }

    std::cout << "MaxTPCSet " << wireReadout.MaxTPCsets() << std::endl;

    tree->Fill();
    f.Write();
}