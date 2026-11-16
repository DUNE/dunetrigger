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
#include <nlohmann/json.hpp>
#include "geo_config.h"

using json = nlohmann::json;


int main(int argc, char** argv) {

    if (argc < 2 || argc > 3) {
        std::cout << "Usage: geo_inspector <geo> [wire_readout]\n\n"
                  << trgtools::shortcut_help() << std::endl;
        return -1;
    }

    std::cout << "Geometry Inspector" << std::endl;

    std::string geo_id = trgtools::resolve_geo(argv[1]);
    if (geo_id != argv[1]) {
        std::cout << "Geo ID: " << geo_id << std::endl;
    }

    std::string const readout_id = argc == 3 ? argv[2] : trgtools::default_wire_readout(geo_id);
    std::cout << "Wire readout: " << readout_id << std::endl;

    json j_geo;

    std::stringstream config(trgtools::service_config(geo_id, readout_id));

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

    j_geo["detector_name"] = geo->DetectorName();

    j_geo["cryostat"] = {
        {"origin", {
            {"x", cryogeo.GetCenter().x()},
            {"y", cryogeo.GetCenter().y()},
            {"z", cryogeo.GetCenter().z()}
        }},
        {"x_range", {
            {"min", cryogeo.MinX()},
            {"max", cryogeo.MaxX()}
        }},
        {"y_range", {
            {"min", cryogeo.MinY()},
            {"max", cryogeo.MaxY()}
        }},
        {"z_range", {
            {"min", cryogeo.MinZ()},
            {"max", cryogeo.MaxZ()}
        }}
    };

    int nTPC = pgeo->NTPC();
    std::cout<<"Total number of TPCs "<< nTPC <<std::endl;
    for (int i=0; i<nTPC; i++) {
        const auto& tpc = pgeo->Cryostat(geo::CryostatID{0}).TPC(i);

        size_t const id = tpc.ID().TPC;
        std::cout << "TPC ID : " << id << std::endl;
        std::cout << " - origin = " << tpc.GetCenter() << std::endl;
        std::cout << " - x-range = [" << tpc.MinX() << ", " << tpc.MaxX() << "]" << std::endl;
        std::cout << " - y-range = [" << tpc.MinY() << ", " << tpc.MaxY() << "]" << std::endl;
        std::cout << " - z-range = [" << tpc.MinZ() << ", " << tpc.MaxZ() << "]" << std::endl;


        j_geo["tpcs"][id] = {
            {"id", id},
            {"origin", {
                {"x", tpc.GetCenter().x()},
                {"y", tpc.GetCenter().y()},
                {"z", tpc.GetCenter().z()}
            }},
            {"x_range", {
                {"min", tpc.MinX()},
                {"max", tpc.MaxX()}
            }},
            {"y_range", {
                {"min", tpc.MinY()},
                {"max", tpc.MaxY()}
            }},
            {"z_range", {
                {"min", tpc.MinZ()},
                {"max", tpc.MaxZ()}
            }},
        };

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


    std::cout << "-----------" << std::endl;

    int nOpDets = pgeo->NOpDets();

    std::cout<<"Total number of OptDet "<<nOpDets<<std::endl;

    for (int i=0; i<nOpDets; i++) {
        const geo::OpDetGeo& od = pgeo->Cryostat(geo::CryostatID{0}).OpDet(i);
        size_t const id = od.ID().OpDet;
        std::cout << "OpDet ID : " << id << std::endl;
        std::cout << " - origin = " << od.GetCenter() << std::endl;
        std::cout << " - width =" << od.Width() << std::endl;
        std::cout << " - height ="  << od.Height() << std::endl;
        std::cout << " - length ="  << od.Length() << std::endl;


        j_geo["opdets"][i] = {
            {"id", id},
            {"origin", {
                {"x", od.GetCenter().x()},
                {"y", od.GetCenter().y()},
                {"z", od.GetCenter().z()}
            }},
            {"width", od.Width()},
            {"height", od.Height()},
            {"length", od.Length()}
        };
    }

    // for (geo::OpDetGeo const& od: pgeo->Iterate<geo::OpDetGeo>(geo::CryostatID{0})) {
    //         // size_t const t = od.ID().OpDet;

    // }




    std::cout << "-----------" << std::endl;
    std::cout << std::setw(4) << j_geo << std::endl;



    // write prettified JSON to another file
    std::string json_geo_file = geo_id+".json";
    std::cout << "Saving json geo to " << json_geo_file << std::endl;
    std::ofstream o(json_geo_file);
    o << std::setw(4) << j_geo << std::endl;

}