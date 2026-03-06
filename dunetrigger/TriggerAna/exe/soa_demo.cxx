// =============================================================================
//  main.cpp – usage demo for SoABuffer with ROOT TTree
//
//  Compile:
//    g++ -std=c++17 main.cpp $(root-config --cflags --libs) -I/path/to/boost -o soa_demo
// =============================================================================

#include "../SoABuffer.hpp"

#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>

#include <iostream>
#include <cassert>

// =============================================================================
//  Domain structs – plain POD, zero modification needed
// =============================================================================

struct Track {
    float  x;        // vertex x  [cm]
    float  y;        // vertex y  [cm]
    float  z;        // vertex z  [cm]
    float  px;       // momentum x [GeV/c]
    float  py;       // momentum y [GeV/c]
    float  pz;       // momentum z [GeV/c]
    float  chi2;     // track fit chi²
    int    n_hits;   // number of detector hits
    int    pdg_id;   // PDG particle code
    bool   is_primary;
};

struct Cluster {
    float  energy;   // deposited energy [GeV]
    float  eta;      // pseudorapidity
    float  phi;      // azimuth [rad]
    int    n_cells;  // number of calorimeter cells
};

// Register field names for C++17 mode (no-op in C++20, names come from PFR).
REGISTER_SOA_FIELD_NAMES(Track,   x, y, z, px, py, pz, chi2, n_hits, pdg_id, is_primary)
REGISTER_SOA_FIELD_NAMES(Cluster, energy, eta, phi, n_cells)

// =============================================================================
//  Helpers
// =============================================================================

/// Simulate one event: fill track and cluster buffers
void simulate_event(SoABuffer<Track>&   tracks,
                    SoABuffer<Cluster>& clusters,
                    TRandom3&           rng,
                    int                 event_id)
{
    const int n_tracks   = rng.Integer(10) + 2;
    const int n_clusters = rng.Integer(6)  + 1;

    for (int t = 0; t < n_tracks; ++t) {
        Track tr;
        tr.x          = (float)rng.Gaus(0, 0.05f);
        tr.y          = (float)rng.Gaus(0, 0.05f);
        tr.z          = (float)rng.Gaus(0, 5.0f);
        tr.px         = (float)rng.Gaus(0, 1.0f);
        tr.py         = (float)rng.Gaus(0, 1.0f);
        tr.pz         = (float)rng.Gaus(0, 10.f);
        tr.chi2       = (float)rng.Exp(1.0);
        tr.n_hits     = (int)(rng.Integer(12) + 3);
        tr.pdg_id     = (rng.Rndm() > 0.5) ? 211 : -211;
        tr.is_primary = (rng.Rndm() > 0.3);
        tracks.push_back(tr);
    }

    for (int c = 0; c < n_clusters; ++c) {
        Cluster cl;
        cl.energy  = (float)rng.Exp(5.0);
        cl.eta     = (float)rng.Uniform(-2.5, 2.5);
        cl.phi     = (float)rng.Uniform(-M_PI, M_PI);
        cl.n_cells = (int)(rng.Integer(20) + 1);
        clusters.push_back(cl);
    }
}

// =============================================================================
//  Write demo
// =============================================================================
void write_demo(const char* filename) {
    std::cout << "\n=== WRITE ===\n";

    // -------------------------------------------------------------------------
    // 1. Create buffers (once, outside the event loop)
    // -------------------------------------------------------------------------
    SoABuffer<Track>   track_buf;
    SoABuffer<Cluster> cluster_buf;

    track_buf  .reserve(64);   // optional performance hint
    cluster_buf.reserve(16);

    // -------------------------------------------------------------------------
    // 2. Create ROOT file + tree, register all branches automatically
    // -------------------------------------------------------------------------
    TFile file(filename, "RECREATE");
    TTree tree("events", "Simulated events");

    // Each field of Track  gets a branch named "trk_<fieldname>"
    // Each field of Cluster gets a branch named "cls_<fieldname>"
    track_buf  .make_branches(&tree, "trk_");
    cluster_buf.make_branches(&tree, "cls_");

    // Print what was registered
    track_buf  .print_summary();
    cluster_buf.print_summary();

    // -------------------------------------------------------------------------
    // 3. Event loop
    // -------------------------------------------------------------------------
    TRandom3 rng(42);
    constexpr int N_EVENTS = 100;

    for (int ev = 0; ev < N_EVENTS; ++ev) {
        // --- clear buffers at the start of each event ---
        track_buf  .clear();
        cluster_buf.clear();

        // --- fill with simulated data ---
        simulate_event(track_buf, cluster_buf, rng, ev);

        // --- fill the tree (buffers are still owned by SoABuffer) ---
        tree.Fill();

        if (ev < 3) {
            std::cout << "  Event " << ev
                      << "  tracks="   << track_buf.size()
                      << "  clusters=" << cluster_buf.size() << '\n';

            // Reconstruct first track as AoS struct for illustration
            if (track_buf.size() > 0) {
                Track t0 = track_buf.get(0);
                std::cout << "    track[0]: px=" << t0.px
                          << " py=" << t0.py
                          << " pdg=" << t0.pdg_id << '\n';
            }
        }
    }

    // -------------------------------------------------------------------------
    // 4. Write and close
    // -------------------------------------------------------------------------
    file.Write();
    file.Close();
    std::cout << "Wrote " << N_EVENTS << " events to " << filename << '\n';
}

// =============================================================================
//  Read-back demo
// =============================================================================
void read_demo(const char* filename) {
    std::cout << "\n=== READ BACK ===\n";

    SoABuffer<Track>   track_buf;
    SoABuffer<Cluster> cluster_buf;

    TFile file(filename, "READ");
    TTree* tree = nullptr;
    file.GetObject("events", tree);
    assert(tree && "Tree not found!");

    // Re-point branch addresses to our new buffer instances
    track_buf  .set_branch_addresses(tree, "trk_");
    cluster_buf.set_branch_addresses(tree, "cls_");

    const Long64_t n_entries = tree->GetEntries();
    std::cout << "Reading " << n_entries << " events\n";

    double total_tracks = 0;
    for (Long64_t ev = 0; ev < n_entries; ++ev) {
        track_buf  .clear();
        cluster_buf.clear();

        tree->GetEntry(ev);
        total_tracks += track_buf.size();
    }

    std::cout << "Average tracks/event: "
              << total_tracks / n_entries << '\n';
    file.Close();
}

// =============================================================================
//  Column-level access demo
// =============================================================================
void column_access_demo() {
    std::cout << "\n=== COLUMN ACCESS ===\n";

    SoABuffer<Track> buf;
    TRandom3 rng(7);

    for (int i = 0; i < 5; ++i) {
        Track tr;
        tr.x=0; tr.y=0; tr.z=0;
        tr.px=(float)rng.Gaus(0,1);
        tr.py=(float)rng.Gaus(0,1);
        tr.pz=(float)rng.Gaus(0,10);
        tr.chi2=1.f; tr.n_hits=8; tr.pdg_id=211; tr.is_primary=true;
        buf.push_back(tr);
    }

    // Direct SIMD-friendly access to the px column (index 3)
    auto& px_vec = buf.column<3>();   // std::vector<float>&
    std::cout << "px values: ";
    for (float v : px_vec) std::cout << v << ' ';
    std::cout << '\n';

    // Field names at compile time
    auto names = SoABuffer<Track>::field_names();
    std::cout << "Field names: ";
    for (auto n : names) std::cout << n << ' ';
    std::cout << '\n';
}

// =============================================================================
//  main
// =============================================================================
int main() {
    const char* fname = "soa_demo.root";

    write_demo(fname);
    read_demo(fname);
    column_access_demo();

    return 0;
}