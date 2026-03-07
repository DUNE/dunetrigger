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

/// Simulate one event: fill track and cluster writers
void simulate_event(SoAWriter<Track>&   tracks,
                    SoAWriter<Cluster>& clusters,
                    TRandom3&           rng)
{
    const int n_tracks   = rng.Integer(10) + 2;
    const int n_clusters = rng.Integer(6)  + 1;

    for (int t = 0; t < n_tracks; ++t) {
        // Set fields on the public staging row, then commit
        tracks.row.x          = (float)rng.Gaus(0, 0.05f);
        tracks.row.y          = (float)rng.Gaus(0, 0.05f);
        tracks.row.z          = (float)rng.Gaus(0, 5.0f);
        tracks.row.px         = (float)rng.Gaus(0, 1.0f);
        tracks.row.py         = (float)rng.Gaus(0, 1.0f);
        tracks.row.pz         = (float)rng.Gaus(0, 10.f);
        tracks.row.chi2       = (float)rng.Exp(1.0);
        tracks.row.n_hits     = (int)(rng.Integer(12) + 3);
        tracks.row.pdg_id     = (rng.Rndm() > 0.5) ? 211 : -211;
        tracks.row.is_primary = (rng.Rndm() > 0.3);
        tracks.push_back();   // commit row → SoA buffer, row is NOT reset automatically
    }

    for (int c = 0; c < n_clusters; ++c) {
        clusters.reset_row();  // explicit reset to avoid stale field values
        clusters.row.energy  = (float)rng.Exp(5.0);
        clusters.row.eta     = (float)rng.Uniform(-2.5, 2.5);
        clusters.row.phi     = (float)rng.Uniform(-M_PI, M_PI);
        clusters.row.n_cells = (int)(rng.Integer(20) + 1);
        clusters.push_back();
    }
}

// =============================================================================
//  Write demo
// =============================================================================
void write_demo(const char* filename) {
    std::cout << "\n=== WRITE ===\n";

    // -------------------------------------------------------------------------
    // 1. Create writers (once, outside the event loop)
    // -------------------------------------------------------------------------
    SoAWriter<Track>   track_writer(64);  // reserves 64 slots in the buffer
    SoAWriter<Cluster> cluster_writer(16);

    // -------------------------------------------------------------------------
    // 2. Create ROOT file + tree, register all branches automatically
    // -------------------------------------------------------------------------
    TFile file(filename, "RECREATE");
    TTree tree("events", "Simulated events");

    track_writer  .make_branches(tree, "trk_");
    cluster_writer.make_branches(tree, "cls_");

    track_writer  .print_summary();
    cluster_writer.print_summary();

    // -------------------------------------------------------------------------
    // 3. Event loop
    // -------------------------------------------------------------------------
    TRandom3 rng(42);
    constexpr int N_EVENTS = 100;

    for (int ev = 0; ev < N_EVENTS; ++ev) {
        // clear() empties the SoA buffer AND zero-initialises the staging row
        track_writer  .clear();
        cluster_writer.clear();

        simulate_event(track_writer, cluster_writer, rng);

        tree.Fill();

        if (ev < 3) {
            std::cout << "  Event " << ev
                      << "  tracks="   << track_writer.size()
                      << "  clusters=" << cluster_writer.size() << '\n';

            if (track_writer.size() > 0) {
                Track t0 = track_writer.buffer().get(0);
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
    track_buf  .set_branch_addresses(*tree, "trk_");
    cluster_buf.set_branch_addresses(*tree, "cls_");

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