#include "dunetrigger/triggeralgs/include/triggeralgs/dbscan/Hit.hpp"

#include <algorithm>

namespace triggeralgs {
namespace dbscan {

//======================================================================
HitSet::HitSet()
{
    hits.reserve(10);
}

//======================================================================
void
HitSet::insert(Hit* h)
{
    // We're typically inserting hits at or near the end, so do a
    // linear scan instead of full binary search. This turns out to be much
    // faster in our case
    auto it = hits.rbegin();
    while (it != hits.rend() && (*it)->time >= h->time) {
        // Don't insert the hit if we already have it
        if (*it == h) {
            return;
        }
        ++it;
    }
    
    if (it == hits.rend() || *it != h) {
        hits.insert(it.base(), h);
    }
}

//======================================================================
Hit::Hit(float _time, int _chan, const triggeralgs::TriggerPrimitive* _prim)
{
    reset(_time, _chan, _prim);
}

//======================================================================

void
Hit::reset(float _time, int _chan, const triggeralgs::TriggerPrimitive* _prim)
{
    time=_time;
    chan=_chan;
    cluster=kUndefined;
    connectedness=Connectedness::kUndefined;
    neighbours.clear();
    if(_prim){
        primitive=*_prim;
    }
}

//======================================================================

// Return true if hit was indeed a neighbour
bool
Hit::add_potential_neighbour(Hit* other, float eps, int minPts)
{
    if (other != this && euclidean_distance_sqr(*this, *other) < eps*eps) {
        neighbours.insert(other);
        if (neighbours.size() + 1 >= minPts) {
            connectedness = Connectedness::kCore;
        }
        // Neighbourliness is symmetric
        other->neighbours.insert(this);
        if (other->neighbours.size() + 1 >= minPts) {
            other->connectedness = Connectedness::kCore;
        }
        return true;
    }
    return false;
}

}
}
// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
