#ifndef SOA_BUFFER_HPP
#define SOA_BUFFER_HPP
// =============================================================================
//  SoABuffer.hpp
//  Structure-of-Arrays storage buffer and write-staging wrapper, auto-reflected
//  from a C++ struct via Boost.PFR.
//
//  Requirements: C++20, Boost >= 1.80 (Boost.PFR), ROOT >= 6.x
//
//  Classes:
//    SoABuffer<Struct>  — raw SoA storage + ROOT TTree interface
//    SoAWriter<Struct>  — SoABuffer + staging row (convenience for writing)
// =============================================================================

#include "FieldNames.hh"

#include <boost/pfr.hpp>
#include <TTree.h>

#include <array>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <vector>

// ---------------------------------------------------------------------------
//  SoABuffer<Struct>
//  ---
//  Wraps one std::vector<T> per field of Struct.  Provides:
//    push_back(const Struct&)              — append one element
//    get(size_t i)                         — reconstruct an AoS element
//    size() / clear() / reserve()
//    make_branches(tree, prefix)           — register vector branches (write)
//    set_branch_addresses(tree, prefix)    — attach to existing branches (read)
// ---------------------------------------------------------------------------
template<trg_concepts::PfrAggregate Struct>
    requires std::is_default_constructible_v<Struct>
          && std::is_copy_constructible_v<Struct>
          && (boost::pfr::tuple_size_v<Struct> > 0)
class SoABuffer {
public:
    static constexpr std::size_t kNFields = boost::pfr::tuple_size_v<Struct>;

    using FieldTuple  = decltype(boost::pfr::structure_to_tuple(std::declval<Struct>()));
    using ArraysTuple = typename trg_detail::Vectorize<FieldTuple>::arrays;
    using PtrsTuple   = typename trg_detail::Vectorize<FieldTuple>::ptrs;

    // ------------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------------

    SoABuffer() : ptrs_(rebind_ptrs()) {}

    SoABuffer(const SoABuffer&)            = delete;
    SoABuffer& operator=(const SoABuffer&) = delete;

    SoABuffer(SoABuffer&& o) noexcept
        : arrays_(std::move(o.arrays_)), ptrs_(rebind_ptrs()) {}

    SoABuffer& operator=(SoABuffer&& o) noexcept {
        if (this != &o) {
            arrays_ = std::move(o.arrays_);
            ptrs_   = rebind_ptrs();
        }
        return *this;
    }

    // ------------------------------------------------------------------
    // Element access
    // ------------------------------------------------------------------

    /// @brief Append a struct as a new row.
    void push_back(const Struct& s) {
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            std::get<I>(arrays_).push_back(boost::pfr::get<I>(s));
        });
    }

    /// @brief Reconstruct a struct from row @p i.
    /// @throws std::out_of_range if @p i >= size().
    Struct get(std::size_t i) const {
        if (i >= size())
            throw std::out_of_range(
                "SoABuffer::get: index " + std::to_string(i) +
                " out of range (size=" + std::to_string(size()) + ")");
        Struct s{};
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            boost::pfr::get<I>(s) = std::get<I>(arrays_)[i];
        });
        return s;
    }

    /// @return Number of rows stored.
    [[nodiscard]] std::size_t size() const { return std::get<0>(arrays_).size(); }

    /// @brief Direct access to the I-th column vector.
    template<std::size_t I>
    auto& column() { return std::get<I>(arrays_); }

    /// @copydoc column()
    template<std::size_t I>
    const auto& column() const { return std::get<I>(arrays_); }

    // ------------------------------------------------------------------
    // Bulk operations
    // ------------------------------------------------------------------

    /// @brief Pre-allocate all column vectors.
    void reserve(std::size_t n) {
        std::apply([n](auto&... vecs) { (vecs.reserve(n), ...); }, arrays_);
    }

    /// @brief Empty all column vectors (capacity is preserved).
    void clear() {
        std::apply([](auto&... vecs) { (vecs.clear(), ...); }, arrays_);
    }

    // ------------------------------------------------------------------
    // ROOT TTree interface
    // ------------------------------------------------------------------

    /// @brief Create one STL-vector branch per field on @p tree.
    /// @param prefix Optional prefix prepended to each branch name.
    void make_branches(TTree& tree, const std::string& prefix = "") {
        constexpr auto names = trg_detail::get_field_names_sv<Struct>();
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            tree.Branch((prefix + names[I]).c_str(), &std::get<I>(arrays_));
        });
    }

    /// @brief Re-point branch addresses to this buffer's vectors.
    /// @throws std::runtime_error if a branch cannot be bound.
    void set_branch_addresses(TTree& tree, const std::string& prefix = "") {
        constexpr auto names = trg_detail::get_field_names_sv<Struct>();
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            const auto name = prefix + names[I];
            check_set_address(tree.SetBranchAddress(name.c_str(), &std::get<I>(ptrs_)), name);
        });
    }

    // ------------------------------------------------------------------
    // Utilities
    // ------------------------------------------------------------------

    /// @brief Returns the compile-time array of field name string_views.
    [[nodiscard]] static constexpr std::array<std::string_view, kNFields>
    field_names() noexcept {
        return trg_detail::get_field_names_sv<Struct>();
    }

    /// @brief Print a short summary of column names and sizes to @p os.
    void print_summary(std::ostream& os = std::cout) const {
        constexpr auto names = trg_detail::get_field_names_sv<Struct>();
        os << "SoABuffer<" << typeid(Struct).name()
           << ">  rows=" << size() << "  fields=" << kNFields << '\n';
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            os << "  [" << I << "] " << names[I]
               << "  size=" << std::get<I>(arrays_).size() << '\n';
        });
    }

private:
    ArraysTuple arrays_;
    PtrsTuple   ptrs_;

    // Re-points each element of ptrs_ at the corresponding vector in arrays_.
    // Called after construction and after any move.
    PtrsTuple rebind_ptrs() noexcept {
        PtrsTuple p;
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            std::get<I>(p) = &std::get<I>(arrays_);
        });
        return p;
    }

    static void check_set_address(Int_t status, const std::string& branch_name) {
        if (status < 0)
            throw std::runtime_error(
                "SoABuffer::set_branch_addresses: failed to bind branch \"" +
                branch_name + "\" (ROOT error code " + std::to_string(status) + ")");
    }
};

// ---------------------------------------------------------------------------
//  SoAWriter<Struct>
//  ---
//  Convenience wrapper around SoABuffer<Struct> with an internal staging row.
//
//  @par Typical write loop
//  @code
//  SoAWriter<Track> writer;
//  writer.make_branches(tree, "trk_");
//
//  for (auto& raw : source) {
//      writer->x  = raw.x;      // fill staging row
//      writer->px = raw.px;
//      writer.push_back();      // commit row → SoA storage
//  }
//  tree.Fill();
//  writer.clear();              // reset storage and staging row
//  @endcode
// ---------------------------------------------------------------------------
template<trg_concepts::PfrAggregate Struct>
    requires std::is_default_constructible_v<Struct>
          && std::is_copy_constructible_v<Struct>
          && (boost::pfr::tuple_size_v<Struct> > 0)
class SoAWriter {
public:
    /// @brief Access the staging row fields via `writer->field`.
    Struct* operator->() noexcept       { return &row_; }
    const Struct* operator->() const noexcept { return &row_; }

    Struct& row() noexcept             { return row_; }
    const Struct& row() const noexcept { return row_; }

    SoAWriter() = default;

    explicit SoAWriter(std::size_t reserve_n) { buffer_.reserve(reserve_n); }

    // ------------------------------------------------------------------
    // Enable / disable
    // ------------------------------------------------------------------

    /// @brief Activate or deactivate this writer (default: enabled).
    /// @note When disabled, make_branches, push_back, commit_and_reset,
    ///       and clear are all no-ops.  Must be called before make_branches.
    void enable(bool e = true) noexcept { enabled_ = e; }

    /// @return @c true if write operations are active.
    [[nodiscard]] bool is_enabled() const noexcept { return enabled_; }

    /// @copydoc is_enabled()
    [[nodiscard]] explicit operator bool() const noexcept { return enabled_; }

    // ------------------------------------------------------------------
    // Core operations
    // ------------------------------------------------------------------

    /// @brief Commit the current staging row to storage.
    /// @note No-op when disabled.
    void push_back() {
        if (enabled_) buffer_.push_back(row_);
    }

    /// @brief Commit staging row and zero-initialise it.
    /// @note No-op when disabled.
    void commit_and_reset() {
        if (enabled_) { buffer_.push_back(row_); reset_row(); }
    }

    /// @brief Empty storage and zero-initialise the staging row.
    /// @note No-op when disabled.
    void clear() {
        if (enabled_) { buffer_.clear(); reset_row(); }
    }

    /// @brief Zero-initialise the staging row without touching the buffer.
    void reset_row() { row_ = Struct{}; }

    // ------------------------------------------------------------------
    // Buffer access
    // ------------------------------------------------------------------

    [[nodiscard]] SoABuffer<Struct>&       buffer() noexcept       { return buffer_; }
    [[nodiscard]] const SoABuffer<Struct>& buffer() const noexcept { return buffer_; }

    /// @return Number of committed rows.
    [[nodiscard]] std::size_t size() const noexcept { return buffer_.size(); }

    // ------------------------------------------------------------------
    // ROOT TTree forwarding
    // ------------------------------------------------------------------

    /// @brief Register branches on @p tree.  No-op when disabled.
    void make_branches(TTree& tree, const std::string& prefix = "") {
        if (enabled_) buffer_.make_branches(tree, prefix);
    }

    void print_summary(std::ostream& os = std::cout) const { buffer_.print_summary(os); }

private:
    bool             enabled_ = true;
    Struct           row_{};
    SoABuffer<Struct> buffer_;
};

#endif // SOA_BUFFER_HPP
