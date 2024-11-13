#ifndef DUNE_COPASI_SOLVER_ISTL_UMFPACK_HH
#define DUNE_COPASI_SOLVER_ISTL_UMFPACK_HH

#if !HAVE_SUITESPARSE_UMFPACK
#error Trying to use umfpack wrapper without including umfpack
#endif

#include <dune/copasi/common/ostream_to_spdlog.hh>

#include <dune/istl/solver.hh>

#include <umfpack.h>

#include <algorithm>
#include <ranges>
#include <version>

// disabled by default due to memory leaks: https://github.com/oneapi-src/oneDPL/pull/1589
#if DUNE_COPASI_ENABLE_PARALLEL_SORT
#if __cpp_lib_execution >= 201603L
#include <execution>
#include <algorithm>
#elif __has_include(<oneapi/dpl/execution>)
#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#endif
#endif

namespace Dune::Copasi::ISTL {

template<Concept::AssembledLinearOperator O, class Alloc = std::allocator<void>>
class UMFPackWapper final : public InverseOperator<typename O::domain_type, typename O::range_type>
{
  using Matrix = typename O::matrix_type;
  using X = typename O::domain_type;
  using Y = typename O::range_type;
  using MI = PDELab::MultiIndex<std::size_t, blockLevel<Matrix>()>;

  static const inline std::map<std::string, double> umf_strategy_map = {
    { "AUTO", UMFPACK_STRATEGY_AUTO },
    { "UNSYMMETRIC", UMFPACK_STRATEGY_UNSYMMETRIC },
    { "UMFPACK_STRATEGY_SYMMETRIC", UMFPACK_STRATEGY_SYMMETRIC },
  };

  static const inline std::map<std::string, double> umf_ordering_map = {
    { "CHOLMOD", UMFPACK_ORDERING_CHOLMOD },
    { "AMD", UMFPACK_ORDERING_AMD },
    { "METIS", UMFPACK_ORDERING_METIS },
    { "BEST", UMFPACK_ORDERING_BEST },
    { "NONE", UMFPACK_ORDERING_NONE },
#if defined(UMFPACK_VER) && (UMFPACK >= UMFPACK_VER_CODE(6, 0))
    { "METIS_GUARD", UMFPACK_ORDERING_METIS_GUARD },
#endif
  };

  // rebind allocators to correct types
  using Int64Alloc = typename std::allocator_traits<Alloc>::template rebind_alloc<int64_t>;
  using DoubleAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<double>;

  struct SymbolicDeleter
  {
    void operator()(void* ptr) const { umfpack_dl_free_symbolic(&ptr); }
  };
  struct NumericDeleter
  {
    void operator()(void* ptr) const { umfpack_dl_free_numeric(&ptr); }
  };

public:
  // note allocator only allocates data for what UMFPack allows to be externally owned. Otherwise,
  // umprfack still allocates some memory internally
  UMFPackWapper(const std::shared_ptr<const O>& op,
                const ParameterTree& config,
                const Alloc& alloc = Alloc())
    : _alloc(alloc)
    , _umf_offsets(Int64Alloc(_alloc))
    , _umf_rows(Int64Alloc(_alloc))
    , _umf_iworkspace(Int64Alloc(_alloc))
    , _umf_dworkspace(DoubleAlloc(_alloc))
    , _umf_values(DoubleAlloc(_alloc))
    , _umf_x(DoubleAlloc(_alloc))
    , _umf_b(DoubleAlloc(_alloc))
    , _report_symbolic{config.get("report_control", false)}
    , _report_numeric{config.get("report_numeric", false)}
    , _report_info{config.get("report_info", true)}
  {
    // set up umfpack configuration options
    umfpack_di_defaults(_umf_control);
    auto set_config = [&](auto id, auto key) {
      _umf_control[id] = config.get(key, _umf_control[id]);
    };
    set_config(UMFPACK_PRL, "verbosity");
    set_config(UMFPACK_DENSE_ROW, "dense_row");
    set_config(UMFPACK_DENSE_COL, "dense_col");
    set_config(UMFPACK_PIVOT_TOLERANCE, "pivot_tolerance");
    set_config(UMFPACK_BLOCK_SIZE, "block_size");
    set_config(UMFPACK_DENSE_ROW, "dense_row");
    set_config(UMFPACK_ALLOC_INIT, "alloc_init");
    set_config(UMFPACK_IRSTEP, "irstep");
    set_config(UMFPACK_FIXQ, "fixq");
    set_config(UMFPACK_AMD_DENSE, "amd_dense");
    set_config(UMFPACK_SYM_PIVOT_TOLERANCE, "sym_pivot_tolerance");
    set_config(UMFPACK_SCALE, "scale");
    set_config(UMFPACK_FRONT_ALLOC_INIT, "front_alloc_init");
    set_config(UMFPACK_DROPTOL, "droptol");
    set_config(UMFPACK_AGGRESSIVE, "aggressive");
    set_config(UMFPACK_SINGLETONS, "singletons");
    _umf_control[UMFPACK_STRATEGY] = umf_strategy_map.at(config.get("strategy", "AUTO"));
    _umf_control[UMFPACK_ORDERING] = umf_ordering_map.at(config.get("ordering", "CHOLMOD"));
    // report about umfpack settings
    if (config.get("report_control", false))
      umfpack_dl_report_control(_umf_control);
    // store matrix into a format that unfpack understand
    storeFlatMatrix(op->getmat());
  }

  ~UMFPackWapper() override {
    if(_umf_numeric_ptr)
      umfpack_dl_free_numeric(&_umf_numeric_ptr);
  }

private:
  template<class M>
  void forEachMatrixEntry(M&& mat, const auto& call_back, MI row = MI(), MI col = MI())
  {
    using MD = std::decay_t<M>;
    if constexpr (IsNumber<MD>::value) {
      std::invoke(call_back, mat, row, col);
    } else { // assume row mayor access!
      for (auto row_it = mat.begin(); row_it != mat.end(); ++row_it) {
        auto this_row = push_back(row, row_it.index());
        for (auto col_it = row_it->begin(); col_it != row_it->end(); ++col_it) {
          forEachMatrixEntry(*col_it, call_back, this_row, push_back(col, col_it.index()));
        }
      }
    }
  }

  void storeFlatMatrix(const Matrix& mat)
  {
    // matrix entries as a (val, row, col) triplet
    struct Entry
    {
      double val;
      MI row;
      MI col;
    };
    using EntryAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<Entry>;
    EntryAlloc entryalloc(_alloc);

    // step 1: store sparse matrix in coordinate format (value, row, col)
    std::size_t nnz = 0;
    if (blockLevel<Matrix>() > 1) {
      forEachMatrixEntry(mat, [&](auto&&, auto&&, auto&&) { ++nnz; });
    } else {
      nnz = mat.nonzeroes();
    }

    auto entries = std::span<Entry>(std::allocator_traits<EntryAlloc>::allocate(entryalloc, nnz), nnz);
    nnz = 0;
    forEachMatrixEntry(mat, [&](double val, MI row, MI col) {
      std::allocator_traits<EntryAlloc>::construct(entryalloc, &entries[nnz++], val, row, col);
    });

    // step 2: flattern row multi-index into a single index
    if (blockLevel<Matrix>() > 1) {

      auto row_compare = [](auto lhs, auto rhs) {
        // note that sort does not need to order the col entry
        return std::ranges::lexicographical_compare(lhs.row, rhs.row);
      };
      // sort row multi-index lexicographically: block ordering of rows -> natural ordering of rows
      sort(
#if DUNE_COPASI_ENABLE_PARALLEL_SORT and (__cpp_lib_execution >= 201603L)
        std::execution::par,
#elif DUNE_COPASI_ENABLE_PARALLEL_SORT and __has_include(<oneapi/dpl/execution>)
        oneapi::dpl::execution::par,
#endif
        entries.begin(),
        entries.end(),
        row_compare);

      // convert multi-index into a single index
      MI last_row;
      std::size_t flat_row = -1;
      for (auto& [val, row, col] : entries) {
        // row is new [[pre: rows are sorted]]: increase flat row count
        if (row != std::exchange(last_row, row))
          ++flat_row;
        // store flat row in first index of row multi-index
        row.resize(1);
        row[0] = flat_row;
      }
      // number of flat rows
      _umf_n = flat_row + 1;
    } else {
      _umf_n = mat.N();
    }

    auto col_row_compare = [](auto lhs, auto rhs) {
      return std::ranges::lexicographical_compare(std::array{ lhs.col, lhs.row },
                                                  std::array{ rhs.col, rhs.row },
                                                  std::ranges::lexicographical_compare);
    };

    // order indices by column (note that sort does need to be order on row entry too)
    sort(
#if DUNE_COPASI_ENABLE_PARALLEL_SORT and (__cpp_lib_execution >= 201603L)
        std::execution::par,
#elif DUNE_COPASI_ENABLE_PARALLEL_SORT and __has_include(<oneapi/dpl/execution>)
        oneapi::dpl::execution::par,
#endif
      entries.begin(),
      entries.end(),
      col_row_compare);

    // initialize umfpack data
    _umf_values.clear();
    _umf_offsets.clear();
    _umf_rows.clear();
    _umf_offsets.reserve(entries.back().col[0]);
    _umf_rows.reserve(entries.size());
    _umf_values.reserve(entries.size());

    // transform col multi-index into a flat offset on the row flat indices
    MI last_col;
    for (std::size_t i = 0; i != entries.size(); ++i) {
      auto [val, row, col] = entries[i];
      // col is new [[pre: cols are sorted]]: store new column offset
      if (col != std::exchange(last_col, col))
        _umf_offsets.push_back(i);
      // store value and flat row index
      _umf_values.push_back(val);
      // [[pre: rows are sorted within column]]
      _umf_rows.push_back(static_cast<int64_t>(row[0]));
      std::allocator_traits<EntryAlloc>::destroy(entryalloc, &entries[i]);
    }
    _umf_m = _umf_offsets.size();
    _umf_offsets.push_back(entries.size());

    std::allocator_traits<EntryAlloc>::deallocate(entryalloc, entries.data(), entries.size());
  }

  void factorize()
  {
    assert(not _umf_offsets.empty());
    assert(not _umf_rows.empty());
    assert(not _umf_values.empty());

    // pre-order matrix columns
    void* symbolic_ptr = nullptr;
    umfpack_dl_symbolic(_umf_n,
                        _umf_m,
                        _umf_offsets.data(),
                        _umf_rows.data(),
                        _umf_values.data(),
                        &symbolic_ptr,
                        _umf_control,
                        _umf_info);
    if (_report_symbolic)
      umfpack_dl_report_symbolic(symbolic_ptr, _umf_control);

    // factorize matrix
    if(_umf_numeric_ptr)
      umfpack_dl_free_symbolic(&_umf_numeric_ptr);
    _umf_numeric_ptr = nullptr;
    umfpack_dl_numeric(_umf_offsets.data(),
                       _umf_rows.data(),
                       _umf_values.data(),
                       symbolic_ptr,
                       &_umf_numeric_ptr,
                       _umf_control,
                       _umf_info);
    // report results
    if (_report_numeric)
      umfpack_dl_report_numeric(_umf_numeric_ptr, _umf_control);
    umfpack_dl_free_symbolic(&symbolic_ptr);
  }

public:
  void apply(X& x, Y& b, InverseOperatorResult& res) override
  {
    // make sure matrix and buffers are setup
    if (not _umf_numeric_ptr)
      factorize();
    _umf_iworkspace.resize(_umf_n);
    _umf_dworkspace.resize(_umf_n * (_umf_control[UMFPACK_IRSTEP] > 0. ? 5 : 1));

    // assign b to working array
    _umf_b.resize(_umf_n);
    std::size_t flat_i = 0;
    PDELab::forEachContainerEntry(b, [&](auto val) { _umf_b[flat_i++] = val; });

    // reserve data for result
    _umf_x.resize(_umf_m);

    // actual solve
    int status = umfpack_dl_wsolve(UMFPACK_A,
                                   _umf_offsets.data(),
                                   _umf_rows.data(),
                                   _umf_values.data(),
                                   _umf_x.data(),
                                   _umf_b.data(),
                                   _umf_numeric_ptr,
                                   _umf_control,
                                   _umf_info,
                                   _umf_iworkspace.data(),
                                   _umf_dworkspace.data());
    // report results
    if (_report_info)
      umfpack_dl_report_info(_umf_control, _umf_info);

    // assign x from working array
    flat_i = 0;
    PDELab::forEachContainerEntry(x, [&](auto& val) { val = _umf_x[flat_i++]; });

    // store result info
    res.iterations = 1;
    res.converged = (status == 0);
  }

  void apply(X& x, Y& b, [[maybe_unused]] double reduction, InverseOperatorResult& res) override
  {
    this->apply(x, b, res);
  }

  SolverCategory::Category category() const override
  {
    return SolverCategory::Category::sequential;
  }

  static auto factory()
  {
    return [](const std::shared_ptr<LinearOperator<X, Y>>& op,
              const ParameterTree& config,
              const Alloc& alloc) -> std::shared_ptr<InverseOperator<X, Y>> {
      // make an allocator for this preconditioner
      using OpAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<UMFPackWapper>;
      OpAlloc opalloc(alloc);

      // cast operator to something the preconditioner expects
      auto aop = std::dynamic_pointer_cast<O>(op);
      if (not aop)
        throw format_exception(InvalidStateException{}, "Linear operator does not hold a matrix!");
      // construct operator instance
      return std::allocate_shared<UMFPackWapper>(opalloc, aop, config, alloc);
    };
  }

  int _umf_n, _umf_m;
  Alloc _alloc;
  std::vector<int64_t, Int64Alloc> _umf_offsets;
  std::vector<int64_t, Int64Alloc> _umf_rows;
  std::vector<int64_t, Int64Alloc> _umf_iworkspace;
  std::vector<double, DoubleAlloc> _umf_dworkspace;
  std::vector<double, DoubleAlloc> _umf_values;
  std::vector<double, DoubleAlloc> _umf_x;
  std::vector<double, DoubleAlloc> _umf_b;
  void* _umf_numeric_ptr = nullptr;
  double _umf_control[UMFPACK_CONTROL];
  double _umf_info[UMFPACK_INFO];
  bool _report_symbolic;
  bool _report_numeric;
  bool _report_info;
};

} // namespace Dune::Copasi::ISTL

#endif // DUNE_COPASI_SOLVER_ISTL_UMFPACK_HH
