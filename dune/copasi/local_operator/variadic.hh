#ifndef DUNE_COPASI_LOCAL_OPERATOR_VARIADIC_HH
#define DUNE_COPASI_LOCAL_OPERATOR_VARIADIC_HH

#include <tuple>
#include <type_traits>

namespace Dune::Copasi {

//! Helper lambda functions which apply functions to local operators
namespace LocalOperatorApply {

// Pattern helepers

static auto patternVolume = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doPatternVolume)
    lop.pattern_volume(args...);
};

static auto patternVolumePostSkeleton = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doPatternVolumePostSkeleton)
    lop.pattern_volume_post_skeleton(args...);
};

static auto patternSkeleton = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doPatternSkeleton)
    lop.pattern_skeleton(args...);
};

static auto patternBoundary = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doPatternBoundary)
    lop.pattern_boundary(args...);
};

// Alpha helepers

static auto alphaVolume = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doAlphaVolume)
    lop.alpha_volume(args...);
};

static auto alphaVolumePostSkeleton = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doAlphaVolumePostSkeleton)
    lop.alpha_volume_post_skeleton(args...);
};

static auto alphaSkeleton = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doAlphaSkeleton)
    lop.alpha_skeleton(args...);
};

static auto alphaBoundary = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doAlphaBoundary)
    lop.alpha_boundary(args...);
};


// Lambda helepers

static auto lambdaVolume = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doLambdaVolume)
    lop.lambda_volume(args...);
};

static auto lambdaVolumePostSkeleton = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doLambdaVolumePostSkeleton)
    lop.lambda_volume_post_skeleton(args...);
};

static auto lambdaSkeleton = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doLambdaSkeleton)
    lop.lambda_skeleton(args...);
};

static auto lambdaBoundary = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doLambdaBoundary)
    lop.lambda_boundary(args...);
};


// Jacobian helepers

static auto jacobianVolume = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doAlphaVolume)
    lop.jacobian_volume(args...);
};

static auto jacobianVolumePostSkeleton = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doAlphaVolumePostSkeleton)
    lop.jacobian_volume_post_skeleton(args...);
};

static auto jacobianSkeleton = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doAlphaSkeleton)
    lop.jacobian_skeleton(args...);
};

static auto jacobianBoundary = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doAlphaBoundary)
    lop.jacobian_boundary(args...);
};


// Jacobian apply helepers

static auto jacobianApplyVolume = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doAlphaVolume)
    lop.jacobian_apply_volume(args...);
};

static auto jacobianApplyVolumePostSkeleton = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doAlphaVolumePostSkeleton)
    lop.jacobian_apply_volume_post_skeleton(args...);
};

static auto jacobianApplySkeleton = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doAlphaSkeleton)
    lop.jacobian_apply_skeleton(args...);
};

static auto jacobianApplyBoundary = [](const auto& lop, auto&... args)
{
  using LOP = std::decay_t<decltype(lop)>;
  if constexpr (LOP::doAlphaBoundary)
    lop.jacobian_apply_boundary(args...);
};

} // namespace LocalOperatorApply

// this class helps to forward different (dynamic) finite elements to different (static) local operators.
// In case that skeleton terms appear on different operators (e.g. inside -> DGP2, outside-> DGP1) such
// local operators should implement one sided skeletons
// TODO: fix this problem by zeoring outside part and call outside part.
// Local operators should only do work for its type of finite element!
template<class FiniteElementMapper, class... LocalOperators>
class VariadicLocalOperator
{
  template<bool... Args>
  bool constexpr static disjunction()
  {
    return std::disjunction_v<std::bool_constant<Args>...>;
  }

  template<bool... Args>
  bool constexpr static conjunction()
  {
    return std::conjunction_v<std::bool_constant<Args>...>;
  }

  template<class P, class F, class... Args>
  void applyLops(P && p, F&& f, Args&&... args) const
  {
    Hybrid::forEach(Dune::range(_integral_size{}),
      [&](auto i){
        if (p(i))
          f(*Hybrid::elementAt(_lops, i), args...);
      }
    );
  }

  struct BilinearMapperVisitor
    : public TypeTree::TreePairVisitor
    , public TypeTree::DynamicTraversal
  {
    BilinearMapperVisitor(const FiniteElementMapper& fe_mapper)
      : _fe_mapper(fe_mapper)
    {}

    template<typename LFSU, typename LFSV, typename TreePath>
    void leaf(LFSU&& lfsu, LFSV&& lfsv, TreePath treePath) const
    {
      _indices.insert(_fe_mapper(lfsu.finiteElement(),lfsv.finiteElement()));
    }

    mutable std::set<std::size_t> _indices;
    const FiniteElementMapper& _fe_mapper;
  };

  struct LinearMapperVisitor
    : public TypeTree::TreeVisitor
    , public TypeTree::DynamicTraversal
  {
    LinearMapperVisitor(const FiniteElementMapper& fe_mapper)
      : _fe_mapper(fe_mapper)
    {}

    template<typename LFSV, typename TreePath>
    void leaf(LFSV&& lfsv, TreePath treePath) const
    {
      _indices.insert(_fe_mapper(lfsv.finiteElement()));
    }

    mutable std::set<std::size_t> _indices;
    const FiniteElementMapper& _fe_mapper;
  };

  template<class LFSU, class LFSV>
  auto indices(const LFSU& lfsu, const LFSV& lfsv) const
  {
    BilinearMapperVisitor visitor(_fe_mapper);
    TypeTree::applyToTreePair(lfsu,lfsv,visitor);
    return visitor._indices;
  }

  // todo: cache indices in the pre volume and skeleton stage
  template<class LFSV>
  auto indices(const LFSV& lfsv) const
  {
    LinearMapperVisitor visitor(_fe_mapper);
    TypeTree::applyToTree(lfsv,visitor);
    return visitor._indices;
  }

public:
  enum {doPatternVolume             = disjunction<LocalOperators::doPatternVolume...>()};
  enum {doPatternVolumePostSkeleton = disjunction<LocalOperators::doPatternVolumePostSkeleton...>()};
  enum {doPatternSkeleton           = disjunction<LocalOperators::doPatternSkeleton...>()};
  enum {doPatternBoundary           = disjunction<LocalOperators::doPatternBoundary...>()};
  enum {doAlphaVolume               = disjunction<LocalOperators::doAlphaVolume...>()};
  enum {doAlphaVolumePostSkeleton   = disjunction<LocalOperators::doAlphaVolumePostSkeleton...>()};
  enum {doAlphaSkeleton             = disjunction<LocalOperators::doAlphaSkeleton...>()};
  enum {doAlphaBoundary             = disjunction<LocalOperators::doAlphaBoundary...>()};
  enum {doLambdaVolume              = disjunction<LocalOperators::doLambdaVolume...>()};
  enum {doLambdaVolumePostSkeleton  = disjunction<LocalOperators::doLambdaVolumePostSkeleton...>()};
  enum {doLambdaSkeleton            = disjunction<LocalOperators::doLambdaSkeleton...>()};
  enum {doLambdaBoundary            = disjunction<LocalOperators::doLambdaBoundary...>()};
  enum {doSkeletonTwoSided          = disjunction<LocalOperators::doSkeletonTwoSided...>()};

  // check that all of doSkeletonTwoSided have the same value
  static_assert(conjunction<LocalOperators::doSkeletonTwoSided...>()
           xor (not disjunction<LocalOperators::doSkeletonTwoSided...>()),
    "doSkeletonTwoSided should be the same for all operators");

  // assume all local operators can be cosntructed with same arguments
  template<class... Args>
  VariadicLocalOperator(Args&&... args)
  {
    Dune::Hybrid::forEach(Dune::range(_integral_size{}), [&](auto i) {
      auto& lop = std::get<i>(_lops);
      using LOP = std::decay_t<decltype(*lop)>;
      lop = std::make_unique<LOP>(args...);
    });
  }

  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume( const LFSU& lfsu, const LFSV& lfsv, LocalPattern& pattern) const
  {
    auto index_set = indices(lfsu,lfsv);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::patternVolume,lfsu,lfsv,pattern);
  }

  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume_post_skeleton
  ( const LFSU& lfsu, const LFSV& lfsv,
    LocalPattern& pattern) const
  {
    auto index_set = indices(lfsu,lfsv);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::patternVolumePostSkeleton,lfsu,lfsv,pattern);
  }

  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_skeleton
  ( const LFSU& lfsu_s, const LFSV& lfsv_s,
    const LFSU& lfsu_n, const LFSV& lfsv_n,
    LocalPattern& pattern_sn,
    LocalPattern& pattern_ns) const
  {
    auto index_set_s = indices(lfsu_s,lfsv_s);
    auto index_set_n = indices(lfsu_n,lfsv_n);

    auto predicate = [&](auto i) {
      auto inside_i = index_set_s.find(i) != index_set_s.end();
      auto outside_i = index_set_n.find(i) != index_set_n.end();

      if (inside_i xor outside_i)
        DUNE_THROW(MathError,
          "Variadic local operator cannot handle skeleton integrals when mappers have different indices");
      return (inside_i and outside_i);
    };

    applyLops(predicate,LocalOperatorApply::patternSkeleton,lfsu_s,lfsv_s,lfsu_n,lfsv_n,pattern_sn,pattern_ns);
  }

  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_boundary
  ( const LFSU& lfsu_s, const LFSV& lfsv_s,
    LocalPattern& pattern_ss) const
  {
    auto index_set = indices(lfsu_s,lfsv_s);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::patternBoundary,lfsu_s,lfsv_s,pattern_ss);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
            typename R>
  void alpha_volume
  ( const EG& eg,
    const LFSU& lfsu, const X& x, const LFSV& lfsv,
    R& r) const
  {
    auto index_set = indices(lfsu,lfsv);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::alphaVolume,eg,lfsu,x,lfsv,r);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
            typename R>
  void alpha_volume_post_skeleton
  ( const EG& eg,
    const LFSU& lfsu, const X& x, const LFSV& lfsv,
    R& r) const
  {
    auto index_set = indices(lfsu,lfsv);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::alphaVolumePostSkeleton,eg,lfsu,x,lfsv,r);
  }

  template<typename IG, typename LFSU, typename X, typename LFSV,
            typename R>
  void alpha_skeleton
  ( const IG& ig,
    const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
    const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
    R& r_s, R& r_n) const
  {
    auto index_set_s = indices(lfsu_s,lfsv_s);
    auto index_set_n = indices(lfsu_n,lfsv_n);

    auto predicate = [&](auto i) {
      auto inside_i = index_set_s.find(i) != index_set_s.end();
      auto outside_i = index_set_n.find(i) != index_set_n.end();

      if (inside_i xor outside_i)
        DUNE_THROW(MathError,
          "Variadic local operator cannot handle skeleton integrals when mappers have different indices");
      return (inside_i and outside_i);
    };

    applyLops(predicate,LocalOperatorApply::alphaSkeleton,ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,r_s,r_n);
  }

  template<typename IG, typename LFSU, typename X, typename LFSV,
            typename R>
  void alpha_boundary
  ( const IG& ig,
    const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
    R& r_s) const
  {
    auto index_set = indices(lfsu_s,lfsv_s);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::patternBoundary,ig,lfsu_s,x_s,lfsv_s,r_s);
  }

  template<typename EG, typename LFSV, typename R>
  void lambda_volume(const EG& eg, const LFSV& lfsv, R& r) const
  {
    auto index_set = indices(lfsv);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::lambdaVolume,eg,lfsv,r);
  }

  template<typename EG, typename LFSV, typename R>
  void lambda_volume_post_skeleton(const EG& eg, const LFSV& lfsv, R& r) const
  {
    auto index_set = indices(lfsv);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::lambdaVolumePostSkeleton,eg,lfsv,r);
  }

  template<typename IG, typename LFSV, typename R>
  void lambda_skeleton(const IG& ig,
                        const LFSV& lfsv_s, const LFSV& lfsv_n,
                        R& r_s, R& r_n) const
  {
    auto index_set_s = indices(lfsv_s);
    auto index_set_n = indices(lfsv_n);

    auto predicate = [&](auto i) {
      auto inside_i = index_set_s.find(i) != index_set_s.end();
      auto outside_i = index_set_n.find(i) != index_set_n.end();

      if (inside_i xor outside_i)
        DUNE_THROW(MathError,
          "Variadic local operator cannot handle skeleton integrals when mappers have different indices");
      return (inside_i and outside_i);
    };

    applyLops(predicate,LocalOperatorApply::lambdaSkeleton,ig,lfsv_s,lfsv_n,r_s,r_n);
  }

  template<typename IG, typename LFSV, typename R>
  void lambda_boundary(const IG& ig, const LFSV& lfsv_s, R& r_s) const
  {
    auto index_set = indices(lfsv_s);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::lambdaBoundary,ig,lfsv_s,r_s);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
            typename Y>
  void jacobian_apply_volume
  ( const EG& eg,
    const LFSU& lfsu, const X& z, const LFSV& lfsv,
    Y& y) const
  {
    auto index_set = indices(lfsv);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::jacobianApplyVolume,eg,lfsu,z,lfsv,y);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
            typename Y>
  void jacobian_apply_volume_post_skeleton
  ( const EG& eg,
    const LFSU& lfsu, const X& z, const LFSV& lfsv,
    Y& y) const
  {
    auto index_set = indices(lfsv);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::jacobianApplyVolumePostSkeleton,eg,lfsu,z,lfsv,y);
  }

  template<typename IG, typename LFSU, typename X, typename LFSV,
            typename Y>
  void jacobian_apply_skeleton
  ( const IG& ig,
    const LFSU& lfsu_s, const X& z_s, const LFSV& lfsv_s,
    const LFSU& lfsu_n, const X& z_n, const LFSV& lfsv_n,
    Y& y_s, Y& y_n) const
  {
    auto index_set_s = indices(lfsu_s,lfsv_s);
    auto index_set_n = indices(lfsu_n,lfsv_n);

    auto predicate = [&](auto i) {
      auto inside_i = index_set_s.find(i) != index_set_s.end();
      auto outside_i = index_set_n.find(i) != index_set_n.end();

      if (inside_i xor outside_i)
        DUNE_THROW(MathError,
          "Variadic local operator cannot handle skeleton integrals when mappers have different indices");
      return (inside_i and outside_i);
    };

    applyLops(predicate,LocalOperatorApply::jacobianApplySkeleton,ig,lfsu_s,z_s,lfsv_s,lfsu_n,z_n,lfsv_n,y_s,y_n);
  }

  template<typename IG, typename LFSU, typename X, typename LFSV,
            typename Y>
  void jacobian_apply_boundary
  ( const IG& ig,
    const LFSU& lfsu_s, const X& z_s, const LFSV& lfsv_s,
    Y& y_s) const
  {
    auto index_set = indices(lfsv_s);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::jacobianApplyBoundary,ig,lfsu_s,z_s,lfsv_s,y_s);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
            typename Y>
  void jacobian_apply_volume
  ( const EG& eg,
    const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv,
    Y& y) const
  {
    auto index_set = indices(lfsu,lfsv);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::jacobianApplyVolume,eg,lfsu,x,z,lfsv,y);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
            typename Y>
  void jacobian_apply_volume_post_skeleton
  ( const EG& eg,
    const LFSU& lfsu, const X& x, const X& z, const LFSV& lfsv,
    Y& y) const
  {
    auto index_set = indices(lfsu,lfsv);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::jacobianApplyVolumePostSkeleton,eg,lfsu,x,z,lfsv,y);
  }

  template<typename IG, typename LFSU, typename X, typename LFSV,
            typename Y>
  void jacobian_apply_skeleton
  ( const IG& ig,
    const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
    const LFSU& lfsu_n, const X& x_n, const X& z_n, const LFSV& lfsv_n,
    Y& y_s, Y& y_n) const
  {
    auto index_set_s = indices(lfsu_s,lfsv_s);
    auto index_set_n = indices(lfsu_n,lfsv_n);

    auto predicate = [&](auto i) {
      auto inside_i = index_set_s.find(i) != index_set_s.end();
      auto outside_i = index_set_n.find(i) != index_set_n.end();

      if (inside_i xor outside_i)
        DUNE_THROW(MathError,
          "Variadic local operator cannot handle skeleton integrals when mappers have different indices");
      return (inside_i and outside_i);
    };

    applyLops(predicate,LocalOperatorApply::jacobianApplySkeleton,ig,lfsu_s,x_s,z_s,lfsv_s,lfsu_n,x_n,z_n,lfsv_n,y_s,y_n);
  }

  template<typename IG, typename LFSU, typename X, typename LFSV,
            typename Y>
  void jacobian_apply_boundary
  ( const IG& ig,
    const LFSU& lfsu_s, const X& x_s, const X& z_s, const LFSV& lfsv_s,
    Y& y_s) const
  {
    auto index_set = indices(lfsu_s,lfsv_s);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::jacobianApplyBoundary,ig,lfsu_s,x_s,z_s,lfsv_s,y_s);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
            typename LocalMatrix>
  void jacobian_volume
  ( const EG& eg,
    const LFSU& lfsu, const X& x, const LFSV& lfsv,
    LocalMatrix& mat) const
  {
    auto index_set = indices(lfsu,lfsv);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::jacobianVolume,eg,lfsu,x,lfsv,mat);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
            typename LocalMatrix>
  void jacobian_volume_post_skeleton
  ( const EG& eg,
    const LFSU& lfsu, const X& x, const LFSV& lfsv,
    LocalMatrix& mat) const
  {
    auto index_set = indices(lfsu,lfsv);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::jacobianVolumePostSkeleton,eg,lfsu,x,lfsv,mat);
  }

  template<typename IG, typename LFSU, typename X, typename LFSV,
            typename LocalMatrix>
  void jacobian_skeleton
  ( const IG& ig,
    const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
    const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
    LocalMatrix& mat_ss, LocalMatrix& mat_sn,
    LocalMatrix& mat_ns, LocalMatrix& mat_nn) const
  {
    auto index_set_s = indices(lfsu_s,lfsv_s);
    auto index_set_n = indices(lfsu_n,lfsv_n);

    auto predicate = [&](auto i) {
      auto inside_i = index_set_s.find(i) != index_set_s.end();
      auto outside_i = index_set_n.find(i) != index_set_n.end();

      if (inside_i xor outside_i)
        DUNE_THROW(MathError,
          "Variadic local operator cannot handle skeleton integrals when mappers have different indices");
      return (inside_i and outside_i);
    };

    applyLops(predicate,LocalOperatorApply::jacobianSkeleton,ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,mat_ss,mat_sn,mat_ns,mat_nn);
  }

  template<typename IG, typename LFSU, typename X, typename LFSV,
            typename LocalMatrix>
  void jacobian_boundary
  ( const IG& ig,
    const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
    LocalMatrix& mat_ss) const
  {
    auto index_set = indices(lfsu_s,lfsv_s);

    auto predicate = [&](auto i) { return index_set.find(i) != index_set.end(); };

    applyLops(predicate,LocalOperatorApply::jacobianBoundary,ig,lfsu_s,x_s,lfsv_s,mat_ss);
  }

  template<class T>
  void setTime (T time)
  {
    Dune::Hybrid::forEach(Dune::range(_integral_size{}), [&](auto i) {
      const auto& lop = std::get<i>(_lops);
      lop->setTime(time);
    });
  }

  auto getTime () const
  {
    return std::get<0>(_lops)->getTime();
  }

  template<class T>
  void preStep (T time, T dt, int stages)
  {
    Dune::Hybrid::forEach(Dune::range(_integral_size{}), [&](auto i) {
      auto& lop = std::get<i>(_lops);
      lop->preStep(time,dt,stages);
    });
  }

  void postStep ()
  {
    Dune::Hybrid::forEach(Dune::range(_integral_size{}), [&](auto i) {
      auto& lop = std::get<i>(_lops);
      lop->postStep();
    });
  }

  template<class T>
  void preStage (T time, int r)
  {
    Dune::Hybrid::forEach(Dune::range(_integral_size{}), [&](auto i) {
      auto& lop = std::get<i>(_lops);
      lop->preStage(time,r);
    });
  }

  int getStage () const
  {
    return std::get<0>(_lops)->getStage();
  }

  void postStage ()
  {
    Dune::Hybrid::forEach(Dune::range(_integral_size{}), [&](auto i) {
      auto& lop = std::get<i>(_lops);
      lop->postStage();
    });
  }

  template<class T>
  T suggestTimestep (T dt) const
  {
    Dune::Hybrid::forEach(Dune::range(_integral_size{}), [&](auto i) {
      const auto& lop = std::get<i>(_lops);
      dt = std::min(dt,lop->suggestTimestep(dt));
    });
    return dt;
  }

  template<class... T>
  void update(T&&... args)
  {
    Dune::Hybrid::forEach(Dune::range(_integral_size{}), [&](auto i) {
      auto& lop = std::get<i>(_lops);
      lop->update(args...); // todo: check whether lop has update method
    });
  }

  template<class Entity>
  const auto& coefficient_mapper_inside(const Entity& entity_inside) const
  {
    // using CoefficientMapper = std::common_type_t<decltype(std::declval<LocalOperators>().coefficient_mapper_inside(std::declval<Entity>()))...>;
    // FIXME!
    return std::get<0>(_lops)->coefficient_mapper_inside(entity_inside);
  }

  template<class Entity>
  const auto& coefficient_mapper_outside(const Entity& entity_outside) const
  {
    // FIXME!
    return std::get<0>(_lops)->coefficient_mapper_outside(entity_outside);
  }

  template<class LFSU, class LFSV>
  const auto& lfs_components(const LFSU& lfsu, const LFSV& lfsv) const
  {
    _lfs_components.clear();
    auto index_set = indices(lfsu,lfsv);
    Dune::Hybrid::forEach(Dune::range(_integral_size{}), [&](auto i) {
      const auto& lop = std::get<i>(_lops);
      if (index_set.find(i) != index_set.end())
      {
        const auto& lfs_components = lop->lfs_components(lfsu,lfsv);
        _lfs_components.insert(_lfs_components.end(),lfs_components.begin(),lfs_components.end());
      }
    });
    return _lfs_components; //todo cache this
  }

private:
  using _integral_size = std::integral_constant<std::size_t,sizeof...(LocalOperators)>;
  mutable std::vector<std::size_t> _lfs_components;
  FiniteElementMapper _fe_mapper;
  std::tuple<std::unique_ptr<LocalOperators>...> _lops;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_VARIADIC_HH