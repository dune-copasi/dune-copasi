#ifndef DUNE_COPASI_OPERATOR_FORWARD_INSTATIONARY_ASSEMBLER_JACOBIAN_HH
#define DUNE_COPASI_OPERATOR_FORWARD_INSTATIONARY_ASSEMBLER_JACOBIAN_HH

#include <dune/copasi/operator/forward/instationary/coefficients.hh>
#include <dune/copasi/operator/forward/instationary/traits.hh>
#include <dune/copasi/operator/operator.hh>

#include <dune/copasi/operator/common/matrixIterator.hh>
#include <dune/copasi/operator/common/indexIterator.hh>
#include <dune/copasi/operator/common/vectorEntry.hh>
#include <dune/copasi/operator/common/matrixEntry.hh>
#include <dune/pdelab/basis/backend/istl.hh> 

// we retain the pdelab interface for compatibility
#include <dune/pdelab/operator/local_assembly/interface.hh>

#include <dune/pdelab/common/local_container.hh>
#include <dune/pdelab/common/local_matrix.hh>
#include <dune/pdelab/common/tree_traversal.hh>
#include <dune/pdelab/common/for_each_entity.hh>
#include <dune/pdelab/common/execution.hh>
#include <dune/pdelab/common/algebra.hh>
#include <dune/pdelab/common/container_entry.hh>

#include <dune/pdelab/concepts/basis.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/common/float_cmp.hh>

#include <array>
#include <functional>
#include <vector>

namespace Dune::Copasi::inline Experimental {

template<class          Coefficients,
         class          Residual,
         PDELab::Concept::Basis TrialBasis,
         PDELab::Concept::Basis TestBasis,
         class          MassLocalOperator,
         class          StiffnessLocalOperator,
         class          JacobianContainer,
         class TimeQuantity = double,
         class DurationQuantity  = double,
         DurationPosition dt_position = DurationPosition::StiffnessNumerator>
class InstationaryJacobianAssembler
  : public Operator<Coefficients, Residual>
{
  using StageCoefficients = typename std::ranges::range_value_t<Coefficients>;
  using StageResidual     = typename std::ranges::range_value_t<Residual>;
  using StageJacobian     = typename std::ranges::range_value_t<typename std::ranges::range_value_t<JacobianContainer>>;

  static_assert(PDELab::Concept::Container<StageCoefficients,TrialBasis>);
  static_assert(PDELab::Concept::Container<StageResidual,    TestBasis>);

  static_assert(std::is_copy_constructible_v<MassLocalOperator>);
  static_assert(std::is_copy_constructible_v<StiffnessLocalOperator>);

  using LocalTestBasis = typename TestBasis::LocalView;
  using LocalTrialBasis = typename TrialBasis::LocalView;

  using IndicesBackend = PDELab::ISTLUniformBackend<int>;
  using IndicesVector = typename TestBasis::template Container<IndicesBackend>;

  // stage jacobian is required to be indexable by local test and trial basis degrees of freedom (see localContainerEntry)
  using LocalJacobian = PDELab::LocalMatrixBuffer<TestBasis, TrialBasis, StageJacobian>;

  using MassFactor = typename InstationaryTraits<
    dt_position>::template MassFactor<DurationQuantity>;
  using StiffnessFactor = typename InstationaryTraits<
    dt_position>::template StiffnessFactor<DurationQuantity>;

  using LocalCoefficients       = PDELab::LocalContainerBuffer<TrialBasis, const StageCoefficients>;
  using LocalMassJacobian       = PDELab::WeightedLocalMatrixView<LocalJacobian, MassFactor>;
  using LocalStiffnessJacobian  = PDELab::WeightedLocalMatrixView<LocalJacobian, StiffnessFactor>;

  using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<typename TrialBasis::EntitySet>;

public:
  InstationaryJacobianAssembler(const TrialBasis& trial,
                                const TestBasis& test,
                                const MassLocalOperator& mass_lop,
                                const StiffnessLocalOperator& stiff_lop)
    : _mass_lop{ mass_lop }
    , _stiff_lop{ stiff_lop }
    , _trial{ trial }
    , _test{ test }
    , _mapper{ _trial.entitySet(), Dune::mcmgElementLayout() }
    , _constrainDOF{ _test.makeContainer(IndicesBackend{}) }
  {
    PDELab::PropertyTree& properties = *this;

    properties["container"].documentation = "Matrix (or tensor) container holding the jacobian of an instationary operator";
    properties["container"] = std::make_shared<JacobianContainer>();

    initConstraints();

  }

  PDELab::ErrorCondition apply(const Coefficients& coefficients, Residual& residual) override
  {
    TRACE_EVENT("dune", "Instationary::Jacobian");
    const JacobianContainer& jac = getJacobianContainer();
    // PDELab::linearTransformation(PDELab::Execution::par, coefficients, jac, residual);
    jac[0][0].mv(coefficients[0], residual[0]);
    return {};
  }

  void linearize(const Coefficients& linearization_point) {
    using namespace Indices;
    const InstationaryCoefficients& icoeff = getInstationaryCoefficients();
    // help compiler to unroll short loops
    if (icoeff.extent(0) == 1 && icoeff.extent(1) == 1)
      linearizeImpl(linearization_point, _1, _1);
    else if (icoeff.extent(0) == 2 && icoeff.extent(1) == 2)
      linearizeImpl(linearization_point, _2, _2);
    else
      linearizeImpl(linearization_point, icoeff.extent(0), icoeff.extent(1));
  }

private:

  void linearizeImpl(const Coefficients& linearization_point, auto steps, auto stages)
  {
    TRACE_EVENT("dune", "Instationary::JacobianLinearize");
    LocalTrialBasis ltrial_in  = _trial.localView();
    LocalTrialBasis ltrial_out = _trial.localView();

    LocalTestBasis ltest_in  = _test.localView();
    LocalTestBasis ltest_out = _test.localView();

  //   assert(this->tableau().stages() == stages);
  //   assert(linearization_point.size() == stages);

    JacobianContainer& jac = getJacobianContainer();
    const InstationaryCoefficients& icoeff = getInstationaryCoefficients();
    const TimeQuantity& time = getTime();
    const DurationQuantity& duration = getDuration();

    PDELab::forEachContainerEntry(PDELab::Execution::par_unseq, jac, []<class T>(T& v){v = T{0.};});

    std::vector<std::vector<bool>> do_mass(stages);
    std::vector<std::vector<bool>> do_stiff(stages);
    for (std::size_t stage = 0; stage != stages; ++stage) {
      do_mass[stage].resize(steps);
      do_stiff[stage].resize(steps);
      for (std::size_t step = 0; step != steps; ++step) {
        do_mass[stage][step] = icoeff.doMass(stage, step);
        do_stiff[stage][step] = icoeff.doStiffness(stage, step);
      }
    }

    const auto& it = _trial.entitySet().template begin<0>();
    if (_trial.entitySet().size(0) == 0)
      return;

    bind(*it, ltest_in, ltrial_in);
    bind(*it, ltest_out, ltrial_out);

    std::vector<LocalCoefficients> llin_in, llin_out;
    for (std::size_t step = 0; step != steps; ++step) {
      llin_in.emplace_back(LocalCoefficients{ ltrial_in, &linearization_point[step] });
      llin_out.emplace_back(LocalCoefficients{ ltrial_in, &linearization_point[step] });
    }

    LocalJacobian ljac_ini{ _test, _trial, jac[0][0] };
    DynamicMatrix<LocalJacobian> ljac_ii(stages, steps, ljac_ini), ljac_io(stages, steps, ljac_ini), ljac_oi(stages, steps, ljac_ini), ljac_oo(stages, steps, ljac_ini);
    for (std::size_t stage = 0; stage != stages; ++stage) {
      for (std::size_t step = 0; step != steps; ++step) {
        ljac_ii[stage][step] = LocalJacobian{ ltest_in,  ltrial_in,  jac[stage][step] };
        ljac_io[stage][step] = LocalJacobian{ ltest_in,  ltrial_out, jac[stage][step] };
        ljac_oi[stage][step] = LocalJacobian{ ltest_out, ltrial_in,  jac[stage][step] };
        ljac_oo[stage][step] = LocalJacobian{ ltest_out, ltrial_out, jac[stage][step] };
      }
    }

    unbind(ltest_in, ltrial_in);
    unbind(ltest_out, ltrial_out);

    auto sub_time_step = [&](std::size_t step) {
      return time + duration * icoeff.timeWeight(step);
    };

    auto mass_weight = [&](std::size_t stage, std::size_t step) {
      return icoeff.massWeight(stage, step) * InstationaryTraits<dt_position>::massFactor(duration);
    };

    auto stiff_weight = [&](std::size_t stage, std::size_t step) {
      return icoeff.stiffnessWeight(stage, step) * InstationaryTraits<dt_position>::stiffnessFactor(duration);
    };

    const auto is_linear = (not PDELab::LocalAssembly::isLinear(_stiff_lop) || not PDELab::LocalAssembly::isLinear(_mass_lop));
    const auto& es = _trial.entitySet();
    const auto& esp = _trial.entitySetPartition();

    PDELab::forEachEntity(PDELab::LocalAssembly::executionPolicy(_stiff_lop), esp,
      [ this, &do_stiff, &do_mass, &stages, &steps, &es,
        ltest_in, ltrial_in, ltest_out, ltrial_out,
        llin_in, llin_out,
        ljac_ii, ljac_io, ljac_oi, ljac_oo,
        mlop = _mass_lop, slop = _stiff_lop,
        is_linear, sub_time_step, mass_weight, stiff_weight
      ] (const auto& entity) mutable
    {
      if (PDELab::LocalAssembly::skipEntity(mlop, entity)) {
        if (not PDELab::LocalAssembly::skipEntity(slop, entity))
          DUNE_THROW(InvalidStateException, "skip methods should yield the same result");
        return;
      }

      bind(entity, ltest_in, ltrial_in);

      if (is_linear)
        for (std::size_t step = 0; step != steps; ++step)
          llin_in[step].load(ltrial_in);

      if (PDELab::LocalAssembly::doVolume(slop) || PDELab::LocalAssembly::doVolume(mlop)) {
        for (std::size_t step = 0; step != steps; ++step) {
          TimeQuantity tp = sub_time_step(step);
          for (std::size_t stage = 0; stage != stages; ++stage) {
            ljac_ii[stage][step].clear(ltest_in, ltrial_in);
            if (do_mass[stage][step]) {
              LocalMassJacobian lmass_jac_ii{ljac_ii[stage][step], mass_weight(stage, step)};
              PDELab::LocalAssembly::jacobianVolume(mlop, tp, ltrial_in, llin_in[step], ltest_in, lmass_jac_ii);
            }
            if (do_stiff[stage][step]) {
              LocalStiffnessJacobian lstiff_jac_ii{ljac_ii[stage][step], stiff_weight(stage, step)};
              PDELab::LocalAssembly::jacobianVolume(slop, tp, ltrial_in, llin_in[step], ltest_in, lstiff_jac_ii);
            }
          }
        }
      }

      if (  PDELab::LocalAssembly::doSkeleton(slop) || PDELab::LocalAssembly::doBoundary(slop)
          || PDELab::LocalAssembly::doSkeleton(mlop) || PDELab::LocalAssembly::doBoundary(mlop)) {

        auto id_in = _mapper.index(entity);

        for (const auto& is : intersections(es, entity)) {
          if (PDELab::LocalAssembly::skipIntersection(mlop, is)) {
            if (not PDELab::LocalAssembly::skipIntersection(slop, is))
              DUNE_THROW(InvalidStateException, "skip methods should yield the same result");
            continue;
          }
          if (is.neighbor() and (PDELab::LocalAssembly::doSkeleton(slop) || PDELab::LocalAssembly::doSkeleton(mlop))) { // interior and periodic cases

            const auto& entity_out = is.outside();
            auto id_out = _mapper.index(entity_out);

            // visit face only the first time we see it
            bool skip = id_in < id_out;

            // we have to be sure that outside entity was not skipped,
            // otherwise, neither side will be visited
            if (skip and not PDELab::LocalAssembly::skipEntity(mlop, entity_out))
              continue;

            bind(entity_out, ltest_out, ltrial_out);

            for (std::size_t step = 0; step != steps; ++step)
              if (is_linear)
                llin_out[step].load(ltrial_out);

            for (std::size_t stage = 0; stage != stages; ++stage) {
              for (std::size_t step = 0; step != steps; ++step) {
                TimeQuantity tp = sub_time_step(step);
                ljac_io[stage][step].clear(ltest_in,  ltrial_out);
                ljac_oi[stage][step].clear(ltest_out, ltrial_in);
                ljac_oo[stage][step].clear(ltest_out, ltrial_out);

                if (do_mass[stage][step]) {
                  LocalMassJacobian lmass_jac_ii{ljac_ii[stage][step], mass_weight(stage, step)};
                  LocalMassJacobian lmass_jac_io{ljac_io[stage][step], mass_weight(stage, step)};
                  LocalMassJacobian lmass_jac_oi{ljac_oi[stage][step], mass_weight(stage, step)};
                  LocalMassJacobian lmass_jac_oo{ljac_oo[stage][step], mass_weight(stage, step)};

                  PDELab::LocalAssembly::jacobianSkeleton(mlop, is, tp, ltrial_in, llin_in[step], ltest_in, ltrial_out, llin_out[step], ltest_out, lmass_jac_ii, lmass_jac_io, lmass_jac_oi, lmass_jac_oo);
                }

                if (do_stiff[stage][step]) {
                  LocalStiffnessJacobian lstiff_jac_ii{ljac_ii[stage][step], stiff_weight(stage, step)};
                  LocalStiffnessJacobian lstiff_jac_io{ljac_io[stage][step], stiff_weight(stage, step)};
                  LocalStiffnessJacobian lstiff_jac_oi{ljac_oi[stage][step], stiff_weight(stage, step)};
                  LocalStiffnessJacobian lstiff_jac_oo{ljac_oo[stage][step], stiff_weight(stage, step)};

                  PDELab::LocalAssembly::jacobianSkeleton(slop, is, tp, ltrial_in, llin_in[step], ltest_in, ltrial_out, llin_out[step], ltest_out, lstiff_jac_ii, lstiff_jac_io, lstiff_jac_oi, lstiff_jac_oo);
                }
                ljac_io[stage][step].fetch_add(ltest_in, ltrial_out);
                ljac_oi[stage][step].fetch_add(ltest_out, ltrial_in);
                ljac_oo[stage][step].fetch_add(ltest_out, ltrial_out);
              }
            }
            unbind(ltest_out, ltrial_out);
          } else if (is.boundary() && (PDELab::LocalAssembly::doBoundary(slop) || PDELab::LocalAssembly::doBoundary(mlop))) {
            for (std::size_t stage = 0; stage != stages; ++stage) {
              for (std::size_t step = 0; step != steps; ++step) {
                TimeQuantity tp = sub_time_step(step);
                if (do_mass[stage][step]) {
                  LocalMassJacobian lmass_jac_ii{ljac_ii[stage][step], mass_weight(stage, step)};
                  PDELab::LocalAssembly::jacobianBoundary(mlop, is, tp, ltrial_in, llin_in[step], ltest_in, lmass_jac_ii);
                }

                if (do_stiff[stage][step]) {
                  LocalStiffnessJacobian lstiff_jac_ii{ljac_ii[stage][step], stiff_weight(stage, step)};
                  PDELab::LocalAssembly::jacobianBoundary(slop, is, tp, ltrial_in, llin_in[step], ltest_in, lstiff_jac_ii);
                }
              }
            }
          }

        }
      }

      for (std::size_t stage = 0; stage != stages; ++stage)
        for (std::size_t step = 0; step != steps; ++step)
          ljac_ii[stage][step].fetch_add(ltest_in, ltrial_in);

      unbind(ltest_in, ltrial_in);
    });

    // compile time decission to apply constraints
    if constexpr( _mass_lop.localAssembleDoConstraints() or  _stiff_lop.localAssembleDoConstraints() )
    {
      // make constrained DOF unit diagional in Jacobian container
      applyConstraints();
    }
  }

  /**
   * @brief Initialize which DOFs have to be constrained
   */
  void initConstraints()
  {
    if constexpr( _mass_lop.localAssembleDoConstraints() )
    {
      _mass_lop.localAssembleConstraints( _test, _constrainDOF );  
    }
    if constexpr( _stiff_lop.localAssembleDoConstraints() )
    {
      _stiff_lop.localAssembleConstraints( _test, _constrainDOF); 
    }
  }

  /**
   * @brief Apply constraints to the constrained DOFs
   *  
   * @details This function will change the Jacobian row vectors associated with the 
   *          constrained degrees of freedom. A row iterator will be used to change 
   *          all off-diagonal values to zero and put a 1 on the diagonal. 
   */
  void applyConstraints()
  {

    JacobianContainer& jac = getJacobianContainer();
    
    using IndexIterator = typename Dune::Copasi::index_iterator<IndicesVector>;
    
    for (IndexIterator it(_constrainDOF); not it.at_end(); ++it)
    {
      if( (*it) == 1 )
      {
        // iterator over the row and change all values to zero 
        for(Dune::Copasi::matrix_row_iterator it_jac(jac[0][0], it.index()); not it_jac.at_end(); ++it_jac)
        {
          *it_jac = 0.0;
        }
        // add a 1.0 to the diagional
        Dune::Copasi::matrixEntry( jac[0][0], it.index() , it.index() ) = double{1.0};

      }
    }
  }

  JacobianContainer& getJacobianContainer() {
    return this->template get<JacobianContainer>("container");
  }

  const InstationaryCoefficients& getInstationaryCoefficients() const {
    return this->template get<InstationaryCoefficients>("instationary_coefficients");
  }

  const TimeQuantity& getTime() const {
    return this->template get<TimeQuantity>("time");
  }

  const DurationQuantity& getDuration() const {
    return this->template get<DurationQuantity>("duration");
  }


  MassLocalOperator _mass_lop;
  StiffnessLocalOperator _stiff_lop;
  TrialBasis _trial;
  TestBasis _test;
  Mapper _mapper;
  IndicesVector _constrainDOF;
};

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_COPASI_OPERATOR_FORWARD_INSTATIONARY_ASSEMBLER_JACOBIAN_HH
