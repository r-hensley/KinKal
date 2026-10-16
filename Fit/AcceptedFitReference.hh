#ifndef KinKal_AcceptedFitReference_hh
#define KinKal_AcceptedFitReference_hh

namespace KinKal {
  template<class KTRAJ> class ParticleTrajectory;

  // Identifies accepted-posterior attachments without magnetic-domain transport.
  // Construction and domain replacement also call updateReference, but do not
  // signify that a new fit has assimilated the current measurement weights.
  // Nested scopes and exceptions restore the prior context; threads are independent.
  // Consumers still need to match the active hit, cached weight and supported time.
  template<class KTRAJ> class AcceptedFitReference {
    public:
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      explicit AcceptedFitReference(PTRAJ const& trajectory, bool hasMagneticDomains) noexcept : previous_(current_) {
        current_ = hasMagneticDomains ? nullptr : &trajectory;
      }
      ~AcceptedFitReference() { current_ = previous_; }
      AcceptedFitReference(AcceptedFitReference const&) = delete;
      AcceptedFitReference& operator=(AcceptedFitReference const&) = delete;

      static bool matches(PTRAJ const& trajectory) noexcept { return current_ == &trajectory; }
    private:
      inline static thread_local PTRAJ const* current_ = nullptr;
      PTRAJ const* const previous_;
  };
}
#endif
