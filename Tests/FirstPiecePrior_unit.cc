// Check endpoint-prior consistency in a single deterministic fit iteration.
#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/Detector/ParameterHit.hh"
#include "KinKal/Fit/Track.hh"

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace {
  using namespace KinKal;
  using Trajectory = KinematicLine;
  using Fit = Track<Trajectory>;

  // Force a smoothed piece boundary without changing the underlying state.
  // Its two sides must therefore have the same parameters and covariance.
  class ZeroEffectXing : public ElementXing<Trajectory> {
    public:
      void updateReference(PTRAJ const& trajectory) override {
        reference_ = trajectory.nearestTraj(time());
      }
      void updateState(MetaIterConfig const&, bool) override {}
      Parameters params() const override { return Parameters(DVEC(), DMAT()); }
      double time() const override { return 1.0; }
      double transitTime() const override { return 0.01; }
      Trajectory const& referenceTrajectory() const override { return *reference_; }
      std::vector<MaterialXing> const& matXings() const override { return materials_; }
      void print(std::ostream& stream, int) const override { stream << "ZeroEffectXing"; }
      bool active() const override { return true; }
    private:
      KTRAJPTR reference_;
      std::vector<MaterialXing> materials_;
  };

  bool close(double actual, double expected) {
    return std::isfinite(actual) && std::abs(actual - expected) <= 1.0e-10 * (1.0 + std::abs(expected));
  }

  bool checkFit(bool withBoundary, double temperature) {
    DVEC seedValues(0.0, 0.2, 0.0, 1.0, 0.0, 100.0);
    DMAT covariance;
    for(size_t parameter = 0; parameter < NParams(); ++parameter)
      covariance(parameter, parameter) = 1.0;
    VEC3 field(0.0, 0.0, 0.0);
    Trajectory seed(Parameters(seedValues, covariance), 0.511, -1, field, TimeRange(0.0, 2.0));
    ParticleTrajectory<Trajectory> reference(seed);

    // Full-rank constraints avoid relying on weakly constrained geometry.
    // Only d0 differs from the seed, so its mean also detects double-counting.
    DVEC measuredValues = seedValues;
    measuredValues[Trajectory::d0_] = 1.0;
    Parameters measurement(measuredValues, covariance);
    ParameterHit<Trajectory>::PMASK mask;
    mask.fill(true);
    Fit::HITCOL hits;
    hits.push_back(std::make_shared<ParameterHit<Trajectory>>(0.0, reference, measurement, mask));
    hits.push_back(std::make_shared<ParameterHit<Trajectory>>(2.0, reference, measurement, mask));
    Fit::EXINGCOL crossings;
    if(withBoundary) crossings.push_back(std::make_shared<ZeroEffectXing>());

    Config config;
    config.schedule_.emplace_back(temperature);
    config.maxniter_ = 0; // canIterate compares this with the zero-based iteration index.
    config.dwt_ = 4.0;
    config.bfcorr_ = false;
    config.ends_ = false;
    UniformBFieldMap fieldMap(field);
    Fit fit(config, fieldMap, seed, hits, crossings);
    if(!fit.hasTraj() || !fit.fitStatus().usable()) {
      std::cerr << "Fit failed: " << fit.fitStatus() << '\n';
      return false;
    }
    if(fit.fitStatus().iter_ != 0) {
      std::cerr << "Expected exactly one fit iteration\n";
      return false;
    }
    if(fit.fitTraj().pieces().size() != (withBoundary ? 2u : 1u)) {
      std::cerr << "Unexpected number of fitted pieces\n";
      return false;
    }

    // Independent scalar normal equations: two unit-variance measurements
    // and the two deweighted endpoint regularizers. At temperature zero,
    // d0 must be 0.8 with variance 0.4. Omitting one prior gives 8/9 and 4/9.
    double varianceScale = (1.0 + temperature) * (1.0 + temperature);
    double priorPrecision = varianceScale / config.dwt_;
    double hitPrecision = 1.0 / varianceScale;
    double expectedVariance = 1.0 / (2.0 * priorPrecision + 2.0 * hitPrecision);
    for(auto const& piece : fit.fitTraj().pieces()) {
      auto const& parameters = piece->params();
      for(size_t row = 0; row < NParams(); ++row) {
        double expectedMean = expectedVariance *
          (2.0 * priorPrecision * seedValues[row] + 2.0 * hitPrecision * measuredValues[row]);
        if(!close(parameters.parameters()[row], expectedMean)) {
          std::cerr << "Incorrect fitted mean for parameter " << row << '\n';
          return false;
        }
        for(size_t column = 0; column < NParams(); ++column) {
          double expected = row == column ? expectedVariance : 0.0;
          if(!close(parameters.covariance()(row, column), expected)) {
            std::cerr << "Incorrect covariance at " << row << ',' << column << '\n';
            return false;
          }
        }
      }
    }
    return true;
  }
}

int main() {
  for(double temperature : {0.0, 1.0}) {
    for(bool withBoundary : {false, true}) {
      if(!checkFit(withBoundary, temperature)) {
        std::cerr << "Failed with temperature=" << temperature
                  << ", material boundary=" << withBoundary << '\n';
        return 1;
      }
    }
  }
  return 0;
}
