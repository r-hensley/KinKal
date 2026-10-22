// Check endpoint-prior consistency in a single deterministic fit iteration.
#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/Detector/ParameterHit.hh"
#include "KinKal/Fit/Track.hh"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
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
      explicit ZeroEffectXing(double when, bool enabled=true) : when_(when), enabled_(enabled) {}
      void updateReference(PTRAJ const& trajectory, MetaIterConfig const&) override {
        reference_ = trajectory.nearestTraj(time());
      }
      void updateState(MetaIterConfig const&, bool) override {}
      Parameters params() const override { return Parameters(DVEC(), DMAT()); }
      double time() const override { return when_; }
      double transitTime() const override { return 0.01; }
      Trajectory const& referenceTrajectory() const override { return *reference_; }
      std::vector<MaterialXing> const& matXings() const override { return materials_; }
      void print(std::ostream& stream, int) const override { stream << "ZeroEffectXing"; }
      bool active() const override { return enabled_; }
    private:
      double when_;
      bool enabled_;
      KTRAJPTR reference_;
      std::vector<MaterialXing> materials_;
  };

  // Noncommuting positive-definite seed and hit covariances; all entries matter.
  void correlatedCovariances(DMAT& seed, DMAT& hit) {
    for(size_t row=0; row<NParams(); ++row) {
      for(size_t column=0; column<=row; ++column) {
        double ur=(row+1)/10.0, uc=(column+1)/10.0;
        double vr=(row%2 ? -1.0 : 1.0)*(row+1)/20.0;
        double vc=(column%2 ? -1.0 : 1.0)*(column+1)/20.0;
        seed(row,column)=(row==column ? 1.0 : 0.0)+ur*uc;
        hit(row,column)=(row==column ? 1.0+row/6.0 : 0.0)+vr*vc;
      }
    }
  }

  // Independent rational Gaussian elimination, evaluated before this native test.
  // C = inverse(2*s/4*inverse(Cseed) + 2/s*inverse(Chit)), s=(1+T)^2.
  // Mean uses those same weights with the seed and two measured parameter vectors.
  void correlatedExpected(double temperature, DVEC& mean, DMAT& covariance) {
    const double means[2][6] = {
      {0.80056298991701125, 0.20302416110985669, 0.002503085026879329, 1.0068708331287866, 0.005369036724964144, 100.01139901418667},
      {0.20109325616885013, 0.2036441872295004, 0.0037400868934346254, 1.0072883744590009, 0.0067338242284246605, 100.01093256168851},
    };
    const double covariances[2][6][6] = {
      {
        {0.40122261224486078, -0.00039655937550590214, 0.0038789579534614666, -0.00017569529243335161, 0.0067754023031404865, 0.00055674822774705963},
        {-0.00039655937550590214, 0.45712062569289513, -5.561771032144374e-05, 0.011717966536762369, 0.0015750397223019015, 0.018510299934292207},
        {0.0038789579534614666, -5.561771032144374e-05, 0.51257553074546414, 0.0020494765459948936, 0.022339595681283013, 0.0059448750032396773},
        {-0.00017569529243335161, 0.011717966536762369, 0.0020494765459948936, 0.57075883918182979, 0.0072048566305232759, 0.040438691173778664},
        {0.0067754023031404865, 0.0015750397223019015, 0.022339595681283013, 0.0072048566305232759, 0.6284288194274249, 0.016194346064220783},
        {0.00055674822774705963, 0.018510299934292207, 0.0059448750032396773, 0.040438691173778664, 0.016194346064220783, 0.69377064378042608},
      },
      {
        {0.40300100158635066, 0.0055556879582245866, 0.0093482435492152176, 0.011697879300573699, 0.015955665469883907, 0.018145052402667568},
        {0.0055556879582245866, 0.42452788253223245, 0.017817859841884415, 0.02628535767921426, 0.030947298251558068, 0.040202220985865768},
        {0.0093482435492152176, 0.017817859841884415, 0.45027844229820396, 0.037395229836142407, 0.049993467218627935, 0.057887545109634661},
        {0.011697879300573699, 0.02628535767921426, 0.037395229836142407, 0.4828178678844986, 0.064827113423738772, 0.083078897403339375},
        {0.015955665469883907, 0.030947298251558068, 0.049993467218627935, 0.064827113423738772, 0.52041643537413662, 0.10023183739464324},
        {0.018145052402667568, 0.040202220985865768, 0.057887545109634661, 0.083078897403339375, 0.10023183739464324, 0.56734629064524911},
      },
    };
    size_t index=temperature==0.0 ? 0u : 1u;
    for(size_t row=0; row<NParams(); ++row) {
      mean[row]=means[index][row];
      for(size_t column=0; column<=row; ++column)
        covariance(row,column)=covariances[index][row][column];
    }
  }

  bool checkFit(size_t boundaryMode, double temperature, bool correlated) {
    DVEC seedValues(0.0, 0.2, 0.0, 1.0, 0.0, 100.0);
    DMAT covariance;
    for(size_t parameter = 0; parameter < NParams(); ++parameter)
      covariance(parameter, parameter) = 1.0;
    DMAT hitCovariance = covariance;
    if(correlated) correlatedCovariances(covariance,hitCovariance);
    VEC3 field(0.0, 0.0, 0.0);
    Trajectory seed(Parameters(seedValues, covariance), 0.511, -1, field, TimeRange(0.0, 2.0));
    ParticleTrajectory<Trajectory> reference(seed);

    // Full-rank constraints avoid relying on weakly constrained geometry.
    // Only d0 differs from the seed, so its mean also detects double-counting.
    DVEC measuredValues = seedValues;
    measuredValues[Trajectory::d0_] = 1.0;
    Parameters measurement(measuredValues, hitCovariance);
    ParameterHit<Trajectory>::PMASK mask;
    mask.fill(true);
    Fit::HITCOL hits;
    hits.push_back(std::make_shared<ParameterHit<Trajectory>>(0.0, reference, measurement, mask));
    hits.push_back(std::make_shared<ParameterHit<Trajectory>>(2.0, reference, measurement, mask));
    Fit::EXINGCOL crossings;
    // 0: none, 1: one active, 2: two active, 3: one inactive crossing.
    if(boundaryMode==1) crossings.push_back(std::make_shared<ZeroEffectXing>(1.0));
    if(boundaryMode==2) {
      crossings.push_back(std::make_shared<ZeroEffectXing>(0.5));
      crossings.push_back(std::make_shared<ZeroEffectXing>(1.5));
    }
    if(boundaryMode==3) crossings.push_back(std::make_shared<ZeroEffectXing>(1.0,false));

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
    if(fit.fitTraj().pieces().size() != (boundaryMode==1 ? 2u : boundaryMode==2 ? 3u : 1u)) {
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
    DVEC expectedMean;
    DMAT expectedCovariance;
    if(correlated) {
      correlatedExpected(temperature,expectedMean,expectedCovariance);
    } else {
      for(size_t row=0; row<NParams(); ++row) {
        expectedMean[row]=expectedVariance *
          (2.0*priorPrecision*seedValues[row]+2.0*hitPrecision*measuredValues[row]);
        expectedCovariance(row,row)=expectedVariance;
      }
    }
    bool passed=true;
    size_t pieceIndex=0;
    for(auto const& piece : fit.fitTraj().pieces()) {
      auto const& parameters=piece->params();
      double maxError=0.0;
      std::cout << "STATE," << temperature << ',' << boundaryMode << ',' << correlated << ',' << pieceIndex;
      for(size_t row=0; row<NParams(); ++row) {
        double actual=parameters.parameters()[row];
        if(!std::isfinite(actual)) return false;
        maxError=std::max(maxError,std::abs(actual-expectedMean[row])/(1.0+std::abs(expectedMean[row])));
        std::cout << ',' << actual;
      }
      for(size_t row=0; row<NParams(); ++row) {
        for(size_t column=0; column<NParams(); ++column) {
          double actual=parameters.covariance()(row,column), expected=expectedCovariance(row,column);
          if(!std::isfinite(actual)) return false;
          maxError=std::max(maxError,std::abs(actual-expected)/(1.0+std::abs(expected)));
          std::cout << ',' << actual;
        }
      }
      std::cout << '\n';
      bool piecePassed=maxError<=1.0e-10;
      std::cout << "CHECK," << temperature << ',' << boundaryMode << ',' << correlated
                << ',' << pieceIndex++ << ',' << piecePassed << ',' << maxError << '\n';
      passed=passed && piecePassed;
    }
    return passed;
  }
}

int main() {
  std::cout << std::setprecision(17);
  size_t cases=0, failures=0;
  for(bool correlated : {false,true}) {
    for(double temperature : {0.0,1.0}) {
      for(size_t mode : {0u,1u,2u,3u}) {
        ++cases;
        if(!checkFit(mode,temperature,correlated)) ++failures;
      }
    }
  }
  std::cout << "SUMMARY," << cases << ',' << failures << '\n';
  return failures==0 ? 0 : 1;
}
