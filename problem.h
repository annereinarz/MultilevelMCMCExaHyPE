#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/SamplingAlgorithms/AMProposal.h"
#include "MUQ/SamplingAlgorithms/MHProposal.h"
#include "MUQ/SamplingAlgorithms/CrankNicolsonProposal.h"
#include "MUQ/SamplingAlgorithms/ParallelizableMIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SubsamplingMIProposal.h"
#include <chrono>
#include <iomanip>

using namespace std::chrono;


//TODO: read in
const int NUM_PARAM = 2;
int count;
const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, 0, ", ", ";\n", "", "", "", ";\n");

int saved_argc; char** saved_argv;

class MySamplingProblem : public AbstractSamplingProblem {
public:
  MySamplingProblem(std::shared_ptr<parcer::Communicator> comm,std::shared_ptr<parcer::Communicator> globalComm,std::shared_ptr<MultiIndex> index)
  : AbstractSamplingProblem(Eigen::VectorXi::Constant(1,NUM_PARAM), Eigen::VectorXi::Constant(1,NUM_PARAM))
{
    this->comm = comm;
    this->globalComm = globalComm;
    this->index = index;
    muq::setCommunicator(comm->GetMPICommunicator());
    muq::init(saved_argc,saved_argv);
}

  virtual ~MySamplingProblem(){
  }

  virtual double LogDensity(std::shared_ptr<SamplingState> const& state) override {
    lastState = state;

    int subgroup_consistent_rank = globalComm->GetRank();
    comm->Bcast(subgroup_consistent_rank, 0);

    std::string inputsfilename = "/tmp/inputs" + std::to_string(subgroup_consistent_rank) + ".txt";
    std::string outputsfilename = "/tmp/outputs" + std::to_string(subgroup_consistent_rank) + ".txt";

    std::ofstream inputsfile(inputsfilename);
    typedef std::numeric_limits<double> dl;
    inputsfile << std::fixed << std::setprecision(dl::digits10);
    for (int i = 0; i < state->state[0].rows(); i++) {
      inputsfile << state->state[0](i) << std::endl;
    }
    inputsfile.close();

    auto start = high_resolution_clock::now();
    // TODO: Need level dependence!
    int status;
    std::string inoutsuffix = inputsfilename + " " + outputsfilename;
    if(index->GetValue(0) == 0) {
      status = system((std::string("cd /hppfs/work/pr83no/ge68wax4/MUQ/ExaHyPE-Tsunami/ApplicationExamples/SWE/SWE_asagi_limited_l0 && ./ExaHyPE-SWE ../SWE_asagi_limited_l0.exahype2 ") + inoutsuffix).c_str());
    } else if(index->GetValue(0) == 1) {
      status = system((std::string("cd /hppfs/work/pr83no/ge68wax4/MUQ/ExaHyPE-Tsunami/ApplicationExamples/SWE/SWE_asagi_limited_l1 && ./ExaHyPE-SWE ../SWE_asagi_limited_l1.exahype2 ") + inoutsuffix).c_str());
    } else if(index->GetValue(0) == 2) {
      status = system((std::string("cd /hppfs/work/pr83no/ge68wax4/MUQ/ExaHyPE-Tsunami/ApplicationExamples/SWE/SWE_asagi_limited_l2 && ./ExaHyPE-SWE ../SWE_asagi_limited_l2.exahype2 ") + inoutsuffix).c_str());
    } else {
      std::cerr << "Unknown model requested by client!" << std::endl;
      exit(-1);
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    std::cout << "Exahype exit status " << status << std::endl;
    std::cout << "Sample took " << duration.count() << std::endl;


    std::vector<double> output(4, 0.0);
    std::ifstream outputsfile(outputsfilename);
    for (int i = 0; i < output.size(); i++) {
      outputsfile >> output[i];
    }
    outputsfile.close();

    std::cout << "Read outputs from exahype:" << std::endl;
    for (int i = 0; i < output.size(); i++) {
      std::cout << output[i] << std::endl;
    }

    double sigma = 1.0;
    double likelihood = 1234;
    if (comm->GetRank() == 0) // Only controller rank does likelihood calculation (muq expects this, so it's fine)
      likelihood = calculateLikelihood(output,globalComm->GetRank(), index->GetValue(0));
    return likelihood;
  };

  virtual std::shared_ptr<SamplingState> QOI() override {
    assert (lastState != nullptr);
    return std::make_shared<SamplingState>(lastState->state, 1.0);
  }

private:
  std::shared_ptr<SamplingState> lastState = nullptr;
  std::shared_ptr<parcer::Communicator> comm;
  std::shared_ptr<parcer::Communicator> globalComm;
  std::shared_ptr<MultiIndex> index;
};


class MyInterpolation : public MIInterpolation {
public:
  std::shared_ptr<SamplingState> Interpolate (std::shared_ptr<SamplingState> const& coarseProposal, std::shared_ptr<SamplingState> const& fineProposal) override {
    return std::make_shared<SamplingState>(coarseProposal->state);
  }
};

class MyMIComponentFactory : public ParallelizableMIComponentFactory {
public:
  MyMIComponentFactory(pt::ptree pt, std::shared_ptr<parcer::Communicator> globalComm) : _pt(pt), _globalComm(globalComm) {}

  virtual std::shared_ptr<MCMCProposal> Proposal (std::shared_ptr<MultiIndex> const& index, std::shared_ptr<AbstractSamplingProblem> const& samplingProblem) override {
    pt::ptree pt;
    pt.put("BlockIndex",0);

    /*Eigen::VectorXd mu = Eigen::VectorXd::Zero(NUM_PARAM);
    Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(NUM_PARAM,NUM_PARAM);

    auto prior = std::make_shared<Gaussian>(mu, cov);

    return std::make_shared<CrankNicolsonProposal>(pt, samplingProblem, prior);*/

    pt.put("AdaptStart", 100);
    pt.put("AdaptSteps", 100);
    pt.put("InitialVariance", 10);
    return std::make_shared<AMProposal>(pt, samplingProblem);

    /*auto mu = Eigen::VectorXd::Zero(NUM_PARAM);
    Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(NUM_PARAM, NUM_PARAM);
    cov *= 0.5;

    auto prior = std::make_shared<Gaussian>(mu, cov);

    return std::make_shared<MHProposal>(pt, samplingProblem, prior);*/
  }

  void SetComm(std::shared_ptr<parcer::Communicator> const& comm) override {
    _comm = comm;
  }

  virtual std::shared_ptr<MultiIndex> FinestIndex() override {
    auto index = std::make_shared<MultiIndex>(1);
    index->SetValue(0, 1);
    return index;
  }

  virtual std::shared_ptr<MCMCProposal> CoarseProposal (std::shared_ptr<MultiIndex> const& fineIndex,
                                                        std::shared_ptr<MultiIndex> const& coarseIndex,
                                                        std::shared_ptr<AbstractSamplingProblem> const& coarseProblem,
                                                        std::shared_ptr<SingleChainMCMC> const& coarseChain) override {
    pt::ptree ptProposal = _pt;
    ptProposal.put("BlockIndex",0);
    return std::make_shared<SubsamplingMIProposal>(ptProposal, coarseProblem, coarseIndex, coarseChain);
  }

  virtual std::shared_ptr<AbstractSamplingProblem> SamplingProblem (std::shared_ptr<MultiIndex> const& index) override {

    return std::make_shared<MySamplingProblem>(_comm,_globalComm,index);
  }

  virtual std::shared_ptr<MIInterpolation> Interpolation (std::shared_ptr<MultiIndex> const& index) override {
    return std::make_shared<MyInterpolation>();
  }

  virtual Eigen::VectorXd StartingPoint (std::shared_ptr<MultiIndex> const& index) override {
    //Starting guess: zero
    Eigen::VectorXd start = Eigen::VectorXd::Ones(NUM_PARAM);
    start(0) = 7.0;
    start(1) = 7.0;
    return start;
  }

private:
  std::shared_ptr<parcer::Communicator> _comm;
  std::shared_ptr<parcer::Communicator> _globalComm;
  pt::ptree _pt;

};
