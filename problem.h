#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/SamplingAlgorithms/MHProposal.h"
#include "MUQ/SamplingAlgorithms/CrankNicolsonProposal.h"
#include "MUQ/SamplingAlgorithms/ParallelizableMIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SubsamplingMIProposal.h"

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
    comm->Barrier();
    lastState = state;

    //TODO pass directly to exahype
    //std::ofstream file("Input/parameters.csv");
    //file << state->state[0].format(CSVFormat);
    //file.close();

    //Write parameters into ExaHyPE readable format
    std::vector<double> param(state->state[0].size());
    for(int i = 0; i < param.size(); i++){
        param[i] = state->state[0][i];
    }

    //Create some debug output
    if(comm->GetRank()==0){
	    //std::cout << "parameter:" << state->state[0].transpose() << std::endl;

	    std::ofstream ost;
	    ost.open("parameters_r" + std::to_string(globalComm->GetRank()) + ".log", std::ios::app);
	    ost << param[0] << ", " << param[1] << std::endl;
	    ost.close();

	    //std::cout << "Sample number: " << count++ << std::endl;
	    //std::cout << "Parameter " << param[0] << " " << param[1] << std::endl;
    }

    //Discard stupid parameters
    if (param[0] >1.0 || param[0] < -0.1 || param[1]>1.0 || param[1]<-0.1)//reject parameters outside domain
        return -24;

    //run forward model
    int level = index->GetValue(0);
    std::cout << "run_exahype with level " << level << " and global communicator number " << globalComm->GetRank()  << std::endl;
    //muq::setCommunicator(comm->GetMPICommunicator());
    auto output = muq::run_exahype(param,globalComm->GetRank(), level);
    /*for(int i = 0; i< output.size(); i++)
	    std::cout << "output " << i << " is " << output[i] << std::endl;*/

    comm->Barrier();
    double sigma = 1.0;
    return calculateLikelihood(output,globalComm->GetRank());// - 0.5/(sigma*sigma)*state->state[0].squaredNorm();
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
  MyMIComponentFactory(std::shared_ptr<parcer::Communicator> globalComm) : _globalComm(globalComm) {}

  virtual std::shared_ptr<MCMCProposal> Proposal (std::shared_ptr<MultiIndex> const& index, std::shared_ptr<AbstractSamplingProblem> const& samplingProblem) override {
    pt::ptree pt;
    pt.put("BlockIndex",0);

    /*Eigen::VectorXd mu = Eigen::VectorXd::Zero(NUM_PARAM);
    Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(NUM_PARAM,NUM_PARAM);

    auto prior = std::make_shared<Gaussian>(mu, cov);

    return std::make_shared<CrankNicolsonProposal>(pt, samplingProblem, prior);*/

    auto mu = Eigen::VectorXd::Zero(NUM_PARAM);
    Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(NUM_PARAM, NUM_PARAM);
    cov *= 0.005;

    auto prior = std::make_shared<Gaussian>(mu, cov);

    return std::make_shared<MHProposal>(pt, samplingProblem, prior);
  }

  void SetComm(std::shared_ptr<parcer::Communicator> const& comm) override {
    _comm = comm;
  }

  virtual std::shared_ptr<MultiIndex> FinestIndex() override {
    auto index = std::make_shared<MultiIndex>(1);
    index->SetValue(0, 1);
    return index;
  }

  virtual std::shared_ptr<MCMCProposal> CoarseProposal (std::shared_ptr<MultiIndex> const& index,
                                                        std::shared_ptr<AbstractSamplingProblem> const& coarseProblem,
                                                        std::shared_ptr<SingleChainMCMC> const& coarseChain) override {
    pt::ptree ptProposal;
    ptProposal.put("BlockIndex",0);
    int subsampling = 5;
    ptProposal.put("Subsampling", subsampling);
    return std::make_shared<SubsamplingMIProposal>(ptProposal, coarseProblem, coarseChain);
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
    start(0) = .6;
    start(1) = .6;
    return start;
  }

private:
  std::shared_ptr<parcer::Communicator> _comm;
  std::shared_ptr<parcer::Communicator> _globalComm;

};
