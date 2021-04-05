#include "InitialData.h"

#if defined(USE_ASAGI)
#include "easi/YAMLParser.h"
#include "easi/ResultAdapter.h"
#include "reader/asagi_reader.h"
#endif
#include <cmath>
#include "muq_globals.h"

using namespace std;
///// 2D /////

#ifdef Dim2

InitialData::InitialData()
	: scenario(){
#if defined(USE_ASAGI)
		std::cout << "Initialising with ASAGI" << std::endl;
		asagiReader = new AsagiReader("");
		parser = new easi::YAMLParser(3, asagiReader);
		model  = parser->parse("data.yaml");
#endif
	}

InitialData::InitialData(int a_scenario, char* filename)
	: scenario(a_scenario){
#if defined(USE_ASAGI)
		std::cout << "Initialising with ASAGI" << std::endl;
		asagiReader = new AsagiReader("");
		parser = new easi::YAMLParser(3, asagiReader);
		model  = parser->parse(filename);
#endif
	}

InitialData::~InitialData(){
#if defined(USE_ASAGI)
	delete asagiReader;
	delete parser;
	delete model;
#endif
}


#if defined(USE_ASAGI)
void InitialData::readAsagiData_nobath(const double* const x,double* Q){
	double displacement[1], bathymetry[1];
	displacement[0] = 0.0; bathymetry[0] = -4.95;

	easi::ArraysAdapter<double> adapter;
	adapter.addBindingPoint("d",displacement);

	easi::Query query(1,3);
	query.x(0,0)=x[0] - muq::param[0];
	query.x(0,1)=x[1] - muq::param[1];
	query.x(0,2)=0;
	model->evaluate(query,adapter);

	Q[0]=max(0.0,-bathymetry[0]); //h = H-b
	if(std::isnan(Q[0]))
		std::cout << "Error when reading bathymetry" << x[0] << "," << x[1] << std::endl;
	Q[1]= 0;
	Q[2]= 0;
	Q[3]= displacement[0]+bathymetry[0]; 
}

void InitialData::readAsagiData_onlybath(const double* const x,double* Q){
	double bathymetry[1], displacement[1];
	bathymetry[0] = 0.0; displacement[0] = 0.0;

	easi::ArraysAdapter<double> adapter;
	adapter.addBindingPoint("b",bathymetry);

	easi::Query query(1,3);
	query.x(0,0)=x[0];
	query.x(0,1)=x[1];
	query.x(0,2)=0;
	model->evaluate(query,adapter);

	easi::ArraysAdapter<double> adapter2;
	adapter2.addBindingPoint("d",displacement);

	easi::Query query2(1,3);
	query2.x(0,0)=x[0] - muq::param[0];
	query2.x(0,1)=x[1] - muq::param[1];
	query2.x(0,2)=0;
	model->evaluate(query2,adapter2);

	Q[3]= displacement[0]+bathymetry[0]; 
}

void InitialData::readAsagiData(const double* const x,double* Q){
	double bathymetry[1], displacement[1];
	bathymetry[0] = 0.0; displacement[0] = 0.0;

	easi::ArraysAdapter<double> adapter;
	adapter.addBindingPoint("b",bathymetry);

	easi::Query query(1,3);
	query.x(0,0)=x[0];
	query.x(0,1)=x[1];
	query.x(0,2)=0;
	model->evaluate(query,adapter);

	easi::ArraysAdapter<double> adapter2;
	adapter2.addBindingPoint("d",displacement);

	easi::Query query2(1,3);
	query2.x(0,0)=x[0] - muq::param[0];
	query2.x(0,1)=x[1] - muq::param[1];
	query2.x(0,2)=0;
	model->evaluate(query2,adapter2);

	Q[0]=std::max(0.0,-bathymetry[0]);
	if(std::isnan(Q[0]))
		std::cout << "Error when reading bathymetry" << x[0] << "," << x[1] << std::endl;
	Q[1]= 0;
	Q[2]= 0;
	Q[3]= displacement[0]+bathymetry[0]; 
}
#endif






void InitialData::getInitialData(const double* const x,double* Q) {
	switch (scenario)
	{
		case 14:
			readAsagiData(x, Q);
			break;
		case 15:
			readAsagiData_nobath(x, Q);
			break;
		case 16:
			readAsagiData_onlybath(x, Q);
			break;
	}
}

#endif
