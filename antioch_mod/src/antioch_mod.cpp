//////-----------------------------------------------------------------------bl-
//////--------------------------------------------------------------------------
//////
////// Antioch - A Gas Dynamics Thermochemistry Library
//////
////// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//////                         Sylvain Plessis, Roy H. Stonger
//////
////// Copyright (C) 2013 The PECOS Development Team
//////
////// This library is free software; you can redistribute it and/or
////// modify it under the terms of the Version 2.1 GNU Lesser General
////// Public License as published by the Free Software Foundation.
//////
////// This library is distributed in the hope that it will be useful,
////// but WITHOUT ANY WARRANTY; without even the implied warranty of
////// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
////// Lesser General Public License for more details.
//////
////// You should have received a copy of the GNU Lesser General Public
////// License along with this library; if not, write to the Free Software
////// Foundation, Inc. 51 Franklin Street, Fifth Floor,
////// Boston, MA  02110-1301  USA
//////
//////-----------------------------------------------------------------------el-
//////
////// $Id$
//////
//////--------------------------------------------------------------------------
//////--------------------------------------------------------------------------
////
////// C++
////#include <iostream>
////#include <cmath>
////
////// Antioch
////#include "antioch/blottner_viscosity.h"
////using namespace Antioch;
////
////template <typename Scalar>
////int test_viscosity( const Scalar mu, const Scalar mu_exact, const Scalar tol )
////{
////  using std::abs;
////
////  int return_flag = 0;
////
////  const double rel_error = abs( (mu - mu_exact)/mu_exact);
////
////  if( rel_error  > tol )
////    {
////      std::cerr << "Error: Mismatch in viscosity" << std::endl
////		<< "mu(T)    = " << mu << std::endl
////		<< "mu_exact = " << mu_exact << std::endl
////		<< "rel_error = " << rel_error << std::endl
////		<< "tol = " << tol << std::endl;
////      return_flag = 1;
////    }
////
////  return return_flag;
////}
////
////template <typename Scalar>
////int tester()
////{
////  const Scalar a = 3.14e-3;
////  const Scalar b = 2.71e-2;
////  const Scalar c = 42.0e-5;
////
////  Antioch::BlottnerViscosity<Scalar> mu(a,b,c);
////
////  std::cout << mu << std::endl;
////
////  const Scalar T = 1521.2;
////
////  // octave gives
////  const Scalar mu_exact = 0.144422234167703;
////
////  int return_flag = 0;
////
////  const Scalar tol = 1.0e-14;
////
////  return_flag = test_viscosity( mu(T), mu_exact, tol );
////
////  const Scalar a2 = 1e-3;
////  const Scalar b2 = 2e-2;
////  const Scalar c2 = 3e-5;
////
////  mu.reset_coeffs( a2, b2, c2 );
////
////  // octave gives
////  const Scalar mu_exact2 = 0.122172495548880;
////
////  return_flag = test_viscosity( mu(T), mu_exact2, tol );
////
////  return return_flag;
////}
////
////int main()
////{
////	std::cout << "Hello... there are errors, but I can still compile" << std::endl << std::endl;
////  return (tester<double>() ||
////          tester<long double>() ||
////          tester<float>());
////}
//
////-----------------------------------------------------------------------bl-
////--------------------------------------------------------------------------
////
//// Antioch - A Gas Dynamics Thermochemistry Library
////
//// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
////                         Sylvain Plessis, Roy H. Stonger
////
//// Copyright (C) 2013 The PECOS Development Team
////
//// This library is free software; you can redistribute it and/or
//// modify it under the terms of the Version 2.1 GNU Lesser General
//// Public License as published by the Free Software Foundation.
////
//// This library is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//// Lesser General Public License for more details.
////
//// You should have received a copy of the GNU Lesser General Public
//// License along with this library; if not, write to the Free Software
//// Foundation, Inc. 51 Franklin Street, Fifth Floor,
//// Boston, MA  02110-1301  USA
////
////-----------------------------------------------------------------------el-
////
//// $Id$
////
////--------------------------------------------------------------------------
////--------------------------------------------------------------------------
//
//// C++
//#include <limits>
//#include <string>
//#include <vector>
//
//// Antioch
//#include "antioch/vector_utils.h"
//
//#include "antioch/antioch_asserts.h"
//#include "antioch/chemical_species.h"
//#include "antioch/chemical_mixture.h"
//#include "antioch/reaction_set.h"
//#include "antioch/read_reaction_set_data.h"
//#include "antioch/cea_evaluator.h"
//#include "antioch/cea_mixture_ascii_parsing.h"
//#include "antioch/kinetics_evaluator.h"
//
//
//template <typename Scalar>
//int checker(const Scalar & theory, const Scalar & computed, const std::string& words)
//{
//
//	int return_flag(0);
//	const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 500;
//
//	const Scalar rel_error = std::abs( (computed - theory)/theory);
//	if( rel_error > tol )
//	{
//		std::cerr << "Error: Mismatch between theory and regression values in test " << words << std::endl;
//		std::cout << std::scientific << std::setprecision(16)
//		<< "theory value        = " << theory    << std::endl
//		<< "computed value      = " << computed  << std::endl
//		<< "relative difference = " << rel_error << std::endl
//		<< "tolerance           = " << tol       << std::endl << std::endl;
//		return_flag = 1;
//	}
//
//	return return_flag;
//}
//
//
//template <typename Scalar>
//int tester(const std::string& input_name)
//{
//	using std::abs;
//
//	std::vector<std::string> species_str_list;
//	const unsigned int n_species = 5;
//	species_str_list.reserve(n_species);
//	species_str_list.push_back( "N2" );
//	species_str_list.push_back( "O2" );
//	species_str_list.push_back( "N" );
//	species_str_list.push_back( "O" );
//	species_str_list.push_back( "NO" );
//
//	Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
//	Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );
//
//	Antioch::CEAThermoMixture<Scalar> cea_mixture( chem_mixture );
//	Antioch::read_cea_mixture_data_ascii( cea_mixture, Antioch::DefaultFilename::thermo_data() );
//	Antioch::CEAEvaluator<Scalar> thermo( cea_mixture );
//
//	Antioch::read_reaction_set_data_xml<Scalar>( input_name, true, reaction_set );
//
//	const Scalar T = 1500.0; // K
//	const Scalar P = 1.0e5; // Pa
//	const Antioch::KineticsConditions<Scalar> conditions(T);
//
//	// Mass fractions
//	std::vector<Scalar> Y(n_species,0.2);
//
//	const Scalar R_mix = chem_mixture.R(Y); // get R_tot in J.kg-1.K-1
//
//	const Scalar rho = P/(R_mix*T); // kg.m-3
//
//	std::vector<Scalar> molar_densities(n_species,0.0);
//	chem_mixture.molar_densities(rho,Y,molar_densities);
//
//	std::vector<Scalar> h_RT_minus_s_R(n_species);
//	std::vector<Scalar> dh_RT_minus_s_R_dT(n_species);
//	Antioch::TempCache<Scalar> temp_cache(T);
//
//	thermo.h_RT_minus_s_R(temp_cache,h_RT_minus_s_R);
//	thermo.dh_RT_minus_s_R_dT(temp_cache,dh_RT_minus_s_R_dT);
//
//	Antioch::KineticsEvaluator<Scalar> kinetics( reaction_set, 0 );
//
//	std::vector<Scalar> omega_dot(n_species);
//	std::vector<Scalar> omega_dot_2(n_species);
//	std::vector<Scalar> omega_dot_3(n_species);
//	std::vector<Scalar> omega_dot_4(n_species);
//	std::vector<Scalar> domega_dot_dT(n_species);
//	std::vector<Scalar> domega_dot_dT_2(n_species);
//
//	std::vector<std::vector<Scalar> > domega_dot_drho_s(n_species);
//	std::vector<std::vector<Scalar> > domega_dot_drho_s_2(n_species);
//	for( unsigned int s = 0; s < n_species; s++ )
//	{
//		domega_dot_drho_s[s].resize(n_species);
//		domega_dot_drho_s_2[s].resize(n_species);
//	}
//
//	// backward compatibility
//
//	kinetics.compute_mass_sources( T , molar_densities, h_RT_minus_s_R, omega_dot);
//
//	kinetics.compute_mass_sources_and_derivs( T , molar_densities, h_RT_minus_s_R, dh_RT_minus_s_R_dT,
//			omega_dot_2, domega_dot_dT, domega_dot_drho_s );
//
//	// kinetics conditions
//
//	kinetics.compute_mass_sources( conditions , molar_densities, h_RT_minus_s_R, omega_dot_3);
//
//	kinetics.compute_mass_sources_and_derivs( conditions , molar_densities, h_RT_minus_s_R, dh_RT_minus_s_R_dT,
//			omega_dot_4, domega_dot_dT_2, domega_dot_drho_s_2 );
//
//	for( unsigned int s = 0; s < n_species; s++)
//	{
//		std::cout << std::scientific << std::setprecision(16)
//		<< "domega_dot_dT(" << chem_mixture.chemical_species()[s]->species() << ") = "
//		<< domega_dot_dT[s] << std::endl;
//	}
//
//	for( unsigned int s = 0; s < n_species; s++)
//	{
//		for( unsigned int t = 0; t < n_species; t++)
//		{
//			std::cout << std::scientific << std::setprecision(16)
//			<< "domega_dot_drho_s(" << chem_mixture.chemical_species()[s]->species()
//			<< ", " << chem_mixture.chemical_species()[t]->species() << ") = "
//			<< domega_dot_drho_s[s][t] << std::endl;
//		}
//	}
//
//	int return_flag = 0;
//
//	// Regression values for omega_dot
//	std::vector<Scalar> omega_dot_reg(n_species);
//	omega_dot_reg[0] =  9.1623705357123753e+04;
//	omega_dot_reg[1] = -3.3462025680272243e+05;
//	omega_dot_reg[2] = -2.1139216712069495e+05;
//	omega_dot_reg[3] =  1.9782018625609628e+05;
//	omega_dot_reg[4] =  2.5656853231019735e+05;
//
//
//	// omega_dots tests
//	for( unsigned int s = 0; s < n_species; s++)
//	{
//		return_flag = checker(omega_dot_reg[s],omega_dot_2[s],"omega dot2 of species " + chem_mixture.chemical_species()[s]->species()) || return_flag;
//		return_flag = checker(omega_dot_reg[s],omega_dot_3[s],"omega dot3 of species " + chem_mixture.chemical_species()[s]->species()) || return_flag;
//		return_flag = checker(omega_dot_reg[s],omega_dot_4[s],"omega dot4 of species " + chem_mixture.chemical_species()[s]->species()) || return_flag;
//	}
//
//	// Regression values for domega_dot_dT
//	std::vector<Scalar> domega_dot_reg_dT(n_species);
//	domega_dot_reg_dT[0] =  1.8014990183270937e+02;
//	domega_dot_reg_dT[1] = -5.2724437115534380e+02;
//	domega_dot_reg_dT[2] = -3.0930094476883017e+02;
//	domega_dot_reg_dT[3] =  3.7972747459781005e+02;
//	domega_dot_reg_dT[4] =  2.7666793949365456e+02;
//
//	for( unsigned int s = 0; s < n_species; s++)
//	{
//		return_flag = checker(domega_dot_reg_dT[s],domega_dot_dT[s],"domega_dot_dT of species "   + chem_mixture.chemical_species()[s]->species()) || return_flag;
//		return_flag = checker(domega_dot_reg_dT[s],domega_dot_dT_2[s],"domega_dot_dT2 of species "  + chem_mixture.chemical_species()[s]->species()) || return_flag;
//	}
//
//	// Regression values for domega_dot_drho_s
//	std::vector<std::vector<Scalar> > domega_dot_reg_drhos(n_species);
//	for( unsigned int s = 0; s < n_species; s++)
//	{
//		domega_dot_reg_drhos[s].resize(n_species);
//	}
//
//	domega_dot_reg_drhos[0][0] = 1.9675775188085109e+04;
//	domega_dot_reg_drhos[0][1] = 1.7226141262419737e+04;
//	domega_dot_reg_drhos[0][2] = 3.2159299284723610e+06;
//	domega_dot_reg_drhos[0][3] = 1.4765214711933021e+05;
//	domega_dot_reg_drhos[0][4] = 2.3225053279918131e+06;
//
//	domega_dot_reg_drhos[1][0] =  8.8927385505978492e+03;
//	domega_dot_reg_drhos[1][1] = -9.9560178070099482e+06;
//	domega_dot_reg_drhos[1][2] = -9.8748760140991123e+06;
//	domega_dot_reg_drhos[1][3] =  4.6143036700500813e+05;
//	domega_dot_reg_drhos[1][4] =  8.3487375168772399e+03;
//
//	domega_dot_reg_drhos[2][0] = -2.2420842426881281e+04;
//	domega_dot_reg_drhos[2][1] = -4.3812843857644886e+06;
//	domega_dot_reg_drhos[2][2] = -6.8343593463263955e+06;
//	domega_dot_reg_drhos[2][3] = -5.4143671040862988e+05;
//	domega_dot_reg_drhos[2][4] = -1.2267997668149246e+06;
//
//	domega_dot_reg_drhos[3][0] = -1.2028166578920147e+04;
//	domega_dot_reg_drhos[3][1] =  4.9713710400172938e+06;
//	domega_dot_reg_drhos[3][2] =  5.7418898143800552e+06;
//	domega_dot_reg_drhos[3][3] = -9.1121284934572734e+05;
//	domega_dot_reg_drhos[3][4] =  1.2431710353864791e+06;
//
//	domega_dot_reg_drhos[4][0] =  5.8804952671184686e+03;
//	domega_dot_reg_drhos[4][1] =  9.3487050114947233e+06;
//	domega_dot_reg_drhos[4][2] =  7.7514156175730915e+06;
//	domega_dot_reg_drhos[4][3] =  8.4356704563001888e+05;
//	domega_dot_reg_drhos[4][4] = -2.3472253340802449e+06;
//
//	for( unsigned int s = 0; s < n_species; s++)
//	{
//		for( unsigned int t = 0; t < n_species; t++)
//		{
//			return_flag = checker(domega_dot_reg_drhos[s][t],domega_dot_drho_s[s][t]  , "domega_dot_drhos of species "
//					+ chem_mixture.chemical_species()[s]->species()
//					+ " with respect to species "
//					+ chem_mixture.chemical_species()[t]->species()) || return_flag;
//			return_flag = checker(domega_dot_reg_drhos[s][t],domega_dot_drho_s_2[s][t], "domega_dot_drhos of species "
//					+ chem_mixture.chemical_species()[s]->species()
//					+ " with respect to species "
//					+ chem_mixture.chemical_species()[t]->species()) || return_flag;
//		}
//	}
//
//	// now some resetting and verifying omega_dot
//	std::string reaction_id("0001");
//	std::vector<std::string> keywords;
//	keywords.push_back("A");
//	Scalar new_value(6e15L); //SI, original is 7e15
//	reaction_set.set_parameter_of_reaction(reaction_id,keywords,new_value);
//	reaction_id = "0002";
//	keywords[0] = "efficiencies";
//	keywords.push_back("O2");
//	new_value = 1.2L;
//	reaction_set.set_parameter_of_reaction(reaction_id,keywords,new_value);
//
//	//recomputing
//	Antioch::set_zero(omega_dot);
//	kinetics.compute_mass_sources( conditions , molar_densities, h_RT_minus_s_R, omega_dot);
//
//	// new values, SI
//	omega_dot_reg[0] =  8.9806036413183845e4;
//	omega_dot_reg[1] = -3.3456693672788515e5;
//	omega_dot_reg[2] = -2.0957449817675504e5;
//	omega_dot_reg[3] =  1.9776686618125900e5;
//	omega_dot_reg[4] =  2.5656853231019735e5;
//
//	for( unsigned int s = 0; s < n_species; s++)
//	{
//		return_flag = checker(omega_dot_reg[s],omega_dot[s] ,"resetted omega dot of species "  + chem_mixture.chemical_species()[s]->species()) || return_flag;
//	}
//
//	// and now a get/set loop
//	reaction_id = "0004";
//	const Scalar val(1.1L);
//	keywords.clear();
//	keywords.push_back("E");
//	//  keywords.push_back("cal/mol"); // no, you only get SI
//	reaction_set.set_parameter_of_reaction(reaction_id,keywords,
//			val * reaction_set.get_parameter_of_reaction(reaction_id,keywords) );
//
//	//recomputing
//	Antioch::set_zero(omega_dot);
//	kinetics.compute_mass_sources( conditions , molar_densities, h_RT_minus_s_R, omega_dot);
//
//	// new values, SI
//	omega_dot_reg[0] =  1.541307374714467399842142e4;
//	omega_dot_reg[1] = -3.345669367278851525665733e5;
//	omega_dot_reg[2] = -1.723780168437354542966397e5;
//	omega_dot_reg[3] =  1.552808795073360031682657e5;
//	omega_dot_reg[4] =  3.362510003171399296965259e5;
//
//	for( unsigned int s = 0; s < n_species; s++)
//	{
//		return_flag = checker(omega_dot_reg[s],omega_dot[s] ,"loop-resetted omega dot of species "  + chem_mixture.chemical_species()[s]->species()) || return_flag;
//	}
//
//	return return_flag;
//}
//
//int main(int argc, char* argv[])
//{
//	// xml file to be read in
//	std::string readXML{"share/xml_inputs/GRIMech30.xml"};
////	std::string readXML{"test/input_files/air_5sp.xml"};
//	// Check command line count.
//	if( argc < 2 )
//	{
//		// TODO: Need more consistent error handling.
//		std::cerr << "Error: Must specify reaction set XML input file." << std::endl;
//		antioch_error();
//	}
//	std::cout << "Read in my own xml file" << std::endl;
//	return (tester<float>(readXML) ||
//			tester<double>(readXML) /* ||
//	          tester<long double>(std::string(argv[1]))*/
//	);
//
//	//	return (tester<float>(std::string(argv[1])) ||
//	//			tester<double>(std::string(argv[1])) /* ||
//	//          tester<long double>(std::string(argv[1]))*/
//	//	);
//}
//

/** Setup an ignition test that is similar to the ignition test in SIMIT, this will require the use of an integrator **/
//include standard libraries
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <vector>
#include <map>
#include <time.h>

// antioch libraries
#include "antioch/vector_utils.h"
#include "antioch/antioch_asserts.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/read_reaction_set_data.h"
#include "antioch/cea_evaluator.h"
#include "antioch/cea_mixture_ascii_parsing.h"
#include "antioch/kinetics_evaluator.h"
#include "antioch/nasa_mixture_parsing_instantiate_macro.h"
#include "antioch/parsing_enum.h"
#include "antioch/ascii_parser.h"
#include "antioch/xml_parser.h"
#include "antioch/chemkin_parser.h"
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_mixture_parsing.h"

// cvode libraries
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <nvector/nvector_serial.h>
#include <cvodes/cvodes.h>
//#include <cvodes/cvodes_sptfqmr.h>
//#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_dense.h>
//#include <cvodes/cvodes_lapack.h>
//#include <cvodes/cvodes_band.h>

static const double Patm{101325.0}; //constant atmospheric pressure in Pa (to convert from atm to Pa)
static const double Tmin{200.0};    // minumum temperatuer allowed
static const double Tmax{10000.0};  // maximum temperature allowed
static const double eps{1.e-6};     // standard error tolerance
static const double tiny{1.e-12};   // tiny number to keep from dividing by zero
//static const double MW_air{0.02885};       // the molecular weight of air (kg/mol)
static const double MW_H{1.007940e-03};         // molecular weight of H (kg/mol)
static const double MW_C{1.201100e-02};         // molecular weight of C (kg/mol)
static const double MW_N{1.400674e-02};         // molecular weight of O (kg/mol)
static const double MW_O{1.599940e-02};         // molecular weight of N (kg/mol)

void setSpecList(std::vector<std::string>& sL)
{
	const unsigned int nspec{53};
	sL.reserve(nspec);
	sL.push_back("CH4");
	sL.push_back("H");
	sL.push_back("O");
	sL.push_back("O2");
	sL.push_back("OH");
	sL.push_back("H2O");
	sL.push_back("HO2");
	sL.push_back("H2O2");
	sL.push_back("C");
	sL.push_back("CH");
	sL.push_back("CH2");
	sL.push_back("CH2(S)");
	sL.push_back("CH3");
	sL.push_back("H2");
	sL.push_back("CO");
	sL.push_back("CO2");
	sL.push_back("HCO");
	sL.push_back("CH2O");
	sL.push_back("CH2OH");
	sL.push_back("CH3O");
	sL.push_back("CH3OH");
	sL.push_back("C2H");
	sL.push_back("C2H2");
	sL.push_back("C2H3");
	sL.push_back("C2H4");
	sL.push_back("C2H5");
	sL.push_back("C2H6");
	sL.push_back("HCCO");
	sL.push_back("CH2CO");
	sL.push_back("HCCOH");
	sL.push_back("N");
	sL.push_back("NH");
	sL.push_back("NH2");
	sL.push_back("NH3");
	sL.push_back("NNH");
	sL.push_back("NO");
	sL.push_back("NO2");
	sL.push_back("N2O");
	sL.push_back("HNO");
	sL.push_back("CN");
	sL.push_back("HCN");
	sL.push_back("H2CN");
	sL.push_back("HCNN");
	sL.push_back("HCNO");
	sL.push_back("HOCN");
	sL.push_back("HNCO");
	sL.push_back("NCO");
	sL.push_back("C3H7");
	sL.push_back("C3H8");
	sL.push_back("CH2CHO");
	sL.push_back("CH3CHO");
	sL.push_back("AR");
	sL.push_back("N2");
	return;
}

typedef struct{
	Antioch::KineticsEvaluator<double> *kinetics;
	Antioch::ChemicalMixture<double> *chem_mixture;
	Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double>> *nasa_thermo;
	Antioch::ReactionSet<double> *reaction_set;
	int nspec;
	int i_fuel=0, i_OH=0, i_O2=0, i_N2=0;  // setup the index locations of the fuel, OH, O2, and N2
	double T =0., P=0., rho=0., Ener=0., Enthalpy=0., Etot=0.; // (units) K, Pa, kg/m^3, J/kg, J/kg
	std::vector<double> Y;               // species mass fraction
	std::vector<std::string> specList;   // species list
	//	std::map<std::string, double> map_specToYi; // map the species name to Yi location
	std::vector<double> MW; //vector to store the molecular weights (kg/mol)
	std::vector<double> Hf; //vector to store heats of formation in (J/kg)
	std::vector<double> w_dot;           // source term for species
	std::vector<double> OH_trace;        // track OH in time
	std::vector<double> P_trace;         // track P in time
	std::vector<double> T_trace;         // track T in time
	std::vector<double> time;            // track time
	std::vector<double> _omega_dot;      // temporary array for the f function
	std::vector<double> _H_RT_minus_S;   // temp array for f function
	std::vector<double> _mole_dens;      // temp arrat for f funtion
} *DataForVODE;

template <typename Scalar=double>
class IgnitionTimes
{
public:
	IgnitionTimes(std::string,std::string,double,double,double);
	~IgnitionTimes();
	int initialize(std::string, std::vector<std::string>);
	static double totalEner(double, std::vector<double>, int, std::vector<double>, double, std::vector<double>, Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double>>*);
	static double T_from_Ener(double, double, std::vector<double>, int, std::vector<double>, std::vector<double>, Antioch::ChemicalMixture<double>*, Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double>>*);
	static int f(realtype t, N_Vector y, N_Vector y_dot, void *f_data);
	static void zeroArray(std::vector<double>&);
	int check_flag(void* flagvalue,const char *funcname, int opt);
	double omega_dot(double dt);
	void normalizeYi(std::vector<Scalar>&);
	void RK_Integrator(Scalar);
	void IgnitionSolver();
	std::vector<Scalar> getYi(void);
	Scalar getRho(void);
	Scalar getP(void);
	Scalar getT(void);
	void setYi(std::vector<Scalar>);
	void setRho(Scalar);
	void setP(Scalar);
	void setT(Scalar);

	//freely access sundial data formats
	DataForVODE data;
	N_Vector y {NULL};             // needed for the values of unknowns in time derivatitve
	N_Vector y_dot {NULL};         // time derivative of the y vector
	void* cvode_mem {NULL}; // cvode_mem needed in ODE solver
//	realtype reltol {1.e-11};// relative tolerance (what is used in Cantera)
//	realtype abstol {1.e-13};// absolute tolerance (what is used in Cantera)
	realtype reltol {1.e-11};// relative tolerance
	realtype abstol {1.e-14};// absolute tolerance
	int cflag{0};           // int for error flagging needed for cvode
	double dt;                           // time step (s)
	double ER;                            // temperature (K)
	std::string fuel;       // name of the fuel
	std::vector<std::string> specList;

private:
	//it seems like when the original calls are used to setup pointed in "data" then they are deleted
	// once we leave the scope of the initialize function
	Antioch::ChemicalMixture<Scalar>* chem_mixture1;
	Antioch::NASAThermoMixture<Scalar, Antioch::NASA7CurveFit<Scalar> >* nasa_mixture1;
	Antioch::NASAEvaluator<Scalar, Antioch::NASA7CurveFit<Scalar>>* nasa_thermo1;
	Antioch::ReactionSet<Scalar>* reaction_set1;
	Antioch::KineticsEvaluator<Scalar>* kinetics1;
};

template<typename Scalar>
IgnitionTimes<Scalar>::IgnitionTimes(std::string filename, std::string fuel, double ER, double P_in, double T_in):
dt{1.e-9}, ER{ER}, fuel{fuel}
{
	data = (DataForVODE) malloc(sizeof *data);
	data->P = P_in*Patm; //given in atm, convert to Pa
	data->T = T_in;      // set initial temperature in K


	//	std::vector<std::string> specList;
	setSpecList(specList);
	data->specList = specList;
	data->nspec = specList.size();
	//reserve space for vectors in data
	data->Hf.resize(data->nspec);
	data->MW.resize(data->nspec);
	data->Y.resize(data->nspec);
	data->w_dot.resize(data->nspec);
	data->_H_RT_minus_S.resize(data->nspec);
	data->_omega_dot.resize(data->nspec);
	data->_mole_dens.resize(data->nspec);
	chem_mixture1 = new Antioch::ChemicalMixture<Scalar>(specList);
	reaction_set1 = new Antioch::ReactionSet<Scalar>(*chem_mixture1);
	nasa_mixture1 = new Antioch::NASAThermoMixture<Scalar, Antioch::NASA7CurveFit<Scalar>>(*chem_mixture1);
	Antioch::read_nasa_mixture_data<Scalar, Antioch::NASA7CurveFit<Scalar>> (*nasa_mixture1, filename, Antioch::XML, true);
	nasa_thermo1 = new Antioch::NASAEvaluator<Scalar, Antioch::NASA7CurveFit<Scalar>>(*nasa_mixture1);
//	nasa_thermo1->initHof();
	bool verbose_read = true;
	Antioch::read_reaction_set_data_xml(filename, verbose_read, *reaction_set1);
	kinetics1 = new Antioch::KineticsEvaluator<Scalar>(*reaction_set1, 0);
	data->chem_mixture = chem_mixture1;
	data->nasa_thermo = nasa_thermo1;
	data->reaction_set = reaction_set1;
	data->kinetics = kinetics1;



	int i = this->initialize(filename, specList);
	if(i!=0)
		std::cout << "Error found in the initializer() function call" << std::endl;
	return;
}

template<typename Scalar>
int IgnitionTimes<Scalar>::initialize(std::string filename, std::vector<std::string> specList){
	const Scalar T_test = 1500.; //K
	const Scalar P_test = 101325.0; // Pa
	const Antioch::KineticsConditions<Scalar> conditions(T_test);

	// Initialize the mass Fractions
	normalizeYi(data->Y); //make sure species mass fractions sum to 1

	//TODO get the mapping to work
	// setup the mapping between species name and mass fraction
	//	data->map_specToYi[data->specList[0]] = data->Y[0];
	//	data->setupMap();
	//	std::map<std::string, double> map1;
	//	for(unsigned i = 0; i<specList.size(); ++i){
	//		map1[specList[i]] = data->Y[i];
	//	}
	//	data->map_specToYi = map1;

	//find the index of the fuel, OH, O2, and N2 in the specList and correspond this to the index location in Y
	auto it = std::find(specList.begin(), specList.end(), fuel);
	if(it == specList.end()){
		std::cout << "No such species in specList when looking for: " << fuel << std::endl;
	}else {data->i_fuel = std::distance(specList.begin(), it);}
	it = std::find(specList.begin(), specList.end(), "OH");
	if(it == specList.end()){
		std::cout << "No such species in specList when looking for: " << "OH" << std::endl;
	}else {data->i_OH = std::distance(specList.begin(), it);}
	it = std::find(specList.begin(), specList.end(), "O2");
	if(it == specList.end()){
		std::cout << "No such species in specList when looking for: " << "O2" << std::endl;
	}else {data->i_O2 = std::distance(specList.begin(), it);}
	it = std::find(specList.begin(), specList.end(), "N2");
	if(it == specList.end()){
		std::cout << "No such species in specList when looking for: " << "N2" << std::endl;
	}else {data->i_N2 = std::distance(specList.begin(), it);}

	// setup the mass fractions based on the equivalence ratio provided
	// initially hard coded for methane (CH4)
	// CxHyOz + (x+y/4-z/2)(O2+3.76N2) --> y/2 H2O + x CO2 + (x+y/4-z/2)3.76 N2
	double x_c {1.0};
	double y_h {4.0};
	double z_o {0.0};
	double MW_fuel = x_c*MW_C + y_h*MW_H + z_o*MW_O;
	double a = x_c + 0.25*y_h + 0.5*z_o; // a = (x + y/4 + z/2)
	//	double AF_stoic = (4.76*a)*(MW_air/MW_fuel); // air-to-fuel ratio at stoichiometric proportions
	//	double AF = AF_stoic/ER; // actual air-to-fuel ratio based on the equivalence ratio provided

	// determine mole fractions of premixed fuel and air at a given ER
	double X_fuel = 1./(1.+((4.762*a)/ER)); // mole fraction of fuel
	double X_air  = 1-X_fuel;              // mole fraction of the air
	double X_O2   = X_air/4.762;           // mole fraction of O2
	double X_N2   = 3.762*X_O2;            // mole fraction of N2
	double MWmix  = X_fuel*MW_fuel + X_O2*2.*MW_O + X_N2*2.*MW_N; // molecular weight of premixed fuel

	// convert mole fractions to mass fraction
	double Yf  = (MW_fuel*X_fuel)/MWmix;  // mass fraction of fuel
	double YO2 = (2.*MW_O*X_O2)/MWmix;    // mass fraction of O2
	double YN2 = (2.*MW_N*X_N2)/MWmix;    // mass fraction of N2

	// set the mass fractions of the fuel, oxygen, and nitrogen
	//	data->map_specToYi[fuel] = Yf;
	//	data->map_specToYi["O2"] = YO2;
	//	data->map_specToYi["N2"] = YN2;
	data->Y[data->i_fuel] = Yf;
	data->Y[data->i_O2]   = YO2;
	data->Y[data->i_N2]   = YN2;
	normalizeYi(data->Y);

	// setup and store vectors for the molcular weight and heat of formation
	std::vector<double> MW1(data->nspec, 0.);
	for(int i=0; i<data->nspec; i++){
		MW1[i] = data->chem_mixture->M(i);
	}
	//	data->chem_mixture[0].M(MW1); // MW should be in kg/mol
	data->MW = MW1;
	std::vector<Scalar> h(data->nspec);
	Scalar TRef = 298.0;
	Antioch::TempCache<Scalar> temp_cache1(TRef);
	data->nasa_thermo->h(temp_cache1,h); // h should be in J/mol
	for(int i=0; i<data->nspec; i++)
		data->Hf[i] = h[i]; //Hf is in J/kg

	const Scalar R_mix = data->chem_mixture->R(data->Y);
	data->rho = data->P/(R_mix*data->T); // set the density in the data

	// set the total enthalpy and enerygy
	data->Ener = totalEner(data->T, data->Y,data->nspec, data->MW, R_mix, data->Hf, data->nasa_thermo);
	data->Enthalpy = data->Ener + data->P/data->rho;

	const Scalar rho = P_test/(R_mix*T_test); // kg/m^3

	std::vector<Scalar> mol_dens(data->nspec,0.0);
	data->chem_mixture->molar_densities(rho,data->Y,mol_dens);

	std::vector<Scalar> h_RT_minus_s_R(data->nspec);
	std::vector<Scalar> dh_RT_minus_s_R_dT(data->nspec);
	Antioch::TempCache<Scalar> temp_cache(T_test);

	data->nasa_thermo->h_RT_minus_s_R(temp_cache,h_RT_minus_s_R);
	data->nasa_thermo->dh_RT_minus_s_R_dT(temp_cache,dh_RT_minus_s_R_dT);

	std::vector<Scalar> omega_dot(data->nspec);
	std::vector<Scalar> omega_dot2(data->nspec);
	std::vector<Scalar> domega_dot_dT(data->nspec);
	std::vector<std::vector<Scalar>> domega_dot_drho_s(data->nspec);
	for( unsigned int s = 0; s < data->nspec; s++ )
	{
		domega_dot_drho_s[s].resize(data->nspec);
	}

//	data->kinetics->compute_mass_sources(T_test, mol_dens, h_RT_minus_s_R, omega_dot);
//	data->kinetics->compute_mass_sources_and_derivs(T_test, mol_dens, h_RT_minus_s_R, dh_RT_minus_s_R_dT, omega_dot2, domega_dot_dT, domega_dot_drho_s);
//
//
//	for( unsigned int s = 0; s < data->nspec; s++)
//	{
//		std::cout << std::scientific << std::setprecision(16)
//		<< "domega_dot_dT(" << data->chem_mixture->chemical_species()[s]->species() << ") = "
//		<< domega_dot_dT[s] << std::endl;
//	}

	//initialize the CVODES parameters
	y =     N_VNew_Serial(data->nspec);
	y_dot = N_VNew_Serial(data->nspec);
	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON); //recommended setup for stiff problems
	if(check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

	std::vector<Scalar> Ytemp = this->getYi();
	normalizeYi(Ytemp);
	for(int i=0; i<data->nspec; i++){
		NV_Ith_S(y,i) = Ytemp[i];
		NV_Ith_S(y_dot,i) = 0.;
	}

	//intialize the CVODE solver
	double t0 = 0.;
	cflag = CVodeInit(cvode_mem,f,t0,y);
	if(check_flag(&cflag, "CvodeInit", 1)) return(1);

	//set the tolerances
	cflag = CVodeSStolerances(cvode_mem, reltol, abstol);
	if (check_flag(&cflag, "CVodeSStolerances", 1)) return(1);

	// set the pointer to the user-defined data
	cflag = CVodeSetUserData(cvode_mem, data);
	if(check_flag(&cflag, "CVodeSetUserData", 1)) return(1);

	//setup the maximum number of steps the linear solver will take (500 is the default) maxStesps < 0 disables it (not recommended)
	int long maxSteps = 100000; //value set in Cantera
	cflag = CVodeSetMaxNumSteps(cvode_mem, maxSteps);
	if(check_flag(&cflag, "CVodeSetMaxSteps", 1)) return(1);

	// try the CVDense call
	int NEQ = data->nspec;
	cflag = CVDense(cvode_mem, NEQ);
	if(check_flag(&cflag, "CVDense", 1)) return (1);

	//try the CVLapackDense function instead of the CVBand call
	/* This seems slow and clunky, cvode only really works with RK1 */
	//	int NEQ = data->nspec;
	//	cflag = CVLapackDense(cvode_mem, NEQ);
	//	if(check_flag(&cflag, "CVLapackDense", 1)) return (1);

	// call CVBand to specify the CVBAND band linear solver
	/* This seemed to result in non-convergence, try something else */
	//	int NEQ = data->nspec*data->nspec;
	//	cflag = CVBand(cvode_mem, NEQ, data->nspec, data->nspec);
	//	if(check_flag(&cflag, "CVBand", 1)) return (1);

	// try the CVSpmgr wgucg uses a preconditioner and GMRES
	//	cflag = CVSpgmr(cvode_mem, PREC_NONE, 10);
	//	if(check_flag(&cflag, "CVSpgmr", 1)) return (1);

	//no try the CVSptfqmr methos
	//	cflag = CVSptfqmr(cvode_mem, PREC_NONE, 10);
	//	if(check_flag(&cflag, "CVSptfqmr", 1)) return (1);



	return 0;
}
// retrieve the protected data
template<typename Scalar>
std::vector<Scalar> IgnitionTimes<Scalar>::getYi(){return this->data->Y;}
template<typename Scalar>
Scalar IgnitionTimes<Scalar>::getRho(){return this->data->rho;}
template<typename Scalar>
Scalar IgnitionTimes<Scalar>::getT(){return this->data->T;}
template<typename Scalar>
Scalar IgnitionTimes<Scalar>::getP(){return this->data->P;}

// modify the protected data
template<typename Scalar>
void IgnitionTimes<Scalar>::setYi(std::vector<Scalar> Y_in){this->data->Y = Y_in;}
template<typename Scalar>
void IgnitionTimes<Scalar>::setRho(Scalar rho_in){this->data->rho = rho_in;}
template<typename Scalar>
void IgnitionTimes<Scalar>::setT(Scalar T_in){this->data->T = T_in;}
template<typename Scalar>
void IgnitionTimes<Scalar>::setP(Scalar P_in){this->data->P = P_in;}

template<typename Scalar>
void IgnitionTimes<Scalar>::normalizeYi(std::vector<Scalar>& Yi){
	double sum = 0.;
	int len = Yi.size();
	for(int i=0; i<len; i++){ //sum the species mass fractions and get rid of negative values
		Yi[i] = Yi[i] < 0 ? 0.0 : Yi[i];
		sum += Yi[i];
	}
	sum += tiny; //make sure sum is not 0
	for(int i=0; i<len; i++) //re-normalize to make sure they sum to 1
		Yi[i] /= sum;
	return;
}

template<typename Scalar>
void IgnitionTimes<Scalar>::zeroArray(std::vector<double> &temp){
	for(int i=0;i<temp.size();i++)
		temp[i] = 0.;
}

//TODO setup the solution method that integrates the entire solution in time,
//this would be equivalent to a flow solver time step situation

//return the total enthalpy in J/kg
template<typename Scalar>
double IgnitionTimes<Scalar>::totalEner(double T, std::vector<double> Yi, int nspec, std::vector<double> MW, double R_mix, std::vector<double> Hof, Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double>> *nasa_thermo){
	Antioch::TempCache<double> tc(T);
	std::vector<double> h(nspec);
	nasa_thermo->h(tc,h); // h in J/kg
	double H = 0.;
	for(int i=0; i<nspec; i++)
		H += Yi[i]*(h[i] - Hof[i]);
	double e = H - T*R_mix;
	return e; // total energy in J/kg
}

//return the total enthalpy in J/kg
template<typename Scalar>
double IgnitionTimes<Scalar>::T_from_Ener(double Tguess, double eref, std::vector<double> Yi, int nspec, std::vector<double> MW, std::vector<double> Hof,Antioch::ChemicalMixture<double> *chem_mixture, Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double>> *nasa_thermo){
	int iter_max = 5000;
	//keep the temperature bounded
	Tguess = Tguess<Tmin ? Tmin : Tguess;
	Tguess = Tguess>Tmax ? Tmax : Tguess;

	double T2 {Tguess};
	double e2, F2, Rmix;
	Rmix = chem_mixture->R(Yi);
	e2 = totalEner(T2, Yi, nspec, MW, Rmix, Hof, nasa_thermo);
	F2 = eref - e2;

	if(std::abs(F2) > eps){
		double T0 = T2;
		double F0 = F2;
		double T1 = T0 + 1.;
		double e1 = totalEner(T1, Yi, nspec, MW, Rmix, Hof, nasa_thermo);
		double F1 = eref - e1;
		int iter = 0;

		for(int it = 0; it < iter_max; it++){
			iter++;
			T2 = T1 - F1*(T1-T0)/(F1-F0+tiny);
			T2 = T2<Tmin ? Tmin : T2;
			T2 = T2>Tmax ? Tmax : T2;
			e2 = totalEner(T2, Yi, nspec, MW, Rmix, Hof, nasa_thermo);
			F2 = eref - e2;
			if(std::abs(F2) <= eps)
				break;
			T0 = T1;
			T1 = T2;
			F0 = F1;
			F1 = F2;
		}

		if(iter > iter_max-1){
			std::cout << "Temperature did not converge in T_from_Enth function call!" << std::endl;
			std::cout << "Rethink the structure of the algorithm or problem definition" << std::endl;
			std::cout << "T_err = " << T2 << " F2 = " << F2 << std::endl;
		}

	}

	return T2;
}

template<typename Scalar>
void IgnitionTimes<Scalar>::IgnitionSolver()
{
	double t0 {0.}; // initial time
	double t_end{0.03}; // end time
	double t  {t0}; // current time
	double dt_max {1.e-5}; // maximum allowable time step
	double dt_min {1.e-13}; // minimum allowable time step
	double dt_tmp {dt_max}; // used to determine new dt
	double dt_old {dt_min}; // use to change dt slowly
	double Yi_min {1.e-5}; //only use mass fractions that are chemically significant
	double safety_fac {0.1}; // safety factor for determining dt
	double OH {0.}; // value of OH
	//start dt at dt_min
	dt = dt_min;
	data->OH_trace.push_back(OH);
	data->time.push_back(t);
	std::cout << "time = " << t << " seconds    dt = " << dt << " s    Y_OH = " << OH << "    T = " << data->T << " K    P = " << data->P << " Pa "<< std::endl;
	dt_old = dt;
	while(t<t_end){
		RK_Integrator(dt);
		OH = data->Y[data->i_OH]; //data->map_specToYi["OH"];
		t += dt;
		data->OH_trace.push_back(OH);
		data->time.push_back(t);
		data->P_trace.push_back(data->P);
		data->T_trace.push_back(data->T);
		// output the steps
		std::cout << "time = " << t << " seconds    dt = " << dt << " s    Y_OH = " << OH << "    T = " << data->T << " K    P = " << data->P << " Pa "<< std::endl;
		// check time stepping for species and sum the source term for energy
		double e_dot = 0.;
		dt_tmp = {dt_max};
		for(int i=0; i<data->nspec; i++){
			if(data->Y[i]>Yi_min) dt_tmp = dt_tmp < std::abs((data->rho*data->Y[i]*safety_fac)/(data->w_dot[i] + tiny)) ? dt_tmp : std::abs((data->rho*data->Y[i]*safety_fac)/(data->w_dot[i] + tiny));
			e_dot += data->w_dot[i] * data->Hf[i];
		}
		e_dot *=-1;
		// now check for energy time-stepping
		dt_tmp = dt_tmp < std::abs((data->rho*data->Ener*safety_fac)/(e_dot + tiny)) ? dt_tmp : std::abs((data->rho*data->Ener*safety_fac)/(e_dot + tiny));
		dt_tmp = dt_tmp < dt_min ? dt_min : dt_tmp; // keep dt within the limits (greater than dt_min)
		dt_tmp = dt_tmp > dt_max ? dt_max : dt_tmp; // keep dt within the limits (less than dt_max)
		// vary dt slowly
		dt_tmp = dt_tmp > dt_old ? (0.5*(safety_fac))*(dt_tmp-dt_old) + dt_old : (0.5*(safety_fac))*(dt_old-dt_tmp) + dt_tmp;

		dt_old = dt_tmp;
		dt = dt_tmp;
	}

	N_VDestroy_Serial(y);
	CVodeFree(&cvode_mem);
	return;
}

template<typename Scalar>
void IgnitionTimes<Scalar>::RK_Integrator(Scalar dt)
{
	// setup a 2nd order Runge-Kutta integration step (Ralston method)
	// y_n+1 = y_n + h*(1/4 *k1 + 3/4 * k2)
	// k1 = f(t_n, y_n)
	// k2 = f(t_n+2/3*h, y_n+2/3*h*k1);
	//	int order = 2;
//	std::vector<double> stp_wt  {0.25, 0.75};  // step weight factor
//	std::vector<double> dt_wt   {1./3., 2./3.};// dt weight factor
//	std::vector<Scalar> Y0 = data->Y; // use to store initial mass fraction
//	double E0 = data->Ener;           // store the initial energy
//	std::vector<Scalar> Ytemp (data->nspec, 0.0); // use to store solution at y_n + 2/3*h*k1
//	std::vector<Scalar> w_dot_temp (data->nspec, 0.0); // use to store source term of t_n, y_n
//
//	// solve for k1
//	double k_1 = omega_dot(1.e-9); // a very small, nearly zero time step to get the source terms, but do not divide by 0
//	for(int i=0; i<data->nspec; i++){
//		w_dot_temp[i] = data->w_dot[i];
//		Ytemp[i] = Y0[i] + dt_wt[1]*dt*data->w_dot[i]/data->rho;
//	}
//	normalizeYi(Ytemp);
//	data->Ener = data->Ener + dt_wt[1]*dt*k_1/data->rho;
//	// decode the enthalpy from the energy
//	Scalar R_mix = data->chem_mixture->R(data->Y);
//	data->P = data->rho*(R_mix*data->T); // Pa
//	data->Enthalpy = data->Ener + data->P/data->rho;
//	data->Y = Ytemp;
//	data->T = T_from_Enth(2000.,data->Enthalpy,data->Y, data->nspec, data->MW, data->nasa_thermo);
//
//	// second step: solve for k2
//	double k_2 = omega_dot(dt_wt[1]*dt);
//	for(int i=0; i<data->nspec; i++){
//		data->Y[i] = Y0[i] + dt*(stp_wt[0]*w_dot_temp[i]/data->rho + stp_wt[1]*data->w_dot[i]/data->rho);
//	}
//	data->Ener = E0 + dt*(stp_wt[0]*k_1/data->rho + stp_wt[1]*k_2/data->rho);
//	// decode the enthalpy from the energy
//	R_mix = data->chem_mixture->R(data->Y);
//	data->P = data->rho*(R_mix*data->T); // Pa
//	data->Enthalpy = data->Ener + data->P/data->rho;
//	data->T = T_from_Enth(2000.,data->Enthalpy,data->Y, data->nspec, data->MW, data->nasa_thermo);



	// setup the first order edition
	std::vector<Scalar> Y0 = data->Y; // use to store initial mass fraction
	double E0 = data->Ener;           // store the initial energy
	double k_2 = omega_dot(dt);
	for(int i=0; i<data->nspec; i++){
		data->Y[i] = Y0[i] + dt*(data->w_dot[i]/data->rho);
	}
	data->Ener = E0 + dt*(k_2/data->rho);
	Scalar R_mix = data->chem_mixture->R(data->Y);
	data->T = T_from_Ener(2000.,data->Ener,data->Y, data->nspec, data->MW, data->Hf, data->chem_mixture, data->nasa_thermo);
	data->P = data->rho*(R_mix*data->T); // Pa
	// decode the enthalpy from the energy
	data->Enthalpy = data->Ener + data->P/data->rho;


	return;
}

template<typename Scalar>
double IgnitionTimes<Scalar>::omega_dot(double dt)
{
	std::vector<Scalar> Y_in = data->Y;
	normalizeYi(Y_in); //normalize Yi to 1
	for(int i=0; i<data->nspec; i++)
		NV_Ith_S(y,i) = Y_in[i];
	/* TRY RE-INITIALIZING THE CVODE SOLVER TO GO FROM t=0 to t=dt*/
	//intialize the CVODE solver
	//	int max_steps = 10;
	double t0 = 0.;
	//	cflag = CVodeInit(cvode_mem,f,t0,y);
	//	if(check_flag(&cflag, "CvodeInit", 1)) return(1);
	//
	//	//set the tolerances
	//	cflag = CVodeSStolerances(cvode_mem, reltol, abstol);
	//	if (check_flag(&cflag, "CVodeSStolerances", 1)) return(1);
	//
	//	// set the pointer to the user-defined data
	//	cflag = CVodeSetUserData(cvode_mem, data);
	//	if(check_flag(&cflag, "CVodeSetUserData", 1)) return(1);
	//
	//	//no try the CVSptfqmr methos
	//	cflag = CVSptfqmr(cvode_mem, PREC_NONE, max_steps);
	//	if(check_flag(&cflag, "CVSptfqmr", 1)) return (1);

	double Hf_in = 0., Etot = 0.;
	for(int i=0; i<data->nspec; i++)
		Hf_in += Y_in[i]*data->Hf[i];
	Etot = data->Ener + Hf_in;
	data->Etot = Etot;

	realtype tret = 0.;
	//Reinit
	cflag = CVodeReInit(cvode_mem, t0, y);
	if(check_flag(&cflag, " CVodeReInit", 1)) return(1);
	//	realtype t_curr = data->time[data->time.size()-1];
	CVode(cvode_mem, dt, y, &tret, CV_NORMAL);

	//try and set the initial step to be dt/2
	//	cflag = CVodeSetInitStep(cvode_mem, dt/4.);
	//	if(check_flag(&cflag, "CvodeSetInitStep",1)) return(1);

	//get the number of steps
	//	long int nst = 0;
	//	cflag = CVodeGetNumSteps(cvode_mem, &nst);
	//	if(check_flag(&cflag, "NumberOfSteps", 1)) return(1);
	//	std::cout<< "Cvode completed step in nstp = " << nst << std::endl;

	std::vector<Scalar> Y_out = this->getYi();
	for(int i=0; i<data->nspec; i++)
		Y_out[i] = NV_Ith_S(y,i);
	normalizeYi(Y_out); //normalize Yi to 1

	//find energy source term
	double wdot_x_Hf = 0.; //species mass source term times heat of formation
//	for(int i = 0; i<data->nspec; i++){
//		data->w_dot[i]=data->rho*(Y_out[i]-Y_in[i])/dt;
		//		data->w_dot[i] = NV_Ith_S(y_dot,i)*data->rho;
//		wdot_x_Hf += data->w_dot[i]*data->Hf[i];
//	}
	wdot_x_Hf *=-1.;

	//try the energy source term as: e_dot = rho*(EsenNew - E0)/dt
	double Hf_new = 0.;
	for(int i = 0; i<data->nspec; i++){
			data->w_dot[i]=data->rho*(Y_out[i]-Y_in[i])/dt;
			wdot_x_Hf += data->w_dot[i]*data->Hf[i];
			Hf_new += Y_out[i]*data->Hf[i];
	}
	wdot_x_Hf *=-1.;
	double EsenNew = Etot - Hf_new;
	double e_dot = data->rho*(EsenNew - data->Ener)/dt;

	double e_dot_diff = std::abs(e_dot - wdot_x_Hf);
	std::cout << "Difference in energy source : e_dot_diff = " << e_dot_diff << std::endl;

	return e_dot;
}

template <typename Scalar>
int IgnitionTimes<Scalar>::f(realtype t, N_Vector y, N_Vector y_dot, void *f_data)
{
	double Tguess = 2000.; // use as a standard guess for determining T
	DataForVODE data1 = (DataForVODE) f_data;
	std::vector<double> Yi;
//		double sum = 0.;
	for(int i=0; i<data1->nspec; i++){ //sum the species mass fractions
		Yi.push_back(NV_Ith_S(y,i));
//				Yi[i] = Yi[i] < 0.0 ? 0.0 : Yi[i]; //according to the documentation, try not to delete the negative values during solution time
//				sum += std::abs(Yi[i]);
	}
//		for(int i=0; i<data1->nspec; i++){ //re-normalize to make sure they sum to 1
//			Yi[i] = Yi[i]/sum;
//			NV_Ith_S(y,i) = Yi[i];
//		}


	double Hf_tot = 0.;
	for(int i=0; i<data1->nspec; i++)
		Hf_tot += Yi[i]*data1->Hf[i];
	double e = data1->Etot - Hf_tot;
	double T = T_from_Ener(Tguess, e, Yi, data1->nspec, data1->MW, data1->Hf, data1->chem_mixture, data1->nasa_thermo);
	Antioch::TempCache<Scalar> temp_cache(T);
	//	const Antioch::KineticsConditions<Scalar> conditions(T);

	// initialize molar densities at constant mass/volume (rho)
	zeroArray(data1->_mole_dens);
	data1->chem_mixture->molar_densities(data1->rho,Yi,data1->_mole_dens);

	//initialize thermodynamic variables need for omega_dot
	zeroArray(data1->_H_RT_minus_S);
	data1->nasa_thermo->h_RT_minus_s_R(temp_cache,data1->_H_RT_minus_S);

	//initialize the source term vector and solve for omega_dot
	zeroArray(data1->_omega_dot);
	data1->kinetics->compute_mass_sources(T, data1->_mole_dens, data1->_H_RT_minus_S, data1->_omega_dot);

	for(int i=0; i<data1->nspec;i++)
		NV_Ith_S(y_dot,i) = data1->_omega_dot[i]/data1->rho;

	return 0;
}

template <typename Scalar>
IgnitionTimes<Scalar>::~IgnitionTimes(){}

template <typename Scalar>
int IgnitionTimes<Scalar>::check_flag(void *flagvalue, const char *funcname, int opt)
{
	int *errflag;

	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */

	if (opt == 0 && flagvalue == NULL) {
		fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
				funcname);
		return(1); }

	/* Check if flag < 0 */

	else if (opt == 1) {
		errflag = (int *) flagvalue;
		if (*errflag < 0) {
			fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
					funcname, *errflag);
			return(1); }}

	/* Check if function returned NULL pointer - no memory allocated */

	else if (opt == 2 && flagvalue == NULL) {
		fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
				funcname);
		return(1); }

	return(0);
}


int main()
{
	auto t_start = std::clock();
	std::string mechFile{"share/xml_inputs/GRIMech30.xml"};  // mechanism filename
	std::string fuel{"CH4"}; // fuel name
	double ER {0.4};         // equivalence ratio
	double Patm {3.02};      // atmospheric pressure
	double T {1000./0.68};   // initial temperature
	IgnitionTimes<double>* IT = new IgnitionTimes<double>(mechFile,fuel,ER,Patm,T);
	double T_out = IT->T_from_Ener(2000.,IT->data->Ener,IT->data->Y,IT->data->nspec,IT->data->MW,IT->data->Hf,IT->data->chem_mixture,IT->data->nasa_thermo);
	std::cout << "Testing the T_from_Enth call with initial setup, initial T = " << T << std::endl;
	std::cout << "Output temperature is, T = " << T_out << " K, with an overall error of : " << 100.*std::abs(T_out-T)/T << " %" << std::endl;
	IT->IgnitionSolver();
	/** Output the results **/
	std::cout << "Done, T = " << IT->getT() << " K and P = " << IT->getP() << " Pa " << std::endl;
	auto t_end = std::clock();
	double duration = (t_end-t_start)/(double) CLOCKS_PER_SEC;
	// output the results to a file
	std::ofstream output;
	output.open("OHandPressureAndTempVsTime.dat");
	output << "TITLE=TimePlot \n";
	output << "variables = t Y_OH P T \n";
	for(unsigned int i=0; i<IT->data->time.size(); i++)
		output << IT->data->time[i] <<"\t"<< IT->data->OH_trace[i] <<"\t"<<IT->data->P_trace[i]<<"\t"<<IT->data->T_trace[i]<<"\n";
	output.close();
	std::cout << "Overall computational time is, t = " << duration << " seconds" << std::endl;
	return 0;
}




