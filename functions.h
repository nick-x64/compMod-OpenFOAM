#ifndef FUNC_H
#define FUNC_H

#include <vector>
#include <string>

/*==============.TXT=FILE=PROCESSINF=FUNCTIONS================================*/

std::vector<std::vector<double>> txtParse(std::string const&);

std::vector<std::pair<double, double>> GetMeshTXT(const std::vector<std::vector<double>>&,
                                                  int colNumber = 1);

/*==============.INP=FILE=PROCESSING=FUNCTIONS================================*/

std::vector<std::vector<double>> inpParse(std::string&);

std::vector<std::pair<double, double>> GetMeshINP(const std::vector<std::vector<double>>&,
                                                  double&);

/*==============.CSV=FILE=PROCESSINF=FUNCTIONS================================*/

std::vector<std::vector<double>> csvParse(std::string const&,
                                          const char deliiter = ';');

std::vector<std::pair<double, double>> GetMeshCSV(const std::vector<std::vector<double>>&,
                                                  double&);

std::vector<double> GetTempCSV(std::string&);

/*==============T=FIELD=FILE=PROCESSINF=FUNCTIONS=============================*/

std::vector<double> GetTempField(std::string&);

/*==============CONCENTRATIONS=FILE=PROCESSING================================*/

std::vector<std::pair<double, double>> concParse(std::string& path,
                                                 std::vector<double>& Temperatures,
                                                 const char delimiter=';');

/*==============GRID=PROCESSING=FUNCTIONS=====================================*/

std::vector<std::pair<double, double>> interMesh(std::vector<std::pair<double, double>>&,
                                                 std::vector<std::pair<double, double>>&);

std::vector<std::pair<double, double>> interMesh(std::vector<std::pair<double, double>>&,
                                                 std::vector<double>&);

/*==============MATH=FUNCTIONS=AND=BOOLEAN=FILTERS============================*/

bool isinX(double x,
           std::vector<std::pair<double, double>>&);

bool isinT(double T,
           const std::vector<std::vector<double>>&);

double Planck(double lam,
              double T);

double integrateTrapezoidal(std::vector<std::pair<double, double>>&);

double integrateTrapezoidal(double(&f)(double, double),
                            double,
                            double,
                            double,
                            double);

/*============================================================================*/

std::vector<std::string> list_dir(std::string *);

std::vector<std::string> get_files_list(std::string,
                                        std::string);

std::vector<double> get_conc_from_conc_list(std::vector<std::string>,
                                            std::string pattern=".*conc.1_m3.p_(.+).csv");

std::vector<double> get_pressure_from_params_list(std::vector<std::string>,
                                                 std::string pattern=".*params.p_(.+).csv");

std::vector<double> get_pressure_from_transport_list(std::vector<std::string>,
                                                    std::string pattern=".*transport.p_(.+).txt");

int getElemIndexFromHeaderCSV(std::string&,
                              std::string,
                              const char delimiter=';');

std::vector<double> calculateRateConst(std::string&,
                                       std::string&,
                                       std::string&,
                                       std::vector<double>,
                                       std::string mode="csv");

void output_results(std::vector<std::vector<double>>&,
                    std::vector<double>&,
                    std::vector<double>&,
                    std::string format="sexp");

void output_results_sexp(std::vector<std::vector<double>>&,
                    std::vector<double>&,
                    std::vector<double>&);

std::vector<std::vector<std::vector<double>>> transportPropertiesParser(
        std::vector<std::string>);

void transportPropertiesOutput(std::vector<std::vector<std::vector<double>>>);

std::vector<std::vector<std::vector<double>>> parametersFWBParser(
        std::vector<std::string>,
        const char delimiter=';');


void densityOutput(std::vector<std::vector<std::vector<double>>>);


#endif // FUNC_H
