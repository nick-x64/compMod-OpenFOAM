#include "stdafx.h"
#include "functions.h"


std::vector<std::string> abs_files = get_files_list("..\\input",
                                                    "abs.p(?=.*.csv)");

std::vector<std::string> conc_files = get_files_list("..\\input",
                                                     "conc.1_m3(?=.*.csv)");

std::string path_txt = "..\\input\\flow_sample.txt";

std::string path_tmp = "..\\input\\T_train_1";

std::vector<double> conc = get_conc_from_conc_list(conc_files);

std::vector<double> Temperatures = GetTempField(path_tmp);

std::vector<std::vector<double>> output_;

//std::vector<std::string> transport_files = get_files_list("..\\input",
//                                                     "transport.p_(?=.*.txt)");

//std::string path_transport = "..\\input\\transport.txt";

//std::vector<std::string> params_files = get_files_list("..\\input",
//                                                     "params.p_(?=.*.csv)");

//std::string path_parameters = "..\\input\\parameters.csv";

//std::vector<double> pressures_params = get_pressure_from_params_list(params_files);

//std::vector<double> pressures_transport = get_pressure_from_transport_list(transport_files);


int main()
{
    Timer t; //time counter

    for (size_t j = 0; j < conc_files.size(); j++){

        output_.push_back(calculateRateConst(conc_files[j],
                                             abs_files[j],
                                             path_txt,
                                             Temperatures,
                                             "csv"));
    }

/*=================CREATING=AN=OUTPUT=FILE====================================*/

        output_results(output_, conc, Temperatures, "sexp");

/*=================TEST=OUTPUT=FUNCTION=======================================*/

//        std::vector<std::vector<double>> out;

//        std::vector<double> con, Temp;

//        con.push_back(0.1);
//        con.push_back(1);
//        con.push_back(2);

//        for (int i = 1323; i <= 1500; i++){Temp.push_back(double(i));}

//        for (size_t i = 0; i < con.size(); i++){
//            std::vector<double> col;
//            for (size_t j = 0; j < Temp.size(); j++){col.push_back(double(42));}
//            out.push_back(col);
//        }

//        output_results(out, con, Temp, "sexp");

/*============================================================================*/

//  Parse the FWB properties output and write them as s-exp. into files

//    transportPropertiesOutput(transportPropertiesParser(transport_files));
//    densityOutput(parametersFWBParser(params_files));

//  Check possible pressure values from files' names

//    for (size_t i = 0; i < pressures_params.size(); i++){std::cout << pressures_params[i] << '\n';}

//    for (size_t i = 0; i < pressures_transport.size(); i++){std::cout << pressures_transport[i] << '\n';}


/*============================================================================*/

//        for (int j = 1300; j <= 1500; j++){std::cout << j << '\n';}

//        for (size_t j = 0; j < output_[0].size(); j++){
//            std::string line;
//            for (size_t i = 0; i < output_.size(); i++){
//                line +=  std::to_string(output_[i][j]);
//                line += '\t';
//           }
//            std::cout << line << '\n';
//        }

    return 0;
};


//         RUNTIME | len(Temperatures) == 240000: 58403.1 seconds
//         RUNTUME | len(Temperatures) == 240000: 63470   seconds
//         RUNTUME | len(Temperatures) == 480000: 115852   seconds
