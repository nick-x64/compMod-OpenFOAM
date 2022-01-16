#include "stdafx.h"
#include "functions.h"
#include <stdexcept>


/*==============.TXT=FILE=PROCESSING=FUNCTIONS================================*/

std::vector<std::vector<double>> txtParse(std::string const& path)
{
    std::fstream tx;

    tx.open(path, std::fstream::in);

    std::vector<std::vector<double>> table;

    try
    {
       if (!tx.is_open())
       {
           throw std::invalid_argument("OpenError: flow.txt file");
       }
       else
       {
//          cout << "Success: .txt file opened" << '\n';

          std::string str;

          std::vector<double> temp;

          while (std::getline(tx, str))
          {
               std::stringstream ss(str);

               while (!ss.eof())
               {
                   ss >> str;

                   temp.push_back(std::stod(str));
               }

               table.push_back(temp);

               temp.clear();

               ss.clear();
          }
       }
    }
    catch (std::invalid_argument& ia)
    {
        std::cerr << ia.what() << '\n';
    }

    tx.close();

//    std::cout << "Success: .txt file closed" << '\n';

    return table;
};

std::vector<std::pair<double, double>> GetMeshTXT(const std::vector<std::vector<double>>& table, int colNumber){

    std::vector<std::pair<double, double>> mesh;

    for (auto& it : table)
    {
        mesh.push_back(std::make_pair(it[0], it[colNumber]));
    }

    return mesh;
}

/*==============.INP=FILE=PROCESSING=FUNCTIONS================================*/

std::vector<std::vector<double>> inpParse(std::string& path_inp) {

    std::ifstream inp(path_inp, std::ifstream::in);

    std::vector<std::vector<double>> table;

    std::string str;

    std::vector<double> NUM; // [0] - NUMBER OF PARAGRAPHS (LINES),
                             // [1] - NUMBER OF VALUES WHITHIN EACH ONE.

    std::vector<double> T, LAM, LINE;

    try {if (!inp.is_open()) {throw std::invalid_argument("OpenError: data.inp file");}

         else {
//            cout << "Success: .inp file opened" << '\n';

/*===============RETRIEVING=FIRST=LINE========================================*/

    for (size_t i = 0; i < 2; ) {

            std::getline(inp, str, ' ');

            try {std::stod(str); NUM.push_back(std::stod(str)); ++i;}

            catch (std::exception& ex) {continue;}
        }

/*===========================================================================*/

    while (table.size() < NUM[0]) {

        int count = 0;

        while (std::getline(inp, str, ' ') && LINE.size() < NUM[1]) {

            try {
                 double val = std::stod(str);
                 ++count;
                 if (count == 2) {T.push_back(val);}
                 if (count > 2) {LINE.push_back(val);}
            }

            catch (std::exception& ex) {continue;}
        }

        if (!std::isnan(LINE[0])) {table.push_back(LINE);}

        LINE.clear();
    }


    while (std::getline(inp, str, ' ')) {

        try {LAM.push_back(std::stod(str));}

        catch (std::exception& ex) {continue;}

    }

    table.push_back(LAM);   // table[table.size() - 2] -- LAM; LAM.size() == 7201

    table.push_back(T);     // table[table.size() - 1] -- T; T.size() == 70

    // Sorting T range, just in case
    std::sort(table[table.size() - 1].begin(),table[table.size() - 1].end());

    // Sorting LAM range, just in case
    std::sort(table[table.size() - 2].begin(),table[table.size() - 2].end());

    }

    } catch (std::invalid_argument& ia) {std::cerr << ia.what() << '\n';}

    inp.close();

//    std::cout << "Success: .inp file closed" << '\n';

    return table;
}

std::vector<std::pair<double, double>> GetMeshINP(const std::vector<std::vector<double>>& table, double& T) {

    std::vector<std::pair<double, double>> mesh;

    //Define a lambda that returns true if the t-value
    //of the given temp range is < the caller's t value
    auto lessThan =
        [](const double temp, double t)
        {return temp < t;};

    try {

        //Find the first table entry whose value is >= caller's T value
        auto iter =
            std::lower_bound(table[table.size() - 1].cbegin(), table[table.size() - 1].cend(), T, lessThan);

        //If the caller's T value is greater than the largest
        //T value in the temp range, we can't interpolate.
        if(iter == table[table.size() - 1].cend()) {
          throw std::invalid_argument("Input T > T_MAX in .inp file");
        }

        //If the caller's T value is less than the smallest
        //T value in the temp range, we can't interpolate.
        if(iter == table[table.size() - 1].cbegin() and T < table[table.size() - 1][0]) {
          throw std::invalid_argument("Input T < T_MIN in .inp file");
        }

        // Interpolation's linear coefficient
        double coef = (T - *(iter - 1)) / (*iter - *(iter - 1));

        // Index of the lower t value in range: table[table.size() - 1][index] < T
        double index = std::distance(table[table.size() - 1].cbegin(), iter) - 1;

        // Fill in new MESH
        for (size_t i = 0; i < table[table.size() - 2].size(); ++i) {

            // if T is not in .INP, but T_MIN <= T <= T_MAX -- interpolate
            if (!isinT(T, table)) {

                double deltaY = table[index + 1][i] - table[index][i];

                if (deltaY >= 0) {mesh.push_back(std::make_pair(
                                  table[table.size() - 2][i],
                                  table[index][i] + coef * deltaY));}

                if (deltaY < 0) {mesh.push_back(std::make_pair(
                                 table[table.size() - 2][i],
                                 table[index][i] + std::pow(coef, -1) * deltaY));}
            } else {

                //if T value coincides with some in .INP -- fill in without interpolation
                //check if coef == 1 -- fill in with the next K(T[index + 1])
                if (coef < 1) {mesh.push_back(std::make_pair(table[table.size() - 2][i],
                                                  table[index][i]));}
                else {mesh.push_back(std::make_pair(table[table.size() - 2][i],
                                     table[index + 1][i]));}
            }
        }
    }
    catch (std::exception& ex) {std::cerr << ex.what() << '\n';}

    // mesh_T[lambda_index].first == LAM[lambda_index]
    // mesh_T[lambda_index].second == K_T[lambda_index]
    return mesh;
}

/*==============.CSV=FILE=PROCESSING=FUNCTIONS================================*/

std::vector<std::vector<double>> csvParse(std::string const& path, const char delimiter)
{
    std::fstream fs;

    fs.open(path, std::fstream::in);

    std::vector<std::vector<double>> table;

    std::vector<double> _TEMP;

    try
    {
       if (!fs.is_open())
       {
           throw std::invalid_argument("OpenError: .csv file");
       }
       else
       {
//          cout << "Success: .csv file opened" << '\n';

          std::string str, header;

          std::vector<double> temp;

          std::stringstream ss;

          std::getline(fs, header, '\n');

          _TEMP = GetTempCSV(header);

          while (std::getline(fs, str, '\n'))
          {
               ss << str;

               while (std::getline(ss, str, delimiter))
               {
                    temp.push_back(std::stod(*(&str)));
               }

               table.push_back(temp);

               temp.clear();

               ss.clear();
          }

          for (auto& row : table){*begin(row) *= 1e-9;}
       }
    }
    catch (std::invalid_argument& ia)
    {
        std::cerr << ia.what() << '\n';
    }

    fs.close();

//    std::cout << "Success: .csv file closed" << '\n';

    table.push_back(_TEMP);


    // K_ABS[lambda][T] == table[lambda][T] [1/cm]
    // T == table[table.size() - 1] [K]
    return table;
};

std::vector<std::pair<double, double>> GetMeshCSV(const std::vector<std::vector<double>>& table, double& T){

    std::vector<std::pair<double, double>> mesh;

    //Define a lambda that returns true if the t-value
    //of the given temp range is < the caller's t value
    auto lessThan =
        [](const double temp, double t)
        {return temp < t;};

    try {

        //Find the first table entry whose value is >= caller's T value
        auto iter =
            std::lower_bound(table[table.size() - 1].cbegin(), table[table.size() - 1].cend(), T, lessThan);

        //If the caller's T value is greater than the largest
        //T value in the temp range, we can't interpolate.
        if(iter == table[table.size() - 1].cend()) {
          throw std::invalid_argument("Input T > T_MAX in .csv file");
        }

        //If the caller's T value is less than the smallest
        //T value in the temp range, we can't interpolate.
        if(iter == table[table.size() - 1].cbegin() and T < table[table.size() - 1][0]) {
          throw std::invalid_argument("Input T < T_MIN in .csv file");
        }

        // Interpolation's linear coefficient
        double coef = (T - *(iter - 1)) / (*iter - *(iter - 1));

        // Index of the lower t value in range: table[table.size() - 1][index] < T
        double index = std::distance(table[table.size() - 1].cbegin(), iter);

        // Fill in new MESH
        for (size_t i = 0; i < table.size() - 1; ++i) {

            // if T is not in .csv, but T_MIN <= T <= T_MAX -- interpolate
            if (!isinT(T, table)) {

                double deltaY = table[i][index + 1] - table[i][index];

                if (deltaY >= 0) {mesh.push_back(std::make_pair(
                                  table[i][0],
                                  100 * (table[i][index] + coef * deltaY)));}

                if (deltaY < 0) {mesh.push_back(std::make_pair(
                                 table[i][0],
                                 100 * (table[i][index] + std::pow(coef, -1) * deltaY)));}
            } else {

                //if T value coincides with some in .csv -- fill in without interpolation
                //check if coef == 1 -- fill in with the next K(T[index + 1])
                if (coef < 1) {mesh.push_back(std::make_pair(table[i][0],
                                                  100 * table[i][index]));}
                else {mesh.push_back(std::make_pair(table[i][0],
                                     100 * table[i][index + 1]));}
            }
        }
    }
    catch (std::exception& ex) {std::cerr << ex.what() << '\n';}

    // mesh_T[lambda_index].first == LAM[lambda_index]
    // mesh_T[lambda_index].second == K_T[lambda_index]
    return mesh;
}

/**
 * Retrieve the temperatures' values from .csv header (with respect to RegEx used)
 *
 * @brief GetColIndForTempCSVHead
 * @param str Header line to be parsed
 * @returns std::vector<double> Temperature values
 */
std::vector<double> GetTempCSV(std::string& str){

    std::regex reg("((\\w+,\\w+\\\\\\w+,\\w+);|T=((\\d+|\\d+\\.\\d*));?)");

    std::vector<double> Temp;

    auto firstIt = std::sregex_iterator(std::begin(str), std::end(str), reg);

    auto lastIt = std::sregex_iterator();

    try{

    for (std::sregex_iterator k = firstIt; k != lastIt; ++k) {

        std::smatch match = *k;

        try {Temp.push_back(std::stod(match[match.size() - 1]));}

        catch (...) {continue;}

//        try {if (std::stod(match[match.size() - 1]) == temp)
//            {return std::distance(firstIt, k);}
//            else{continue;}
//        } catch (...) {continue;}
        }
    }

    catch (std::exception& e){std::cerr << e.what() << '\n';}

    return Temp;
};

/*==============T=FIELD=FILE=PROCESSING=FUNCTIONS=============================*/

/**
 * Retrieve temperature values from OpenFOAM "T" file.
 *
 * @brief GetTempField
 * @param pathTemp Absolute PATH to "T" file.
 * @return std::vector<double> Temperature values
 */
std::vector<double> GetTempField(std::string& pathTemp)
{
std::ifstream Temperatures(pathTemp, std::ios_base::in);

std::vector<double> tempVal;

try
    {
       if (!Temperatures.is_open())
       {
           throw invalid_argument("OpenError: T field file");
       }
       else
       {
//          cout << "Success: temp file opened" << '\n';

          std::string str;

          while (str != "internalField   nonuniform List<scalar> ")
          {
              getline(Temperatures, str, '\n');
          }

          getline(Temperatures, str, '\n');

          int size = std::stoi(str);

          getline(Temperatures, str, '\n');

          for (int i = 0; i < size; ++i)
          {
              getline(Temperatures, str, '\n');

              tempVal.push_back(std::stod(str));
          }
       }
    }
    catch (invalid_argument& ia)
    {
        cerr << ia.what() << endl;
    }

    Temperatures.close();

//    std::cout << "Success: temp file closed" << '\n';

    return tempVal;
}

/*==============CONCENTRATIONS=FILE=PROCESSING================================*/

/**
  Retrieveing S2 concentrations from FWB file "con.1_m3.p_(...).csv"
 * @brief conсentraionParse
 * @param path - path to the contcentraions' file from FWB's ter-reactor's output
 * @return std::vector<std::pair<double, double>> which contains concentrations according to the T field's temperature
 */
std::vector<std::pair<double, double>> concParse(std::string& path, std::vector<double>& Temperatures, const char delimiter) {

    std::fstream cs(path, std::fstream::in);

    //S2 concentrations range
    std::vector<std::pair<double, double>> conc;

    try
    {
       if (!cs.is_open())
       {
           throw std::invalid_argument("OpenError: con.csv file");
       }
       else
       {
    //          cout << "Success: con.csv file opened" << '\n';

          std::string str;

          std::vector<double> temp;

          std::stringstream ss;

          std::getline(cs, str, '\n');

          int elem_index = getElemIndexFromHeaderCSV(str, "S2");

          while (std::getline(cs, str, '\n'))
          {
               ss << str;

               while (std::getline(ss, str, delimiter))
               {
                    temp.push_back(std::stod(*(&str)));
               }

               conc.push_back(std::make_pair(temp[0], temp[elem_index]));

               temp.clear();

               ss.clear();
          }
       }
    }
    catch (std::invalid_argument& ia)
    {
        std::cerr << ia.what() << '\n';
    }

    cs.close();

    //    std::cout << "Success: con.csv file closed" << '\n';

    std::vector<std::pair<double, double>> concMesh = interMesh(conc, Temperatures);

    concMesh.erase(std::unique(concMesh.begin(), concMesh.end(),
             [](const std::pair<double, double>& a, const std::pair<double, double>& b)
             {return std::abs(a.first - b.first) < 1e-11;}), concMesh.end());

    return concMesh;
}


/*==============GRID=PROCESSING=FUNCTIONS=====================================*/

std::vector<std::pair<double, double>> interMesh(std::vector<std::pair<double, double>>& m1,
                                                 std::vector<std::pair<double, double>>& m2) {
    Interpolator grid1(m1);

    std::vector<std::pair<double, double>> m1n;

    for (auto& it: m1) {m1n.push_back(it);}

    for (size_t i = 0; i < m2.size(); ++i) {

        double xn = m2[i].first;

        if (!isinX(xn, m1)) {m1n.push_back(std::make_pair(xn, grid1.findValue(xn)));}

    }

    std::sort(std::begin(m1n), std::end(m1n),
              [](std::pair<double, double> a, std::pair<double, double> b) {
        return (b.first - a.first) > 1e-11;
    });

    return m1n;
}

std::vector<std::pair<double, double>> interMesh(std::vector<std::pair<double, double>>& m1,
                                                 std::vector<double>& m2) {
    Interpolator grid1(m1);

    std::vector<std::pair<double, double>> m1n;

    for (auto& it: m1) {m1n.push_back(it);}

    for (size_t i = 0; i < m2.size(); ++i) {

        double xn = m2[i];

        if (!isinX(xn, m1)) {m1n.push_back(std::make_pair(xn, grid1.findValue(xn)));}

    }

    std::sort(std::begin(m1n), std::end(m1n),
              [](std::pair<double, double> a, std::pair<double, double> b) {
        return (b.first - a.first) > 1e-11;
    });

    return m1n;
}

/*==============MATH=FUNCTIONS=AND=BOOLEAN=FILTERS============================*/

bool isinX(double x, std::vector<std::pair<double, double>>& mesh) {

    std::vector<double> valuesX;

    for (auto& it: mesh) {valuesX.push_back(it.first);};

    if (std::find(std::begin(valuesX), std::end(valuesX), x) != std::end(valuesX)){return true;}

    return false;
}

bool isinT(double T, const std::vector<std::vector<double>>& table) {

    int count = 0;

    for (auto& it: table[table.size() - 1]) {if (it == T) {++count;}}

    if (count > 0) {return true;}

    return false;
}

double Planck(double lam, double T)
{
    return 2 * 6.626e-34 * std::pow(3e8, 2) * std::pow(lam, -5) /
            std::expm1(6.626e-34 * 3e8 * std::pow(lam * T * 1.38e-23, -1));
};

double integrateTrapezoidal(std::vector<std::pair<double, double>>& F)
{
    double sum = 0;

    for (std::size_t i = 1; i < F.size(); ++i)
    {
        sum += ((F[i - 1].second + F[i].second) * (F[i].first - F[i - 1].first));

    }
    return sum / 2;
}

double integrateTrapezoidal(double(&f)(double, double), double min, double max, double N, double T){

    double h = (max - min) / N;

    double sum = f(min, T) + f(max, T);

    for (int i = 1; i < N; ++i){sum += 2 * f(min + i * h, T);}

    return sum * (h / 2);
}

/*============================================================================*/

std::vector<std::string> list_dir(std::string *path) {
   struct dirent *entry;
   std::vector<std::string> ent_list;
   DIR *dir = opendir(path->c_str());

   if (dir == NULL) {
      return ent_list;
   }
   while ((entry = readdir(dir)) != NULL) {
//   cout << entry->d_name << endl;
   ent_list.push_back((entry->d_name));
   }
   closedir(dir);
   return ent_list;
}


std::vector<std::string> get_files_list(
        std::string input_dir,
        std::string mask){

    std::regex reg(mask);

    std::vector<std::string> files_list;

    std::vector<std::string> input_entries = list_dir(&input_dir);

    for (size_t i = 0; i < input_entries.size(); i++){

        if (std::regex_search(input_entries[i], reg)){

            files_list.push_back(input_dir + "\\" + input_entries[i]);
        }
    }
    return files_list;
}

std::vector<double> get_conc_from_conc_list(
        std::vector<std::string> conc_files,
        std::string pattern){

    std::regex re(pattern);

    std::smatch match;

    std::vector<double> conc;

    for (size_t i = 0; i < conc_files.size(); i++){
        std::regex_search(conc_files[i], match, re);
        conc.push_back(std::stod(match[match.size() - 1]));}

    return conc;
}


std::vector<double> get_pressure_from_params_list(
        std::vector<std::string> params_files,
        std::string pattern){

    std::regex re(pattern);

    std::smatch match;

    std::vector<double> pressure;

    for (size_t i = 0; i < params_files.size(); i++){
        std::regex_search(params_files[i], match, re);
        pressure.push_back(std::stod(match[match.size() - 1]));}

    return pressure;
}


std::vector<double> get_pressure_from_transport_list(
        std::vector<std::string> params_files,
        std::string pattern){

    std::regex re(pattern);

    std::smatch match;

    std::vector<double> pressure;

    for (size_t i = 0; i < params_files.size(); i++){
        std::regex_search(params_files[i], match, re);
        pressure.push_back(std::stod(match[match.size() - 1]));}

    return pressure;
}


int getElemIndexFromHeaderCSV(std::string& line,
                              std::string elem,
                              const char delimiter){

    std::stringstream ss;

    std::string str;

    int count = -1;

    ss << line;

    while(std::getline(ss, str, delimiter)){

        count++;

        if (str == elem){return count;}
    }
    return count;
}

std::vector<double> calculateRateConst(
        std::string& conc_files,
        std::string& path_csv,
        std::string& path_txt,
        std::vector<double> Temperatures,
        std::string mode){

    std::vector<std::pair<double, double>> ConcentrationS2 = concParse(conc_files, Temperatures);

    std::vector<std::vector<double>> kabs;

    std::vector<std::pair<double, double>> mesh2;

    std::vector<double> column;

    std::cout << "Processing:\t" << path_csv << '\n';

    std::cout << "Processing:\t" << conc_files << '\n';

//    for (size_t i = 1; i < ConcentrationS2.size(); ++i){
//        std::cout << ConcentrationS2[i].first << '\t' << ConcentrationS2[i].second << '\n';}

    if (mode == "inp") {

        thread t1([&path_csv, &kabs]() {kabs = inpParse(path_csv);});

        thread t2([&path_txt, &mesh2]() {mesh2 = GetMeshTXT(txtParse(path_txt));});

        t1.join();

        t2.join();
    }

    if (mode == "csv") {

        thread t1([&path_csv, &kabs]() {kabs = csvParse(path_csv);});

        thread t2([&path_txt, &mesh2]() {mesh2 = GetMeshTXT(txtParse(path_txt));});

        t1.join();

        t2.join();
    }

/*==============ITERATING=OVER=TEMPERATURE=RANGE==============================*/

    for (size_t i = 1; i < Temperatures.size(); ++i){

        double T = Temperatures[i];

        vector<pair<double, double>> mesh1, mesh1new, mesh2new;


    /*=============GETTING=MESH=FROM=KABS=========================================*/

        if (mode == "inp") {mesh1 = GetMeshINP(kabs, T);}

        if (mode == "csv") {mesh1 = GetMeshCSV(kabs, T);}

    /*=============INTERPOLATION==================================================*/

        mesh1new = interMesh(mesh1, mesh2);

        mesh2new = interMesh(mesh2, mesh1);

    /*=====================ERASING=REDUNDANT=POINTS===============================*/

        mesh1new.erase(std::unique(mesh1new.begin(), mesh1new.end(),
                 [](const std::pair<double, double>& a, const std::pair<double, double>& b)
                 {return std::abs(a.first - b.first) < 1e-11;}), mesh1new.end());

        mesh2new.erase(std::unique(mesh2new.begin(), mesh2new.end(),
                 [](const std::pair<double, double>& a, const std::pair<double, double>& b)
                 {return std::abs(a.first - b.first) < 1e-11;}), mesh2new.end());

    /*=========================SORTING=CHECK======================================*/

        if (!std::is_sorted(std::begin(mesh1new), std::end(mesh1new),
                       [](std::pair<double, double>& a, std::pair<double, double>& b)
                       {return a.first < b.first;}))
                  {std::sort(std::begin(mesh1new), std::end(mesh1new),
                  [](std::pair<double, double>& a, std::pair<double, double>& b)
                  {return a.first < b.first;});}

        if (!std::is_sorted(std::begin(mesh2new), std::end(mesh2new),
                            [](std::pair<double, double>& a, std::pair<double, double>& b)
                            {return a.first < b.first;}))
                  {std::sort(std::begin(mesh2new), std::end(mesh2new),
                  [](std::pair<double, double>& a, std::pair<double, double>& b)
                  {return a.first < b.first;});}

    /*==========================CALCULATION=======================================*/

        std::vector<std::pair<double, double>> func;

        for (size_t i = 0; i < mesh1.size(); ++i) {

            double mlt = mesh1[i].second * mesh2new[i].second *
                         mesh1[i].first * 5.0307e33;
            if (mlt > 1e-10 && !std::isnan(mlt)) {func.push_back(
                            std::make_pair(mesh1[i].first, mlt));}
        }

        func.pop_back();

    /*=====================ERASING=REDUNDANT=POINTS===============================*/

        func.erase(std::unique(func.begin(), func.end(),
                 [](const std::pair<double, double>& a, const std::pair<double, double>& b)
                 {return std::abs(a.first - b.first) < 1e-11;}), func.end());

    /*=======================INTEGRATION==========================================*/

        double result = 4 * M_PI * integrateTrapezoidal(func) / ConcentrationS2[i].second;

        column.push_back(result);

//        std::cout << result << '\n';

    }

    return column;
}

void output_results(std::vector<std::vector<double>>& output_,
                    std::vector<double>& conc,
                    std::vector<double>& Temperatures,
                    std::string format){

    /*=================CREATING=AN=OUTPUT=FILE====================================*/

            try {

                if (format == "csv"){

                    std::ofstream data;

                    data.open("..\\output\\photorate.dat", std::ios_base::out);

                    std::string header = "T[K]/P[atm];";

                    for (size_t i = 0; i < conc.size(); i++){header += std::to_string(conc[i]) + ";";}

                    data << header << '\n';

                    for (size_t i = 0; i < Temperatures.size(); i++){

                        data << Temperatures[i] << ';';

                        for(size_t j = 0; j < conc.size(); j++){data << std::to_string(output_[j][i]) + ";";}

                        data << "\n";

                        data.close();   //closing output file
                    }
                } else if (format == "sexp"){

                    std::ofstream data;

                    data.open("..\\output\\photorate.dat", std::ios_base::out);

                    // change output format settings with member functions
                    data.setf(std::ios::scientific, std::ios::floatfield); // set fixed floating format

                    data << "(\n";

                    for (size_t i = 0; i < conc.size(); i++){

                        data.precision(4); // for fixed format, two decimal places

                        data << '(' << std::setw(11) << std::uppercase
                             << conc[i] * 101325 << " " << "((" << std::setprecision(6) << Temperatures[0]
                             << "  "  << output_[i][0] << ")\n";

                        for(size_t j = 1; j < Temperatures.size() - 1; j++){

                            data << "\t\t\t  (" << Temperatures[j]
                                 << "  " << output_[i][j] << ")\n";

                        }
                        data << "\t\t\t  (" << Temperatures[Temperatures.size() - 1]
                             << "  " << output_[i][Temperatures.size() - 1] << ")))\n";
                    }

                    data << "\n)";

                    data.close();   //closing output file

                } else {throw invalid_argument("Invalid output format. Choose from: sexp, csv.");}

            } catch (invalid_argument& ia){cerr << ia.what() << endl;}

}

void output_results_sexp(std::vector<std::vector<double>>& output_,
                    std::vector<double>& conc,
                    std::vector<double>& Temperatures){

    /*=================CREATING=AN=OUTPUT=FILE====================================*/

            std::ofstream data;

            data.open("..\\output\\phptorate.dat", std::ios_base::out);

//            // save the current settings
//            ios::fmtflags old_settings = data.flags(); //save previous format flags
//            int old_precision = data.precision(); // save previous precision setting


            // change output format settings with member functions
            data.setf(std::ios::scientific, std::ios::floatfield); // set fixed floating format
            data.precision(2); // for fixed format, two decimal places

            data << "(\n";

            for (size_t i = 0; i < conc.size(); i++){

                data << '(' << std::setw(9) << conc[i] << " " << "((" << Temperatures[0] << "  " << output_[i][0] << ")\n";

                for(size_t j = 1; j < Temperatures.size() - 1; j++){

                    data << "\t\t\t(" << std::setw(9) << Temperatures[j] << "  " << output_[i][j] << ")\n";

                }
                data << "\t\t\t(" << Temperatures[Temperatures.size() - 1] << "  " << output_[i][Temperatures.size() - 1] << ")))\n";
            }

            data << "\n)";

            data.close();   //closing output file

}

std::vector<std::vector<std::vector<double>>> transportPropertiesParser(
        std::vector<std::string> transport_files){

    std::vector<std::vector<std::vector<double>>> tables;

    for (auto path: transport_files){

        std::fstream tr;

        tr.open(path, std::fstream::in);

        std::vector<std::vector<double>> table;

        try
        {
           if (!tr.is_open())
           {
               throw std::invalid_argument("OpenError: transport.txt file");
           }
           else
           {
              std::string str, header;

              std::vector<double> temp;

              std::getline(tr, header);

              while (std::getline(tr, str))
              {


                   std::stringstream ss(str);

                   while (!ss.eof())
                   {
                       ss >> str;

                       temp.push_back(std::stod(str));
                   }


                   if (temp.size() == 9){table.push_back(temp);}


                   temp.clear();

                   ss.clear();
              }

           }
        }
        catch (std::invalid_argument& ia)
        {
            std::cerr << ia.what() << '\n';
        }

        tr.close();

        tables.push_back(table);


    // Print the contents
    //    for (size_t i = 0; i < table.size(); i++){
    //        for (size_t j = 0; j < table[i].size(); j++){
    //            std::cout << table[i][j] <<  '\t';
    //        }
    //        std::cout << '\n';
    //    }

    }

    return tables;

}


void transportPropertiesOutput(
        std::vector<std::vector<std::vector<double>>> transportTables){

    // Creating Enthalpy file
            std::ofstream data;

            data.open("..\\output\\EnthalpyTable", std::ios_base::out);

            // change output format settings with member functions
            data.setf(std::ios::scientific, std::ios::floatfield); // set fixed floating format

            data << "(\n";

            for (size_t i = 0; i < transportTables[0].size(); i++){

                data.precision(4); // for fixed format, two decimal places

                data << '(' << std::setw(11) << std::uppercase
                     << transportTables[0][i][0] << " ( ";

                for (size_t j = 0; j < transportTables.size(); j++){
                    data << '(' << std::setprecision(6) << transportTables[j][i][1] * 101325
                         << ' '  << transportTables[j][i][2] * 28.0134 << ") ";
                }
                data << "))\n";
            }

            data << ")";

            data.close();   //closing output file


    // Creating HeatCapacity file
            data.open("..\\output\\HeatCapacityTable", std::ios_base::out);

            // change output format settings with member functions
            data.setf(std::ios::scientific, std::ios::floatfield); // set fixed floating format

            data << "(\n";

            for (size_t i = 0; i < transportTables[0].size(); i++){

                data.precision(4); // for fixed format, two decimal places

                data << '(' << std::setw(11) << std::uppercase
                     << transportTables[0][i][0] << " ( ";

                for (size_t j = 0; j < transportTables.size(); j++){
                    data << '(' << std::setprecision(6) << transportTables[j][i][1] * 101325
                         << ' '  << transportTables[j][i][3] * 28.0134 / 1000 << ") ";
                }
                data << "))\n";
            }

            data << ")";

            data.close();   //closing output file


    // Create MolecularWeight file
            data.open("..\\output\\MolecularWeightTable", std::ios_base::out);

            // change output format settings with member functions
            data.setf(std::ios::scientific, std::ios::floatfield); // set fixed floating format

            data << "(\n";

            for (size_t i = 0; i < transportTables[0].size(); i++){

                data.precision(4); // for fixed format, two decimal places

                data << '(' << std::setw(11) << std::uppercase
                     << transportTables[0][i][0] << " ( ";

                for (size_t j = 0; j < transportTables.size(); j++){
                    data << '(' << std::setprecision(6) << transportTables[j][i][1] * 101325
                         << ' '  << transportTables[j][i][7] << ") ";
                }
                data << "))\n";
            }

            data << ")";

            data.close();   //closing output file


    // Create Viscosity file
            data.open("..\\output\\ViscosityTable", std::ios_base::out);

            // change output format settings with member functions
            data.setf(std::ios::scientific, std::ios::floatfield); // set fixed floating format

            data << "(\n";

            for (size_t i = 0; i < transportTables[0].size(); i++){

                data.precision(4); // for fixed format, two decimal places

                data << '(' << std::setw(11) << std::uppercase
                     << transportTables[0][i][0] << " ( ";

                for (size_t j = 0; j < transportTables.size(); j++){
                    data << '(' << std::setprecision(6) << transportTables[j][i][1] * 101325
                         << ' '  << transportTables[j][i][4] << ") ";
                }
                data << "))\n";
            }

            data << ")";

            data.close();   //closing output file


    // Create ThermalConductivity file
            data.open("..\\output\\ThermalConductivityTable", std::ios_base::out);

            // change output format settings with member functions
            data.setf(std::ios::scientific, std::ios::floatfield); // set fixed floating format

            data << "(\n";

            for (size_t i = 0; i < transportTables[0].size(); i++){

                data.precision(4); // for fixed format, two decimal places

                data << '(' << std::setw(11) << std::uppercase
                     << transportTables[0][i][0] << " ( ";

                for (size_t j = 0; j < transportTables.size(); j++){
                    data << '(' << std::setprecision(6) << transportTables[j][i][1] * 101325
                         << ' '  << transportTables[j][i][5] << ") ";
                }
                data << "))\n";
            }

            data << ")";

            data.close();   //closing output file
}


std::vector<std::vector<std::vector<double>>> parametersFWBParser(
        std::vector<std::string> params_files,
        const char delimiter){

    std::vector<std::vector<std::vector<double>>> tables;

    for (auto path: params_files){

        std::fstream fs;

        fs.open(path, std::fstream::in);

        std::vector<std::vector<double>> table;

        try
        {
           if (!fs.is_open())
           {
               throw std::invalid_argument("OpenError: (parameters from FWB).csv file");
           }
           else
           {
              std::string str, header;

              std::vector<double> temp;

              std::stringstream ss;

              std::getline(fs, header, '\n');

              while (std::getline(fs, str, '\n'))
              {
                   ss << str;

                   while (std::getline(ss, str, delimiter))
                   {
                        temp.push_back(std::stod(*(&str)));
                   }

                   table.push_back(temp);

                   temp.clear();

                   ss.clear();
              }
           }
        }
        catch (std::invalid_argument& ia)
        {
            std::cerr << ia.what() << '\n';
        }

        fs.close();

        tables.push_back(table);

     }

    return tables;
}


void densityOutput(std::vector<std::vector<std::vector<double>>> parametersTables){


            std::ofstream data;

            data.open("..\\output\\DensityTable", std::ios_base::out);

            // change output format settings with member functions
            data.setf(std::ios::scientific, std::ios::floatfield); // set fixed floating format

            data << "(\n";

            for (size_t i = 0; i < parametersTables[0].size(); i++){

                data.precision(4); // for fixed format, two decimal places

                data << '(' << std::setw(11) << std::uppercase
                     << parametersTables[0][i][0] << " ( ";

                for (size_t j = 0; j < parametersTables.size(); j++){
                    data  << '(' << std::setprecision(6) << parametersTables[j][i][2] * 101325
                     << ' '  << parametersTables[j][i][4] << ") ";
                }
                data << "))\n";
            }

            data << ")";

            data.close();   //closing output file
}






/*===========FLOW=FILE=GENERATING=============================================*/

//    double te = 5000;

//    std::ofstream flow("..\\output\\flow_my.txt", std::ios_base::out);

//    std::vector<std::vector<double>> table_csv = csvParse("..\\input\\absS2_my.csv", ';');

//    std::vector<std::pair<double, double>> m_csv = GetMeshCSV(table_csv, te);

//    for (size_t i = 20; i < 801; ++i) {

//        flow << " " << i * 1e-9 << " " << Planck(i * 1e-9, te) / 1e-9
//             << " " << NAN << "  " << NAN << '\n';
//    }

//    flow.close();


//    std::vector<double> lam;

//    for (int i = 0; i < 781; ++i) {lam.push_back((20 + i) * 1e-9);}

//        std::ofstream _flow;

//        _flow.open("..\\output\\flow_sample_" + std::to_string(5000) + ".txt", std::ios_base::out);

//        for (int j = 0; j < 781; ++j) {

//            double lam = (20 + j) * 1e-9;

//            double _planck = Planck(lam, 5000);

//            if (_planck <= 1e-45) {_flow << ' ' << lam << ' ' << 0 << ' ' << NAN << ' ' << NAN << '\n';}

//            else {_flow << ' ' << lam << ' ' << _planck << ' ' << NAN << ' ' << NAN << '\n';}
//        }

//        _flow.close();


///*============================================================================*/


////    std::ofstream _flow;

////    _flow.open("..\\output\\flow_sample.txt", std::ios_base::out);

////    for (int i = 0; i < 781; ++i) {

////        double lam = (20 + i) * 1e-9;

////        double _planck = Planck(lam, 1323);

////        if (_planck < 1e-45){_flow << ' ' << lam << ' ' << 0 << ' ' << NAN << ' ' << NAN << '\n';}

////        else{_flow << ' ' << lam << ' ' << _planck << ' ' << NAN << ' ' << NAN << '\n';}
