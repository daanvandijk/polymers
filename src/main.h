#ifndef _MAIN_
#define _MAIN_
#include <vector>
#include <string>
#include <algorithm>
#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>
#include "cereal_matrix.h"
#include "polymer.h"
#include "time_integration.h"
#include "analysis.h"
#include "experiment.h"

// display help information 
void help();

std::ostream& operator << (std::ostream &os, const data_header &s);

// the interesting functions
void perform_experiments(const char *experiments_path);
void list_experiments(const char *experiments_path);
void perform_experiment(data_header const &exp, polymer_trial &trial, int i);
analysis::ensemble* read_experiment(const char *path);
polymer_trial* read_experiment_raw(const char* path);

std::vector<data_header> get_experiments(const char *experiments_path);
#endif
