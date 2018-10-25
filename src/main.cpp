#include <iostream>
#include <string>
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <thread>
#include <algorithm>
#include <sstream>
#include "main.h"

using namespace std;
using namespace cereal;

std::ostream& operator << (std::ostream &os, const data_header &s) {
  // alpha_p 
  return os << 
    "k = " << s.p.k <<
    ", a = " << s.p.a <<
    ", C = " << s.p.C <<
    ", method = " << s.p.method <<
    " | " << s.p.title;
}

void help() {
  cout << "Polymer simulations\n";
  cout << "~Daan van Dijk 2017-2018\n\n";
  cout << "usage: main [json path] [mode]\n";
  cout << "Modes:\n";
  cout << "  experiment : perform predefined experiment\n";
  cout << "  read [title] : perform analysis of trials\n";
  cout << "  list : list experiments\n";
  cout << "  tex : list experiments in latex table format\n";
}

void perform_experiments(const char *path) {
  vector<data_header> experiments = get_experiments(path);

  // loop through experiments
  int num_threads = std::max((int) thread::hardware_concurrency()-2, 1);
  printf("Number of threads available on system: %i, will spawn %i threads\n", 
      (int) thread::hardware_concurrency(), num_threads);

  // check if experiment needs to be performed...
  auto experiment_left = [&] (const data_header& exp) {
    char filename[255]; 
    sprintf(filename, "%s.dat", exp.p.title.c_str());
    try {
      ifstream file(filename);
      if (file.good()) {
        ifstream is(filename);
        data_header header;
        cereal::PortableBinaryInputArchive archive(is);
        archive(header);

        if (header.p == exp.p) return false;
      }
    }
    catch (const cereal::Exception &e) { printf("Failed to read file\n"); }
    return true;
  };

  vector<data_header> todo(experiments.size());
  auto it = copy_if(experiments.begin(), experiments.end(), todo.begin(), experiment_left);
  todo.resize(distance(todo.begin(), it));

  printf("Experiments left todo: %lu/%lu\n", todo.size(), experiments.size());

  for (auto &exp : todo) {
    printf("%s...\n", exp.p.title.c_str());

    char filename[255]; 
    sprintf(filename, "%s.dat", exp.p.title.c_str());

    // generate polymer
    exp.Ri = MatrixXd::Zero(3, exp.p.Np);
    polymer::generate(exp.p.mu_l, 
        exp.p.sigma_l, 
        exp.p.mu_theta, 
        exp.p.sigma_theta, 
        exp.Ri);
    exp.size = exp.p.Ntrials;

    // save experiment header
    ofstream os(filename);
    cereal::PortableBinaryOutputArchive archive(os);
    archive(exp);

    int i = 0;
    // perform experiments 
    while (i < exp.p.Ntrials) {
      int N = std::min(
          exp.p.Ntrials - i, // number of trials left
          num_threads // number of threads available
          );
      thread *threads = new thread[N];
      polymer_trial *trials = new polymer_trial[N];

      // start threads
      for (int j = 0; j < N; j++) {
        threads[j] = thread(&perform_experiment, ref(exp), ref(trials[j]), j);
        i++;
      }
      printf("Experiment: %i/%i\n", i+1, exp.p.Ntrials);
      // wait till threads are finished
      for (int j = 0; j < N; j++) { 
        threads[j].join();
        // save experiment trial 
        archive(trials[j]);
      }
      delete[] threads;
      delete[] trials;
    }
  }
}

void perform_experiment(data_header const & exp, polymer_trial &trial, int i) {
  trial.R = MatrixXd::Zero(3 * exp.p.datapoints, exp.p.Np);
  trial.time = MatrixXd::Zero(1, exp.p.datapoints);
  time_integration::one_polymer(exp.p, exp.Ri, trial.R, trial.time, i == 0);
}

void list_experiments(const char *path) {
  vector<data_header> experiments = get_experiments(path);

  for (auto &exp : experiments) {
    std::cout << exp << std::endl;
  }
}

polymer_trial* read_experiment_raw(const char* path) {
  //ifstream is(path);
  //data_header header;
  //cereal::PortableBinaryInputArchive archive(is);
  //archive(header);
  //polymer_trial *trial = new polymer_trial;
  //archive(trial);
  //is.close();
  //return trial;
  return NULL;
}

analysis::ensemble* read_experiment(const char *path) {
  ifstream is(path);
  data_header header;
  cereal::PortableBinaryInputArchive archive(is);
  archive(header);
  analysis::ensemble *ensemble = new analysis::ensemble(header.size);

  // output parameter in json format to cout
  stringstream out;
  {
    cereal::JSONOutputArchive arx( out );
    arx(cereal::make_nvp("parameters", header.p)); 
  }
  //std::cout << out.str() << std::endl;

  // We have to read the trials two times.
  // This is because we need to calculate variances,
  // and in order to do so, we need to first calculate mean values. 
  streampos pos = is.tellg();
  for (int i = 0; i < header.size; i++) {
    polymer_trial trial;
    archive(trial);

    // do analysis
    ensemble->add_trial(trial.R, trial.time);
  }
  // read trials again...
  is.seekg(pos);
  for (int i = 0; i < header.size; i++) {
    polymer_trial trial;
    archive(trial);

    // do analysis
    ensemble->add_trial(trial.R, trial.time);
  }
  is.close();

  return ensemble;
}

void tex(const char *path) {
  auto experiments = get_experiments(path);
  cout << "\\begin{table}[h]" << endl;
  cout << "\\centering" << endl;
  cout << "\\footnotesize" << endl;
  cout << "\\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|l|l|l|l|}" << endl;
  cout << "\\hline" << endl;
  cout << "Name & $N_p$ & $N_t$ & $dt$ & $a$ & $k$ & $C$ & $k_BT$ & $l$ & $\\mu_l$ & $\\sigma_l$ & $\\mu_\\theta$ & $\\sigma_\\theta$ & Trials \\\\ \\hline" << endl;

  for (auto &e: experiments) {
    if (e.p.tex_title == "") continue;
    cout.precision(3);
    cout << 
      e.p.tex_title << " & " << 
      ((double) e.p.Np) << " & " << 
      ((double) e.p.Nt) << " & " << 
      e.p.dt << " & " <<
      e.p.a << " & " <<
      e.p.k << " & " <<
      e.p.C << " & " <<
      e.p.k_BT << " & " <<
      e.p.l << " & " <<
      e.p.mu_l << " & " <<
      e.p.sigma_l << " & " <<
      e.p.mu_theta << " & " <<
      e.p.sigma_theta << " & " <<
      e.p.Ntrials << " \\\\" << endl;
  }
  cout << "\\hline" << endl; 
  cout << "\\end{tabular}" << endl;
  cout << "\\caption{Experiment parameters}" << endl;
  cout << "\\label{table:experiment parameters}" << endl;
  cout << "\\end{table}" << endl;
}

int main(int argc, char **argv) {
  if (argc <= 2) 
    help();
  else if (argc > 2) { 
    // path to json file containing experiment parameters
    const char *path = argv[1]; 
    if (strcmp(argv[2], "experiment") == 0) 
      perform_experiments(path);
    else if (strcmp(argv[2], "read") == 0 && argc > 3) {
      auto ensemble = read_experiment(argv[3]);
      free(ensemble);
    }
    else if (strcmp(argv[2], "list") == 0) 
      list_experiments(path);
    else if (strcmp(argv[2], "tex") == 0)
      tex(path);
    else
      help();
  }

  return 0;
}

vector <data_header> get_experiments(const char* path) {
  vector <data_header> headers;
  vector<experiment_parameters> parameters;

  ifstream handle(path);
  cereal::JSONInputArchive archive(handle);
  archive(parameters);
  handle.close();
  for (auto &p: parameters) {
    data_header h;
    h.p = p;
    headers.push_back(h);
  }

  return headers;
}
