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
    " | " << s.p.title;
}

void help() {
  cout << "Polymer simulations\n";
  cout << "~Daan van Dijk 2017-2018\n\n";
  cout << "usage: main [mode]\n";
  cout << "Modes:\n";
  cout << "  experiment : perform predefined experiment\n";
  cout << "  read [title] : perform analysis of trials\n";
  cout << "  list : list experiments\n";
  cout << "  tex : list experiments in latex table format\n";
}

void perform_experiments() {
  vector<data_header> experiments = get_experiments();

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
  time_integration::one_polymer(exp.p.Nt, 
      exp.p.dt, 
      exp.p.a, 
      exp.p.k, 
      exp.p.l, 
      exp.p.k_BT, 
      exp.p.C,
      exp.Ri, 
      exp.p.datapoints, 
      trial.R, 
      trial.time,
      exp.p.space,
      i == 0); 
}

void list_experiments() {
  vector<data_header> experiments = get_experiments();

  ofstream handle("experiments.json");
  cereal::JSONOutputArchive archive(handle);
  for (auto &exp : experiments) {
    std::cout << exp << std::endl;
    archive(CEREAL_NVP(exp.p));
  }
  handle.close();
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

void tex() {
  auto experiments = get_experiments();
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
  if (argc == 1) 
    help();
  else if (argc > 1) {
    if (strcmp(argv[1], "experiment") == 0) 
      perform_experiments();
    else if (strcmp(argv[1], "read") == 0 && argc > 2) {
      auto ensemble = read_experiment(argv[2]);
      free(ensemble);
    }
    else if (strcmp(argv[1], "list") == 0) 
      list_experiments();
    else if (strcmp(argv[1], "tex") == 0)
      tex();
    else
      help();
  }

  return 0;
}

vector <data_header> get_experiments() {
  vector<data_header> experiments;

  data_header header_base; 
  header_base.p.Ntrials = 100;
  header_base.p.Np = 64;
  header_base.p.Nt = 1e5;
  header_base.p.a = 0.0;
  header_base.p.k = 1e5;
  header_base.p.l = 1e-7;
  header_base.p.k_BT = 1e-7;
  header_base.p.C = 0.0;
  header_base.p.dt = 1e-7;
  header_base.p.datapoints = 300;
  header_base.p.mu_l = header_base.p.l;
  header_base.p.sigma_l = header_base.p.mu_l / 5;
  header_base.p.mu_theta = 0.0;
  header_base.p.sigma_theta = 0.0;
  header_base.p.space = "log";
  header_base.p.tex_title = "";


  // Standard Rouse model for different dt
  {
    vector<double> list = {1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5}; 
    for (auto dt : list) {
      data_header header = header_base; 
      header.p.dt = dt;

      char title[1024];
      sprintf(title, "Rouse dt %.1e", dt);
      header.p.title = title;
      if (abs(dt - 1e-8) < 1e-10 || abs(dt - 1e-5) < 1e-10)
        header.p.tex_title = "Rouse model";
      experiments.push_back(header);
    }
  }
  
  // vary l
  {
    vector<double> list = {1e-9, 1e-8, 1e-7, 1e-6, 1e-5};
    for (auto l : list) {
      data_header header = header_base; 
      header.p.l = l;

      char title[1024];
      sprintf(title, "Rouse l %.1e", l);
      header.p.title = title;
      experiments.push_back(header);
    }
  }

  // vary k_BT
  {
    vector<double> list = {1e-11, 1e-9, 1e-7, 1e-6, 1e-5, 1e-4};
    for (auto k_BT : list) {
      data_header header = header_base; 
      header.p.k_BT = k_BT;

      char title[1024];
      sprintf(title, "Rouse k_BT %.1e", k_BT);
      header.p.title = title;
      experiments.push_back(header);
    }
  }

  // Vandebroek-Vanderzande
  {
    data_header header = header_base; 
    header.p.dt = 5e-6;
    header.p.C = 1e-6;
    header.p.Ntrials = 20;
    header.p.Nt = 1e6;
    header.p.Np = 16;

    char title[1024];
    sprintf(title, "Vanderzande-Vandebroek C %.1e", header.p.C);
    header.p.title = title;
    experiments.push_back(header);
  }

  // Vandebroek-Vanderzande
  {
    data_header header = header_base; 
    header.p.dt = 1e-7;
    header.p.C = 1e-6;
    header.p.Ntrials = 100;
    header.p.Nt = 1e5;
    header.p.Np = 64;

    char title[1024];
    sprintf(title, "Vanderzande-Vandebroek 2");
    header.p.title = title;
    header.p.tex_title = "zande-broek";
    experiments.push_back(header);
  }

  // Semiflexible, let's choose small k
  {
    vector<double> list = {1e-6, 0.0};

    for (auto C : list) {
      data_header header = header_base; 
      header.p.C = C;
      //header.p.Np = 16;
      //header.p.Nt = 1e6;
      header.p.a = 10;
      header.p.k = 1e-6;
      header.p.Ntrials = 100;
      //header.p.dt = 5e-6;
      if (C == 0) {
        header.p.dt = 1e-8;
        //header.p.Np = 64;
      }
      //header.p.k_BT = 1e-2;

      char title[1024];
      sprintf(title, "Semiflexible C %.1e ", header.p.C);
      header.p.title = title;
      header.p.tex_title = "Semiflexible";
      experiments.push_back(header);
    }
  }

  return experiments;
}
