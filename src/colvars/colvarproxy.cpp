// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <fstream>
#include <list>
#include <utility>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvar.h"
#include "colvarbias.h"
#include "colvarmodule_utils.h"


colvarproxy::colvarproxy()
{
  colvars = NULL;
  // By default, simulation engines allow to immediately request atoms
  engine_ready_ = true;
  b_simulation_running = true;
  b_simulation_continuing = false;
  b_delete_requested = false;
  version_int = -1;
  features_hash = 0;
  config_queue_ = reinterpret_cast<void *>(new std::list<std::pair<std::string, std::string> >);
}


colvarproxy::~colvarproxy()
{
  close_output_streams();
  if (colvars != NULL) {
    delete colvars;
    colvars = NULL;
  }
  delete reinterpret_cast<std::list<std::pair<std::string, std::string> > *>(config_queue_);
}


bool colvarproxy::io_available()
{
  return true;
}


int colvarproxy::reset()
{
  if (cvm::debug()) {
    cvm::log("colvarproxy::reset()\n");
  }
  int error_code = COLVARS_OK;
  total_force_requested = false;
  return error_code;
}


int colvarproxy::request_deletion()
{
  return cvm::error("Error: \"delete\" command is only available in VMD; "
                    "please use \"reset\" instead.\n",
                    COLVARS_NOT_IMPLEMENTED);
}


void colvarproxy::add_config(std::string const &cmd, std::string const &conf)
{
  reinterpret_cast<std::list<std::pair<std::string, std::string> > *>(config_queue_)->push_back(std::make_pair(cmd, conf));
}


int colvarproxy::setup()
{
  return COLVARS_OK;
}


int colvarproxy::parse_module_config()
{
  int error_code = COLVARS_OK;
  // Read any configuration queued up for Colvars
  std::list<std::pair<std::string, std::string> > *config_queue = reinterpret_cast<std::list<std::pair<std::string, std::string> > *>(config_queue_);
  while (config_queue->size() > 0) {
    std::pair<std::string, std::string> const &p = config_queue->front();
    if (p.first == "config") {
      error_code |= colvars->read_config_string(p.second);
    } else if (p.first == "configfile") {
      error_code |= colvars->read_config_file(p.second.c_str());
    } else {
      error_code |= cvm::error(std::string("Error: invalid keyword \"") +
                               p.first +
                               std::string("\" in colvarproxy::setup()\n"),
                               COLVARS_BUG_ERROR);
    }
    config_queue->pop_front();
  }
  return error_code;
}


int colvarproxy::update_input()
{
  return COLVARS_OK;
}


int colvarproxy::update_output()
{
  return COLVARS_OK;
}


int colvarproxy::end_of_step()
{
  if (cached_alch_lambda_changed) {
    send_alch_lambda();
    cached_alch_lambda_changed = false;
  }
  return COLVARS_OK;
}


int colvarproxy::post_run()
{
  int error_code = COLVARS_OK;
  if (colvars->output_prefix().size()) {
    error_code |= colvars->write_restart_file(cvm::output_prefix()+".colvars.state");
    error_code |= colvars->write_output_files();
  }
  error_code |= flush_output_streams();
  return error_code;
}


void colvarproxy::log(std::string const &message)
{
  fprintf(stdout, "colvars: %s", message.c_str());
}


void colvarproxy::error(std::string const &message)
{
  // TODO handle errors?
  colvarproxy::log(message);
}


void colvarproxy::add_error_msg(std::string const &message)
{
  std::istringstream is(message);
  std::string line;
  while (std::getline(is, line)) {
    error_output += line+"\n";
  }
}


void colvarproxy::clear_error_msgs()
{
  error_output.clear();
}


std::string const & colvarproxy::get_error_msgs()
{
  return error_output;
}


int colvarproxy::get_version_from_string(char const *version_string)
{
  std::string const v(version_string);
  std::istringstream is(v.substr(0, 4) + v.substr(5, 2) + v.substr(8, 2));
  int newint;
  is >> newint;
  return newint;
}


