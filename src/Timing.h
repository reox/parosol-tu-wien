/*
 * ParOSol: a parallel FE solver for trabecular bone modeling
 * Copyright (C) 2011, Cyril Flaig
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TIMING_H
#define TIMING_H

#include <map>
#include <mpi.h>
#include <string>
#include <sys/time.h>

#define COUTTIME(msg)                                                          \
  "avg: " << msg.avg << " min: " << msg.min << " max: " << msg.max

struct t_timing {
  double min;
  double avg;
  double max;
};

/*! A simple class to time some parts of the code. This object can start and
stop different timers.
 *
 * \verbatim
t.start("a"); bla();t.stop("a"); foo(); t.Restart("a"); foobar(); t.Stop("a");
t_timing t_elapsed = t.ElapsedTime("a");
\endverbatim
 */

class Timer {
public:
  typedef struct timeval t_time;
  /**
   * @brief Constructor
   *
   * @param a_comm MPI Communicator
   */

  Timer(MPI_Comm a_comm) : _comm(a_comm) { MPI_Comm_size(a_comm, &Size); }

  /**
   * @brief Starts a named timer. It uses a Barrier to start the timer on all
   * cpu at the same time.
   *
   * @param timer name of the timer
   */

  void Start(std::string timer) {
    t_time now;
    _elapsed[timer] = 0;
    MPI_Barrier(_comm);
    gettimeofday(&now, NULL);
    _timers[timer] = now;
  }

  /**
   * @brief Restarts a stopped timer.
   *
   * Restart should only be used after the timer is stopped.
   *
   * @param timer name of the timer
   */
  void Restart(std::string timer) {
    t_time now;
    if (_elapsed.find(timer) == _elapsed.end()) {
      _elapsed[timer] = 0;
    }
    gettimeofday(&now, NULL);
    _timers[timer] = now;
  }

  /**
   * @brief Stops a running timer.
   *
   * @param timer name of the timer
   */
  void Stop(std::string timer) {
    t_time now;
    gettimeofday(&now, NULL);
    t_time start = _timers[timer];
    _timers.erase(timer);
    double res = difftime(&start, &now);
    _elapsed[timer] = _elapsed[timer] + res;
    return;
  }

  /**
   * @brief Computes the time that is elapsed between staring and stopping the
   * timer.
   *
   * It uses MPI_Allreduce to compute the timings
   *
   * @param timer name of the timer
   * @return t_timing struct with the timings.
   */
  t_timing ElapsedTime(std::string timer) {
    double res = _elapsed[timer];
    double min = 0, max = 0, avg = 0;
    MPI_Allreduce(&res, &max, 1, MPI_DOUBLE, MPI_MAX, _comm);
    MPI_Allreduce(&res, &min, 1, MPI_DOUBLE, MPI_MIN, _comm);
    MPI_Allreduce(&res, &avg, 1, MPI_DOUBLE, MPI_SUM, _comm);
    avg = avg / Size;
    t_timing t;
    t.min = min;
    t.max = max;
    t.avg = avg;
    return t;
  }

private:
  std::map<std::string, t_time> _timers;
  std::map<std::string, double> _elapsed;
  MPI_Comm _comm;
  int Size;

  /// Compute the time between two timing objects.
  double difftime(t_time *start, t_time *end) {
    double res = 0;
    res = end->tv_sec - start->tv_sec;
    res += (((double)end->tv_usec) - start->tv_usec) / 1e6;
    return res;
  }
};
#endif
