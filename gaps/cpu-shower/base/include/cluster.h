#ifndef cluster_h_
#define cluster_h_

#include "vec4.h"

class cluster {
  /**
   * @class cluster
   * @brief A colour-singlet composite particle made from two quarks/antiquarks.
   *
   * A cluster stores the indices i1 and i2 of its two constituents in the
   * event record, rather than pointers or copies. This keeps the design safe
   * for both CPU and GPU — no dynamic allocation or pointer arithmetic.
   *
   * The combined momentum is passed in at construction; the event is not
   * referenced here so that event.h can include cluster.h without a circular
   * dependency.
   *
   * pid = 81 is the PDG code for a generic cluster.
   */

 private:
  // ---------------------------------------------------------------------------
  // member variables

  // Indices of the two constituent particles in the event record
  int i1 = -1;
  int i2 = -1;

  // Combined four-momentum (p1 + p2)
  vec4 mom;

  // Cluster pid (81 = generic cluster)
  int pid = 81;

 public:
  // ---------------------------------------------------------------------------
  // constructor

  cluster() = default;

  cluster(int i1, int i2, vec4 mom) : i1(i1), i2(i2), mom(mom), pid(81) {}

  // ---------------------------------------------------------------------------
  // getters

  int get_i1() const {
    /**
     * @brief get the index of the first constituent particle
     */

    return i1;
  }

  int get_i2() const {
    /**
     * @brief get the index of the second constituent particle
     */

    return i2;
  }

  int get_pid() const {
    /**
     * @brief get the PDG code of the cluster
     */

    return pid;
  }

  vec4 get_mom() const {
    /**
     * @brief get the combined four-momentum of the cluster
     */

    return mom;
  }
};

struct cluster_list {
  /**
   * @struct cluster_list
   * @brief A fixed-size array of clusters for a single event, held separately
   * from the event record.
   *
   * Moving cluster storage out of event reduces the event struct size by
   * ~3.4 KB (max_particles * sizeof(cluster)), keeping cluster data out of
   * the event record during the shower phase where it is not needed.
   *
   * A single cluster_list is constructed in hadronisation::run() and passed by
   * reference to form_clusters, fission_clusters, and decay_clusters. It goes
   * out of scope when run() returns.
   *
   * max_particles is defined in base.h, which is always in scope before this
   * header is included (following the same convention as event.h).
   */

  // Cluster array and occupancy counter
  cluster data[max_particles];
  int n = 0;

  // ---------------------------------------------------------------------------
  // getters

  cluster get_cluster(int i) const {
    /**
     * @brief get the cluster at index i
     *
     * @param i the index of the cluster to get
     * @return the cluster at index i
     */

    return data[i];
  }

  int get_size() const {
    /**
     * @brief get the number of clusters in the list
     * @return the number of clusters
     */
    return n;
  }

  // ---------------------------------------------------------------------------
  // setters

  void add_cluster(int i1, int i2, vec4 mom) {
    /**
     * @brief add a new cluster to the list
     *
     * @param i1 the index of the first constituent particle
     * @param i2 the index of the second constituent particle
     * @param mom the combined four-momentum of the cluster
     */

    data[n++] = cluster(i1, i2, mom);
  }

  void set_cluster(int i, cluster c) {
    /**
     * @brief overwrite the cluster at index i
     *
     * @param i the index of the cluster to overwrite
     * @param c the new cluster
     */

    data[i] = c;
  }

  // Clear the list at the start of form_clusters.
  void reset() {
    /**
     * @brief reset the cluster list
     *
     * This sets the number of clusters to zero, effectively clearing the list.
     */
    n = 0;
  }
};

#endif  // cluster_h_
