#ifndef cluster_cuh_
#define cluster_cuh_

#include "event.cuh"
#include "vec4.cuh"

class cluster {
  /**
   * @class cluster
   * @brief A colour-singlet composite particle made from two quarks/antiquarks.
   *
   * A cluster stores the indices i1 and i2 of its two constituents in the event
   * record, together with their combined four-momentum. The momentum is summed
   * from the constituents by the constructor, so callers only supply the two
   * indices (and the event that holds them) — they never build the sum by hand.
   *
   * The stored momentum is a snapshot taken at construction. If a constituent's
   * momentum later changes, rebuild the cluster (add_cluster / set_cluster) to
   * refresh it. Including event.cuh is safe here — event.cuh does not include
   * cluster.cuh, so there is no circular dependency.
   */

 private:
  // ---------------------------------------------------------------------------
  // member variables

  // Indices of the two constituent particles in the event record
  int i1 = -1;
  int i2 = -1;

  // Combined four-momentum (p_i1 + p_i2), set at construction
  vec4 mom;

 public:
  // ---------------------------------------------------------------------------
  // constructor

  __host__ __device__ cluster() {}

  __host__ __device__ cluster(int i1, int i2, const event& ev)
      : i1(i1),
        i2(i2),
        mom(ev.get_particle(i1).get_mom() + ev.get_particle(i2).get_mom()) {}

  // ---------------------------------------------------------------------------
  // getters

  __host__ __device__ int get_i1() const {
    /**
     * @brief get the index of the first constituent particle
     */

    return i1;
  }

  __host__ __device__ int get_i2() const {
    /**
     * @brief get the index of the second constituent particle
     */

    return i2;
  }

  __host__ __device__ vec4 get_mom() const {
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
   * Moving cluster storage out of event keeps cluster data out of the event
   * record during the shower phase where it is not needed.
   *
   * Each hadronisation kernel receives a parallel array of cluster_list objects
   * (one per event) indexed by the same thread ID used to index events.
   *
   * max_particles is defined in base.cuh, which is always in scope before this
   * header is included (following the same convention as event.cuh).
   */

  // Cluster array and occupancy counter
  cluster data[max_particles];
  int n = 0;

  // ---------------------------------------------------------------------------
  // getters

  __host__ __device__ cluster get_cluster(int i) const {
    /**
     * @brief get the cluster at index i
     *
     * @param i the index of the cluster to get
     * @return the cluster at index i
     */

    return data[i];
  }

  __host__ __device__ int get_size() const {
    /**
     * @brief get the number of clusters in the list
     * @return the number of clusters
     */
    return n;
  }

  // ---------------------------------------------------------------------------
  // setters

  __device__ void add_cluster(int i1, int i2, const event& ev) {
    /**
     * @brief add a new cluster to the list
     *
     * The momentum is summed from the constituents by the cluster constructor.
     *
     * @param i1 the index of the first constituent particle
     * @param i2 the index of the second constituent particle
     * @param ev the event holding the constituent particles
     */

    data[n++] = cluster(i1, i2, ev);
  }

  __device__ void set_cluster(int i, int i1, int i2, const event& ev) {
    /**
     * @brief overwrite the cluster at index i
     *
     * The momentum is summed from the constituents by the cluster constructor.
     *
     * @param i  the index of the cluster to overwrite
     * @param i1 the index of the first constituent particle
     * @param i2 the index of the second constituent particle
     * @param ev the event holding the constituent particles
     */

    data[i] = cluster(i1, i2, ev);
  }

  // Clear the list at the start of form_clusters.
  __device__ void reset() {
    /**
     * @brief reset the cluster list
     *
     * This sets the number of clusters to zero, effectively clearing the list.
     */
    n = 0;
  }
};

#endif  // cluster_cuh_
