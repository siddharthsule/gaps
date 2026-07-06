#include "hadronisation.cuh"

__global__ void form_clusters(event* events, cluster_list* cls,
                              hadronisation* had, int n) {
  /**
   * @brief form clusters from the particles in the event
   *
   * A cluster is a composite particle made from two particles. It is a colour
   * singlet, formed by a quark and an antiquark. We loop over all pairs of
   * particles in the event, and if they are colour connected, we form a
   * cluster from them. We set the first particle to the cluster, and remove
   * the second particle from the event.
   *
   * @param ev the event to form clusters from
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Event Preamble
  event& ev = events[idx];
  cluster_list& cl = cls[idx];
  cl.reset();
  // ---------------------------------------------

  // Loop over all particle pairs in the event
  for (int i = 2; i < ev.get_size(); i++) {
    for (int j = i + 1; j < ev.get_size(); j++) {
      // Skip non-partons
      if (!ev.get_particle(i).is_parton() || !ev.get_particle(j).is_parton())
        continue;

      // Do nothing if the partons are not color connected
      if (!ev.get_particle(i).is_color_connected(ev.get_particle(j))) continue;

      // Create a cluster from the two particles (q first, then qbar)
      if (ev.get_particle(i).get_pid() > 0) {
        cl.add_cluster(
            i, j, ev.get_particle(i).get_mom() + ev.get_particle(j).get_mom());
      } else {
        cl.add_cluster(
            j, i, ev.get_particle(j).get_mom() + ev.get_particle(i).get_mom());
      }

      // Skip to the next iteration of the inner loop
      break;
    }
  }
}
