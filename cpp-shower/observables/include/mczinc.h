#ifndef mczinc_h_
#define mczinc_h_

#include "event.h"

// breit frame event shapes - simple calculations so only one
void calculate_zboson_obs(event& ev);
void calculate_lepton_obs(event& ev);

#endif  // breitframe_h_