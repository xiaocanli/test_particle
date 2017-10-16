#ifndef LEN_MAX
#define LEN_MAX 200
#endif

/* Light speed */
#ifndef c0
#define c0 3.0E10
#endif

/* Normalized light speed */
#ifndef c0_norm
#define c0_norm 1.0
#endif

/* Spatial normalization */
#ifndef L0
#define L0 6.96E8
#define iL0 1/L0
#endif

/* sqrt(kT0/m_p), T0 = 1E6 K, CGS unit */
#ifndef VTHE_NORM
#define VTHE_NORM 9.08E6
#endif

#ifndef CHARGE_MASS_PROTON
#define CHARGE_MASS_PROTON 9.58E3
#endif

#ifndef CHARGE_MASS_ELECTRON_NORM
#define CHARGE_MASS_ELECTRON_NORM 1.0
#endif

/* Rest energy of proton in MeV */
#ifndef REST_ENE_PROTON
#define REST_ENE_PROTON 938.272046
#endif

/* Rest energy of electron in MeV */
#ifndef REST_ENE_ELECTRON
#define REST_ENE_ELECTRON 0.511
#endif

typedef enum {false, true} bool;
