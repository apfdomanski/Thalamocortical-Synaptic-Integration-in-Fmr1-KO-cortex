#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _NMDA_fdsexp2s_reg();
extern void _ca_reg();
extern void _cad_reg();
extern void _distr_reg();
extern void _fdsexp2s_reg();
extern void _kca_reg();
extern void _km_reg();
extern void _kv_reg();
extern void _na_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," NMDA_fdsexp2s.mod");
fprintf(stderr," ca.mod");
fprintf(stderr," cad.mod");
fprintf(stderr," distr.mod");
fprintf(stderr," fdsexp2s.mod");
fprintf(stderr," kca.mod");
fprintf(stderr," km.mod");
fprintf(stderr," kv.mod");
fprintf(stderr," na.mod");
fprintf(stderr, "\n");
    }
_NMDA_fdsexp2s_reg();
_ca_reg();
_cad_reg();
_distr_reg();
_fdsexp2s_reg();
_kca_reg();
_km_reg();
_kv_reg();
_na_reg();
}
