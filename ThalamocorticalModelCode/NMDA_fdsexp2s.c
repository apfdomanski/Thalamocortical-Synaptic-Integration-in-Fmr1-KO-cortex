/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__NMDA_FDSExp2Syn
#define _nrn_initial _nrn_initial__NMDA_FDSExp2Syn
#define nrn_cur _nrn_cur__NMDA_FDSExp2Syn
#define _nrn_current _nrn_current__NMDA_FDSExp2Syn
#define nrn_jacob _nrn_jacob__NMDA_FDSExp2Syn
#define nrn_state _nrn_state__NMDA_FDSExp2Syn
#define _net_receive _net_receive__NMDA_FDSExp2Syn 
#define _f_rates _f_rates__NMDA_FDSExp2Syn 
#define rates rates__NMDA_FDSExp2Syn 
#define state state__NMDA_FDSExp2Syn 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define e _p[0]
#define tau1 _p[1]
#define tau2 _p[2]
#define mg _p[3]
#define f _p[4]
#define tau_F _p[5]
#define d1 _p[6]
#define tau_D1 _p[7]
#define d2 _p[8]
#define tau_D2 _p[9]
#define i _p[10]
#define g _p[11]
#define A _p[12]
#define B _p[13]
#define Fraction_active _p[14]
#define factor _p[15]
#define DA _p[16]
#define DB _p[17]
#define DFraction_active _p[18]
#define _g _p[19]
#define _tsav _p[20]
#define _nd_area  *_ppvar[0]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_rates();
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "rates", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
#define total total_NMDA_FDSExp2Syn
 double total = 0;
#define usetable usetable_NMDA_FDSExp2Syn
 double usetable = 1;
#define vmax vmax_NMDA_FDSExp2Syn
 double vmax = 100;
#define vmin vmin_NMDA_FDSExp2Syn
 double vmin = -120;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "d2", 0, 1,
 "d1", 0, 1,
 "f", 0, 1e+009,
 "tau_D2", 1e-009, 1e+009,
 "tau_D1", 1e-009, 1e+009,
 "tau_F", 1e-009, 1e+009,
 "tau2", 1e-009, 1e+009,
 "tau1", 1e-009, 1e+009,
 "usetable_NMDA_FDSExp2Syn", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "vmin_NMDA_FDSExp2Syn", "mV",
 "vmax_NMDA_FDSExp2Syn", "mV",
 "total_NMDA_FDSExp2Syn", "umho",
 "e", "mV",
 "tau1", "ms",
 "tau2", "ms",
 "mg", "mM",
 "f", "1",
 "tau_F", "ms",
 "d1", "1",
 "tau_D1", "ms",
 "d2", "1",
 "tau_D2", "ms",
 "A", "umho",
 "B", "umho",
 "i", "nA",
 "g", "umho",
 0,0
};
 static double A0 = 0;
 static double B0 = 0;
 static double Fraction_active0 = 0;
 static double delta_t = 0.01;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "vmin_NMDA_FDSExp2Syn", &vmin_NMDA_FDSExp2Syn,
 "vmax_NMDA_FDSExp2Syn", &vmax_NMDA_FDSExp2Syn,
 "total_NMDA_FDSExp2Syn", &total_NMDA_FDSExp2Syn,
 "usetable_NMDA_FDSExp2Syn", &usetable_NMDA_FDSExp2Syn,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"NMDA_FDSExp2Syn",
 "e",
 "tau1",
 "tau2",
 "mg",
 "f",
 "tau_F",
 "d1",
 "tau_D1",
 "d2",
 "tau_D2",
 0,
 "i",
 "g",
 0,
 "A",
 "B",
 "Fraction_active",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 21, _prop);
 	/*initialize range parameters*/
 	e = 0;
 	tau1 = 2.23;
 	tau2 = 150;
 	mg = 1.2;
 	f = 0.917;
 	tau_F = 94;
 	d1 = 0.416;
 	tau_D1 = 380;
 	d2 = 0.975;
 	tau_D2 = 9200;
  }
 	_prop->param = _p;
 	_prop->param_size = 21;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _net_receive(Point_process*, double*, double);
 static void _net_init(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _NMDA_fdsexp2s_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 21, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 5;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 NMDA_FDSExp2Syn C:/Users/ad15419/Dropbox/aleks paper/2018/Box Docs/NatCommsSubmission/FFImodelCode/NMDA Depressing synapses with spiking big run/NMDA_fdsexp2s.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_Fraction_active;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_rates(double);
static int rates(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_rates(double);
 static int _slist1[2], _dlist1[2];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   DA = - A / tau1 ;
   DB = - B / tau2 ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 DA = DA  / (1. - dt*( ( - 1.0 ) / tau1 )) ;
 DB = DB  / (1. - dt*( ( - 1.0 ) / tau2 )) ;
  return 0;
}
 /*END CVODE*/
 static int state () {_reset=0;
 {
    A = A + (1. - exp(dt*(( - 1.0 ) / tau1)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - A) ;
    B = B + (1. - exp(dt*(( - 1.0 ) / tau2)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - B) ;
   }
  return 0;
}
 static double _mfac_rates, _tmin_rates;
 static void _check_rates();
 static void _check_rates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_mg;
  if (!usetable) {return;}
  if (_sav_mg != mg) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rates =  vmin ;
   _tmax =  vmax ;
   _dx = (_tmax - _tmin_rates)/200.; _mfac_rates = 1./_dx;
   for (_i=0, _x=_tmin_rates; _i < 201; _x += _dx, _i++) {
    _f_rates(_x);
    _t_Fraction_active[_i] = Fraction_active;
   }
   _sav_mg = mg;
  }
 }

 static int rates(double _lv){ _check_rates();
 _n_rates(_lv);
 return 0;
 }

 static void _n_rates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rates(_lv); return; 
}
 _xi = _mfac_rates * (_lv - _tmin_rates);
 if (isnan(_xi)) {
  Fraction_active = _xi;
  return;
 }
 if (_xi <= 0.) {
 Fraction_active = _t_Fraction_active[0];
 return; }
 if (_xi >= 200.) {
 Fraction_active = _t_Fraction_active[200];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 Fraction_active = _t_Fraction_active[_i] + _theta*(_t_Fraction_active[_i+1] - _t_Fraction_active[_i]);
 }

 
static int  _f_rates (  double _lv ) {
   Fraction_active = 1.0 / ( 1.0 + exp ( 0.062 * - _lv ) * ( mg / 3.57 ) ) ;
    return 0; }
 
static double _hoc_rates(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
  _r = 1.;
 rates (  *getarg(1) );
 return(_r);
}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{    _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   _args[1] = 1.0 + ( _args[1] - 1.0 ) * exp ( - ( t - _args[4] ) / tau_F ) ;
   _args[2] = 1.0 - ( 1.0 - _args[2] ) * exp ( - ( t - _args[4] ) / tau_D1 ) ;
   _args[3] = 1.0 - ( 1.0 - _args[3] ) * exp ( - ( t - _args[4] ) / tau_D2 ) ;
   _args[4] = t ;
       if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = A;
    double __primary = (A + _args[0] * factor * _args[1] * _args[2] * _args[3] ) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau1 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - __primary );
    A += __primary;
  } else {
 A = A + _args[0] * factor * _args[1] * _args[2] * _args[3]  ;
     }
     if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = B;
    double __primary = (B + _args[0] * factor * _args[1] * _args[2] * _args[3] ) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau2 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - __primary );
    B += __primary;
  } else {
 B = B + _args[0] * factor * _args[1] * _args[2] * _args[3]  ;
     }
 total = total + _args[0] * _args[1] * _args[2] * _args[3] ;
   _args[1] = _args[1] + f ;
   _args[2] = _args[2] * d1 ;
   _args[3] = _args[3] * d2 ;
   } }
 
static void _net_init(Point_process* _pnt, double* _args, double _lflag) {
    _args[1] = 1.0 ;
   _args[2] = 1.0 ;
   _args[3] = 1.0 ;
   _args[4] = t ;
   }
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 ();
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  A = A0;
  B = B0;
  Fraction_active = Fraction_active0;
 {
   double _ltp ;
 rates ( _threadargscomma_ v ) ;
   total = 0.0 ;
   if ( tau1 / tau2 > 0.9999 ) {
     tau1 = 0.9999 * tau2 ;
     }
   A = 0.0 ;
   B = 0.0 ;
   _ltp = ( tau1 * tau2 ) / ( tau2 - tau1 ) * log ( tau2 / tau1 ) ;
   factor = - exp ( - _ltp / tau1 ) + exp ( - _ltp / tau2 ) ;
   factor = 1.0 / factor ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   rates ( _threadargscomma_ v ) ;
   g = ( B - A ) ;
   i = g * Fraction_active * ( v - e ) ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 { error =  state();
 if(error){fprintf(stderr,"at line 94 in file NMDA_fdsexp2s.mod:\n	SOLVE state METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 }}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(A) - _p;  _dlist1[0] = &(DA) - _p;
 _slist1[1] = &(B) - _p;  _dlist1[1] = &(DB) - _p;
   _t_Fraction_active = makevector(201*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "NMDA_fdsexp2s.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "NMDA channel\n"
  "This is an adapted version of Exp2Syn.\n"
  "Adapted by Kevin M Biddell similar to as described by wolf et al 2006\n"
  "4/21/07\n"
  "verified 3/29/2012 KMB \n"
  "kevin.biddell@gmail.com\n"
  "\n"
  "Voltage dependence of Mg2+ block:\n"
  "	Jahr & Stevens 1990. J Neurosci 10: 1830.\n"
  "	Jahr & Stevens 1990. J Neurosci 10: 3178.\n"
  "\n"
  "Two state kinetic scheme synapse described by rise time tau1,\n"
  "and decay time constant tau2. The normalized peak condunductance is 1.\n"
  "Decay time MUST be greater than rise time.\n"
  "\n"
  "The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is\n"
  " A = a*exp(-t/tau1) and\n"
  " G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))\n"
  "	where tau1 < tau2\n"
  "\n"
  "If tau2-tau1 -> 0 then we have a alphasynapse.\n"
  "and if tau1 -> 0 then we have just single exponential decay.\n"
  "\n"
  "The factor is evaluated in the\n"
  "initial block such that an event of weight 1 generates a\n"
  "peak conductance of 1.\n"
  "\n"
  "Because the solution is a sum of exponentials, the\n"
  "coupled equations can be solved as a pair of independent equations\n"
  "by the more efficient cnexp method.\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS NMDA_FDSExp2Syn\n"
  "	RANGE tau1, tau2, g, e, i, mg, f, tau_F, d1, tau_D1, d2, tau_D2\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	GLOBAL total,vmin, vmax\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(umho) = (micromho)\n"
  "	(mM) = (milli/liter)\n"
  "}	\n"
  "\n"
  "PARAMETER {\n"
  "	e	= 0    (mV)	: reversal potential\n"
  "	tau1	= 2.23  (ms)<1e-9,1e9>	: changed by KMB 5/07/07 from 2.82 to 20.77% less (wolf et 2006)\n"
  "	tau2	= 150   (ms)<1e-9,1e9> 	: changed by KMB 4/23/07 from 160 to 52.7% less (wolf et al 2006)\n"
  "	mg	= 1.2    (mM)	: external magnesium concentration\n"
  "	vmin = -120	(mV)\n"
  "	vmax = 100	(mV)\n"
  "	f = 0.917 (1) < 0, 1e9 >    : facilitation\n"
  "        tau_F = 94 (ms) < 1e-9, 1e9 >\n"
  "        d1 = 0.416 (1) < 0, 1 >     : fast depression\n"
  "        tau_D1 = 380 (ms) < 1e-9, 1e9 >\n"
  "        d2 = 0.975 (1) < 0, 1 >     : slow depression\n"
  "        tau_D2 = 9200 (ms) < 1e-9, 1e9 >\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	i (nA)\n"
  "	g (umho)\n"
  "	factor\n"
  "	total (umho)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	A (umho)\n"
  "	B (umho)\n"
  "	Fraction_active		: fraction free of Mg2+ block\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	LOCAL tp\n"
  "	rates(v)\n"
  "	total = 0\n"
  "	if (tau1/tau2 > 0.9999) {\n"
  "		tau1 = 0.9999*tau2\n"
  "	}\n"
  "	A = 0\n"
  "	B = 0\n"
  "	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)\n"
  "	factor = -exp(-tp/tau1) + exp(-tp/tau2)\n"
  "	factor = 1/factor\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	rates(v)\n"
  "	SOLVE state METHOD cnexp\n"
  "        g = (B-A) 	\n"
  "        i = g*Fraction_active*(v - e) : N.B. C deals with fraction free of Mg++ block	\n"
  "}\n"
  "\n"
  "\n"
  "DERIVATIVE state {\n"
  "	A' = -A/tau1\n"
  "	B' = -B/tau2\n"
  "}\n"
  "	\n"
  "PROCEDURE rates(v(mV)) {\n"
  "	TABLE Fraction_active\n"
  "	DEPEND mg\n"
  "	FROM vmin TO vmax WITH 200\n"
  "\n"
  "	: from Jahr & Stevens\n"
  "\n"
  "	Fraction_active = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))\n"
  "}\n"
  "\n"
  "\n"
  "NET_RECEIVE(weight (umho), F, D1, D2, tsyn (ms)) {\n"
  "INITIAL {\n"
  ": these are in NET_RECEIVE to be per-stream\n"
  "        F = 1\n"
  "        D1 = 1\n"
  "        D2 = 1\n"
  "        tsyn = t\n"
  ": this header will appear once per stream\n"
  ": printf(\"t\\t t-tsyn\\t F\\t D1\\t D2\\t amp\\t newF\\t newD1\\t newD2\\n\")\n"
  "}\n"
  "\n"
  "        F = 1 + (F-1)*exp(-(t - tsyn)/tau_F)\n"
  "        D1 = 1 - (1-D1)*exp(-(t - tsyn)/tau_D1)\n"
  "        D2 = 1 - (1-D2)*exp(-(t - tsyn)/tau_D2)\n"
  ": printf(\"%g\\t%g\\t%g\\t%g\\t%g\\t%g\", t, t-tsyn, F, D1, D2, weight*F*D1*D2)\n"
  "        tsyn = t\n"
  "\n"
  "	state_discontinuity(A, A + weight*factor*F*D1*D2)\n"
  "	state_discontinuity(B, B + weight*factor*F*D1*D2)\n"
  "	total = total+weight*F*D1*D2\n"
  "\n"
  "        F = F + f\n"
  "        D1 = D1 * d1\n"
  "        D2 = D2 * d2\n"
  ": printf(\"\\t%g\\t%g\\t%g\\n\", F, D1, D2)\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  ;
#endif
