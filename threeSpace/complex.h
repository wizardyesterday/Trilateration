matutil.h                                                                                           0100600 0000764 0000144 00000002200 07676241537 011745  0                                                                                                    ustar   chris                           users                                                                                                                                                                                                                  /*************************************************************************
 file name: matutil.h
*************************************************************************/

#ifndef __MATUTIL_H__
#define __MATUTIL_H__

/*************************************************************************
 defines    
*************************************************************************/
#define MAT_MAX_SIZE (50 + 1)

/*************************************************************************
 structure definitions     
*************************************************************************/     

/*************************************************************************
 function prototypes     
*************************************************************************/
void matInit(void);
void matCopy(int m,int n,float **a,float **b);
void matTranspose(int m,int n,float **a,float **aT);
void vecDotProduct(int n,float *v1,float *v2,float *v1DotV2);
void matAdd(int m,int n,float **a,float **b,float **c);
void matMultiply(int m,int n,int p,float **a,float **b,float **c);
void matInverse(int n,float **a,float **aInv);

#endif /* __MATUTIL_H__ */
                                                                                                                                                                                                                                                                                                                                                                                                nr.h                                                                                                0100600 0000764 0000144 00000075357 07675067760 010737  0                                                                                                    ustar   chris                           users                                                                                                                                                                                                                  #ifndef _NR_H_
#define _NR_H_

#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX {float r,i;} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

#ifndef _ARITHCODE_DECLARE_T_
typedef struct {
	unsigned long *ilob,*iupb,*ncumfq,jdif,nc,minint,nch,ncum,nrad;
} arithcode;
#define _ARITHCODE_DECLARE_T_
#endif /* _ARITHCODE_DECLARE_T_ */

#ifndef _HUFFCODE_DECLARE_T_
typedef struct {
	unsigned long *icod,*ncod,*left,*right,nch,nodemax;
} huffcode;
#define _HUFFCODE_DECLARE_T_
#endif /* _HUFFCODE_DECLARE_T_ */

#include <stdio.h>

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

void addint(double **uf, double **uc, double **res, int nf);
void airy(float x, float *ai, float *bi, float *aip, float *bip);
void amebsa(float **p, float y[], int ndim, float pb[],	float *yb,
	float ftol, float (*funk)(float []), int *iter, float temptr);
void amoeba(float **p, float y[], int ndim, float ftol,
	float (*funk)(float []), int *iter);
float amotry(float **p, float y[], float psum[], int ndim,
	float (*funk)(float []), int ihi, float fac);
float amotsa(float **p, float y[], float psum[], int ndim, float pb[],
	float *yb, float (*funk)(float []), int ihi, float *yhi, float fac);
void anneal(float x[], float y[], int iorder[], int ncity);
double anorm2(double **a, int n);
void arcmak(unsigned long nfreq[], unsigned long nchh, unsigned long nradd,
	arithcode *acode);
void arcode(unsigned long *ich, unsigned char **codep, unsigned long *lcode,
	unsigned long *lcd, int isign, arithcode *acode);
void arcsum(unsigned long iin[], unsigned long iout[], unsigned long ja,
	int nwk, unsigned long nrad, unsigned long nc);
void asolve(unsigned long n, double b[], double x[], int itrnsp);
void atimes(unsigned long n, double x[], double r[], int itrnsp);
void avevar(float data[], unsigned long n, float *ave, float *var);
void balanc(float **a, int n);
void banbks(float **a, unsigned long n, int m1, int m2, float **al,
	unsigned long indx[], float b[]);
void bandec(float **a, unsigned long n, int m1, int m2, float **al,
	unsigned long indx[], float *d);
void banmul(float **a, unsigned long n, int m1, int m2, float x[], float b[]);
void bcucof(float y[], float y1[], float y2[], float y12[], float d1,
	float d2, float **c);
void bcuint(float y[], float y1[], float y2[], float y12[],
	float x1l, float x1u, float x2l, float x2u, float x1,
	float x2, float *ansy, float *ansy1, float *ansy2);
void beschb(double x, double *gam1, double *gam2, double *gampl,
	double *gammi);
float bessi(int n, float x);
float bessi0(float x);
float bessi1(float x);
void bessik(float x, float xnu, float *ri, float *rk, float *rip,
	float *rkp);
float bessj(int n, float x);
float bessj0(float x);
float bessj1(float x);
void bessjy(float x, float xnu, float *rj, float *ry, float *rjp,
	float *ryp);
float bessk(int n, float x);
float bessk0(float x);
float bessk1(float x);
float bessy(int n, float x);
float bessy0(float x);
float bessy1(float x);
float beta(float z, float w);
float betacf(float a, float b, float x);
float betai(float a, float b, float x);
float bico(int n, int k);
void bksub(int ne, int nb, int jf, int k1, int k2, float ***c);
float bnldev(float pp, int n, long *idum);
float brent(float ax, float bx, float cx,
	float (*f)(float), float tol, float *xmin);
void broydn(float x[], int n, int *check,
	void (*vecfunc)(int, float [], float []));
void bsstep(float y[], float dydx[], int nv, float *xx, float htry,
	float eps, float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []));
void caldat(long julian, int *mm, int *id, int *iyyy);
void chder(float a, float b, float c[], float cder[], int n);
float chebev(float a, float b, float c[], int m, float x);
void chebft(float a, float b, float c[], int n, float (*func)(float));
void chebpc(float c[], float d[], int n);
void chint(float a, float b, float c[], float cint[], int n);
float chixy(float bang);
void choldc(float **a, int n, float p[]);
void cholsl(float **a, int n, float p[], float b[], float x[]);
void chsone(float bins[], float ebins[], int nbins, int knstrn,
	float *df, float *chsq, float *prob);
void chstwo(float bins1[], float bins2[], int nbins, int knstrn,
	float *df, float *chsq, float *prob);
void cisi(float x, float *ci, float *si);
void cntab1(int **nn, int ni, int nj, float *chisq,
	float *df, float *prob, float *cramrv, float *ccc);
void cntab2(int **nn, int ni, int nj, float *h, float *hx, float *hy,
	float *hygx, float *hxgy, float *uygx, float *uxgy, float *uxy);
void convlv(float data[], unsigned long n, float respns[], unsigned long m,
	int isign, float ans[]);
void copy(double **aout, double **ain, int n);
void correl(float data1[], float data2[], unsigned long n, float ans[]);
void cosft(float y[], int n, int isign);
void cosft1(float y[], int n);
void cosft2(float y[], int n, int isign);
void covsrt(float **covar, int ma, int ia[], int mfit);
void crank(unsigned long n, float w[], float *s);
void cyclic(float a[], float b[], float c[], float alpha, float beta,
	float r[], float x[], unsigned long n);
void daub4(float a[], unsigned long n, int isign);
float dawson(float x);
float dbrent(float ax, float bx, float cx,
	float (*f)(float), float (*df)(float), float tol, float *xmin);
void ddpoly(float c[], int nc, float x, float pd[], int nd);
int decchk(char string[], int n, char *ch);
void derivs(float x, float y[], float dydx[]);
float df1dim(float x);
void dfour1(double data[], unsigned long nn, int isign);
void dfpmin(float p[], int n, float gtol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []));
float dfridr(float (*func)(float), float x, float h, float *err);
void dftcor(float w, float delta, float a, float b, float endpts[],
	float *corre, float *corim, float *corfac);
void dftint(float (*func)(float), float a, float b, float w,
	float *cosint, float *sinint);
void difeq(int k, int k1, int k2, int jsf, int is1, int isf,
	int indexv[], int ne, float **s, float **y);
void dlinmin(float p[], float xi[], int n, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float[]));
double dpythag(double a, double b);
void drealft(double data[], unsigned long n, int isign);
void dsprsax(double sa[], unsigned long ija[], double x[], double b[],
	unsigned long n);
void dsprstx(double sa[], unsigned long ija[], double x[], double b[],
	unsigned long n);
void dsvbksb(double **u, double w[], double **v, int m, int n, double b[],
	double x[]);
void dsvdcmp(double **a, int m, int n, double w[], double **v);
void eclass(int nf[], int n, int lista[], int listb[], int m);
void eclazz(int nf[], int n, int (*equiv)(int, int));
float ei(float x);
void eigsrt(float d[], float **v, int n);
float elle(float phi, float ak);
float ellf(float phi, float ak);
float ellpi(float phi, float en, float ak);
void elmhes(float **a, int n);
float erfcc(float x);
float erff(float x);
float erffc(float x);
void eulsum(float *sum, float term, int jterm, float wksp[]);
float evlmem(float fdt, float d[], int m, float xms);
float expdev(long *idum);
float expint(int n, float x);
float f1(float x);
float f1dim(float x);
float f2(float y);
float f3(float z);
float factln(int n);
float factrl(int n);
void fasper(float x[], float y[], unsigned long n, float ofac, float hifac,
	float wk1[], float wk2[], unsigned long nwk, unsigned long *nout,
	unsigned long *jmax, float *prob);
void fdjac(int n, float x[], float fvec[], float **df,
	void (*vecfunc)(int, float [], float []));
void fgauss(float x, float a[], float *y, float dyda[], int na);
void fill0(double **u, int n);
void fit(float x[], float y[], int ndata, float sig[], int mwt,
	float *a, float *b, float *siga, float *sigb, float *chi2, float *q);
void fitexy(float x[], float y[], int ndat, float sigx[], float sigy[],
	float *a, float *b, float *siga, float *sigb, float *chi2, float *q);
void fixrts(float d[], int m);
void fleg(float x, float pl[], int nl);
void flmoon(int n, int nph, long *jd, float *frac);
float fmin(float x[]);
void four1(float data[], unsigned long nn, int isign);
void fourew(FILE *file[5], int *na, int *nb, int *nc, int *nd);
void fourfs(FILE *file[5], unsigned long nn[], int ndim, int isign);
void fourn(float data[], unsigned long nn[], int ndim, int isign);
void fpoly(float x, float p[], int np);
void fred2(int n, float a, float b, float t[], float f[], float w[],
	float (*g)(float), float (*ak)(float, float));
float fredin(float x, int n, float a, float b, float t[], float f[], float w[],
	float (*g)(float), float (*ak)(float, float));
void frenel(float x, float *s, float *c);
void frprmn(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []));
void ftest(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *f, float *prob);
float gamdev(int ia, long *idum);
float gammln(float xx);
float gammp(float a, float x);
float gammq(float a, float x);
float gasdev(long *idum);
void gaucof(int n, float a[], float b[], float amu0, float x[], float w[]);
void gauher(float x[], float w[], int n);
void gaujac(float x[], float w[], int n, float alf, float bet);
void gaulag(float x[], float w[], int n, float alf);
void gauleg(float x1, float x2, float x[], float w[], int n);
void gaussj(float **a, int n, float **b, int m);
void gcf(float *gammcf, float a, float x, float *gln);
float golden(float ax, float bx, float cx, float (*f)(float), float tol,
	float *xmin);
void gser(float *gamser, float a, float x, float *gln);
void hpsel(unsigned long m, unsigned long n, float arr[], float heap[]);
void hpsort(unsigned long n, float ra[]);
void hqr(float **a, int n, float wr[], float wi[]);
void hufapp(unsigned long index[], unsigned long nprob[], unsigned long n,
	unsigned long i);
void hufdec(unsigned long *ich, unsigned char *code, unsigned long lcode,
	unsigned long *nb, huffcode *hcode);
void hufenc(unsigned long ich, unsigned char **codep, unsigned long *lcode,
	unsigned long *nb, huffcode *hcode);
void hufmak(unsigned long nfreq[], unsigned long nchin, unsigned long *ilong,
	unsigned long *nlong, huffcode *hcode);
void hunt(float xx[], unsigned long n, float x, unsigned long *jlo);
void hypdrv(float s, float yy[], float dyyds[]);
fcomplex hypgeo(fcomplex a, fcomplex b, fcomplex c, fcomplex z);
void hypser(fcomplex a, fcomplex b, fcomplex c, fcomplex z,
	fcomplex *series, fcomplex *deriv);
unsigned short icrc(unsigned short crc, unsigned char *bufptr,
	unsigned long len, short jinit, int jrev);
unsigned short icrc1(unsigned short crc, unsigned char onech);
unsigned long igray(unsigned long n, int is);
void iindexx(unsigned long n, long arr[], unsigned long indx[]);
void indexx(unsigned long n, float arr[], unsigned long indx[]);
void interp(double **uf, double **uc, int nf);
int irbit1(unsigned long *iseed);
int irbit2(unsigned long *iseed);
void jacobi(float **a, int n, float d[], float **v, int *nrot);
void jacobn(float x, float y[], float dfdx[], float **dfdy, int n);
long julday(int mm, int id, int iyyy);
void kendl1(float data1[], float data2[], unsigned long n, float *tau, float *z,
	float *prob);
void kendl2(float **tab, int i, int j, float *tau, float *z, float *prob);
void kermom(double w[], double y, int m);
void ks2d1s(float x1[], float y1[], unsigned long n1,
	void (*quadvl)(float, float, float *, float *, float *, float *),
	float *d1, float *prob);
void ks2d2s(float x1[], float y1[], unsigned long n1, float x2[], float y2[],
	unsigned long n2, float *d, float *prob);
void ksone(float data[], unsigned long n, float (*func)(float), float *d,
	float *prob);
void kstwo(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *d, float *prob);
void laguer(fcomplex a[], int m, fcomplex *x, int *its);
void lfit(float x[], float y[], float sig[], int ndat, float a[], int ia[],
	int ma, float **covar, float *chisq, void (*funcs)(float, float [], int));
void linbcg(unsigned long n, double b[], double x[], int itol, double tol,
	 int itmax, int *iter, double *err);
void linmin(float p[], float xi[], int n, float *fret,
	float (*func)(float []));
void lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
	 float *f, float stpmax, int *check, float (*func)(float []));
void load(float x1, float v[], float y[]);
void load1(float x1, float v1[], float y[]);
void load2(float x2, float v2[], float y[]);
void locate(float xx[], unsigned long n, float x, unsigned long *j);
void lop(double **out, double **u, int n);
void lubksb(float **a, int n, int *indx, float b[]);
void ludcmp(float **a, int n, int *indx, float *d);
void machar(int *ibeta, int *it, int *irnd, int *ngrd,
	int *machep, int *negep, int *iexp, int *minexp, int *maxexp,
	float *eps, float *epsneg, float *xmin, float *xmax);
void matadd(double **a, double **b, double **c, int n);
void matsub(double **a, double **b, double **c, int n);
void medfit(float x[], float y[], int ndata, float *a, float *b, float *abdev);
void memcof(float data[], int n, int m, float *xms, float d[]);
int metrop(float de, float t);
void mgfas(double **u, int n, int maxcyc);
void mglin(double **u, int n, int ncycle);
float midexp(float (*funk)(float), float aa, float bb, int n);
float midinf(float (*funk)(float), float aa, float bb, int n);
float midpnt(float (*func)(float), float a, float b, int n);
float midsql(float (*funk)(float), float aa, float bb, int n);
float midsqu(float (*funk)(float), float aa, float bb, int n);
void miser(float (*func)(float []), float regn[], int ndim, unsigned long npts,
	float dith, float *ave, float *var);
void mmid(float y[], float dydx[], int nvar, float xs, float htot,
	int nstep, float yout[], void (*derivs)(float, float[], float[]));
void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,
	float *fc, float (*func)(float));
void mnewt(int ntrial, float x[], int n, float tolx, float tolf);
void moment(float data[], int n, float *ave, float *adev, float *sdev,
	float *var, float *skew, float *curt);
void mp2dfr(unsigned char a[], unsigned char s[], int n, int *m);
void mpadd(unsigned char w[], unsigned char u[], unsigned char v[], int n);
void mpdiv(unsigned char q[], unsigned char r[], unsigned char u[],
	unsigned char v[], int n, int m);
void mpinv(unsigned char u[], unsigned char v[], int n, int m);
void mplsh(unsigned char u[], int n);
void mpmov(unsigned char u[], unsigned char v[], int n);
void mpmul(unsigned char w[], unsigned char u[], unsigned char v[], int n,
	int m);
void mpneg(unsigned char u[], int n);
void mppi(int n);
void mprove(float **a, float **alud, int n, int indx[], float b[],
	float x[]);
void mpsad(unsigned char w[], unsigned char u[], int n, int iv);
void mpsdv(unsigned char w[], unsigned char u[], int n, int iv, int *ir);
void mpsmu(unsigned char w[], unsigned char u[], int n, int iv);
void mpsqrt(unsigned char w[], unsigned char u[], unsigned char v[], int n,
	int m);
void mpsub(int *is, unsigned char w[], unsigned char u[], unsigned char v[],
	int n);
void mrqcof(float x[], float y[], float sig[], int ndata, float a[],
	int ia[], int ma, float **alpha, float beta[], float *chisq,
	void (*funcs)(float, float [], float *, float [], int));
void mrqmin(float x[], float y[], float sig[], int ndata, float a[],
	int ia[], int ma, float **covar, float **alpha, float *chisq,
	void (*funcs)(float, float [], float *, float [], int), float *alamda);
void newt(float x[], int n, int *check,
	void (*vecfunc)(int, float [], float []));
void odeint(float ystart[], int nvar, float x1, float x2,
	float eps, float h1, float hmin, int *nok, int *nbad,
	void (*derivs)(float, float [], float []),
	void (*rkqs)(float [], float [], int, float *, float, float,
	float [], float *, float *, void (*)(float, float [], float [])));
void orthog(int n, float anu[], float alpha[], float beta[], float a[],
	float b[]);
void pade(double cof[], int n, float *resid);
void pccheb(float d[], float c[], int n);
void pcshft(float a, float b, float d[], int n);
void pearsn(float x[], float y[], unsigned long n, float *r, float *prob,
	float *z);
void period(float x[], float y[], int n, float ofac, float hifac,
	float px[], float py[], int np, int *nout, int *jmax, float *prob);
void piksr2(int n, float arr[], float brr[]);
void piksrt(int n, float arr[]);
void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
	float ***c, float **s);
float plgndr(int l, int m, float x);
float poidev(float xm, long *idum);
void polcoe(float x[], float y[], int n, float cof[]);
void polcof(float xa[], float ya[], int n, float cof[]);
void poldiv(float u[], int n, float v[], int nv, float q[], float r[]);
void polin2(float x1a[], float x2a[], float **ya, int m, int n,
	float x1, float x2, float *y, float *dy);
void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
void powell(float p[], float **xi, int n, float ftol, int *iter, float *fret,
	float (*func)(float []));
void predic(float data[], int ndata, float d[], int m, float future[], int nfut);
float probks(float alam);
void psdes(unsigned long *lword, unsigned long *irword);
void pwt(float a[], unsigned long n, int isign);
void pwtset(int n);
float pythag(float a, float b);
void pzextr(int iest, float xest, float yest[], float yz[], float dy[],
	int nv);
float qgaus(float (*func)(float), float a, float b);
void qrdcmp(float **a, int n, float *c, float *d, int *sing);
float qromb(float (*func)(float), float a, float b);
float qromo(float (*func)(float), float a, float b,
	float (*choose)(float (*)(float), float, float, int));
void qroot(float p[], int n, float *b, float *c, float eps);
void qrsolv(float **a, int n, float c[], float d[], float b[]);
void qrupdt(float **r, float **qt, int n, float u[], float v[]);
float qsimp(float (*func)(float), float a, float b);
float qtrap(float (*func)(float), float a, float b);
float quad3d(float (*func)(float, float, float), float x1, float x2);
void quadct(float x, float y, float xx[], float yy[], unsigned long nn,
	float *fa, float *fb, float *fc, float *fd);
void quadmx(float **a, int n);
void quadvl(float x, float y, float *fa, float *fb, float *fc, float *fd);
float ran0(long *idum);
float ran1(long *idum);
float ran2(long *idum);
float ran3(long *idum);
float ran4(long *idum);
void rank(unsigned long n, unsigned long indx[], unsigned long irank[]);
void ranpt(float pt[], float regn[], int n);
void ratint(float xa[], float ya[], int n, float x, float *y, float *dy);
void ratlsq(double (*fn)(double), double a, double b, int mm, int kk,
	double cof[], double *dev);
double ratval(double x, double cof[], int mm, int kk);
float rc(float x, float y);
float rd(float x, float y, float z);
void realft(float data[], unsigned long n, int isign);
void rebin(float rc, int nd, float r[], float xin[], float xi[]);
void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	int ic1, int jc1, int jcf, int kc, float ***c, float **s);
void relax(double **u, double **rhs, int n);
void relax2(double **u, double **rhs, int n);
void resid(double **res, double **u, double **rhs, int n);
float revcst(float x[], float y[], int iorder[], int ncity, int n[]);
void reverse(int iorder[], int ncity, int n[]);
float rf(float x, float y, float z);
float rj(float x, float y, float z, float p);
void rk4(float y[], float dydx[], int n, float x, float h, float yout[],
	void (*derivs)(float, float [], float []));
void rkck(float y[], float dydx[], int n, float x, float h,
	float yout[], float yerr[], void (*derivs)(float, float [], float []));
void rkdumb(float vstart[], int nvar, float x1, float x2, int nstep,
	void (*derivs)(float, float [], float []));
void rkqs(float y[], float dydx[], int n, float *x,
	float htry, float eps, float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []));
void rlft3(float ***data, float **speq, unsigned long nn1,
	unsigned long nn2, unsigned long nn3, int isign);
float rofunc(float b);
void rotate(float **r, float **qt, int n, int i, float a, float b);
void rsolv(float **a, int n, float d[], float b[]);
void rstrct(double **uc, double **uf, int nc);
float rtbis(float (*func)(float), float x1, float x2, float xacc);
float rtflsp(float (*func)(float), float x1, float x2, float xacc);
float rtnewt(void (*funcd)(float, float *, float *), float x1, float x2,
	float xacc);
float rtsafe(void (*funcd)(float, float *, float *), float x1, float x2,
	float xacc);
float rtsec(float (*func)(float), float x1, float x2, float xacc);
void rzextr(int iest, float xest, float yest[], float yz[], float dy[], int nv);
void savgol(float c[], int np, int nl, int nr, int ld, int m);
void score(float xf, float y[], float f[]);
void scrsho(float (*fx)(float));
float select(unsigned long k, unsigned long n, float arr[]);
float selip(unsigned long k, unsigned long n, float arr[]);
void shell(unsigned long n, float a[]);
void shoot(int n, float v[], float f[]);
void shootf(int n, float v[], float f[]);
void simp1(float **a, int mm, int ll[], int nll, int iabf, int *kp,
	float *bmax);
void simp2(float **a, int n, int l2[], int nl2, int *ip, int kp, float *q1);
void simp3(float **a, int i1, int k1, int ip, int kp);
void simplx(float **a, int m, int n, int m1, int m2, int m3, int *icase,
	int izrov[], int iposv[]);
void simpr(float y[], float dydx[], float dfdx[], float **dfdy,
	int n, float xs, float htot, int nstep, float yout[],
	void (*derivs)(float, float [], float []));
void sinft(float y[], int n);
void slvsm2(double **u, double **rhs);
void slvsml(double **u, double **rhs);
void sncndn(float uu, float emmc, float *sn, float *cn, float *dn);
double snrm(unsigned long n, double sx[], int itol);
void sobseq(int *n, float x[]);
void solvde(int itmax, float conv, float slowc, float scalv[],
	int indexv[], int ne, int nb, int m, float **y, float ***c, float **s);
void sor(double **a, double **b, double **c, double **d, double **e,
	double **f, double **u, int jmax, double rjac);
void sort(unsigned long n, float arr[]);
void sort2(unsigned long n, float arr[], float brr[]);
void sort3(unsigned long n, float ra[], float rb[], float rc[]);
void spctrm(FILE *fp, float p[], int m, int k, int ovrlap);
void spear(float data1[], float data2[], unsigned long n, float *d, float *zd,
	float *probd, float *rs, float *probrs);
void sphbes(int n, float x, float *sj, float *sy, float *sjp, float *syp);
void splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a);
void splin2(float x1a[], float x2a[], float **ya, float **y2a, int m, int n,
	float x1, float x2, float *y);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void spread(float y, float yy[], unsigned long n, float x, int m);
void sprsax(float sa[], unsigned long ija[], float x[], float b[],
	unsigned long n);
void sprsin(float **a, int n, float thresh, unsigned long nmax, float sa[],
	unsigned long ija[]);
void sprspm(float sa[], unsigned long ija[], float sb[], unsigned long ijb[],
	float sc[], unsigned long ijc[]);
void sprstm(float sa[], unsigned long ija[], float sb[], unsigned long ijb[],
	float thresh, unsigned long nmax, float sc[], unsigned long ijc[]);
void sprstp(float sa[], unsigned long ija[], float sb[], unsigned long ijb[]);
void sprstx(float sa[], unsigned long ija[], float x[], float b[],
	unsigned long n);
void stifbs(float y[], float dydx[], int nv, float *xx,
	float htry, float eps, float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []));
void stiff(float y[], float dydx[], int n, float *x,
	float htry, float eps, float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []));
void stoerm(float y[], float d2y[], int nv, float xs,
	float htot, int nstep, float yout[],
	void (*derivs)(float, float [], float []));
void svbksb(float **u, float w[], float **v, int m, int n, float b[],
	float x[]);
void svdcmp(float **a, int m, int n, float w[], float **v);
void svdfit(float x[], float y[], float sig[], int ndata, float a[],
	int ma, float **u, float **v, float w[], float *chisq,
	void (*funcs)(float, float [], int));
void svdvar(float **v, int ma, float w[], float **cvm);
void toeplz(float r[], float x[], float y[], int n);
void tptest(float data1[], float data2[], unsigned long n, float *t, float *prob);
void tqli(float d[], float e[], int n, float **z);
float trapzd(float (*func)(float), float a, float b, int n);
void tred2(float **a, int n, float d[], float e[]);
void tridag(float a[], float b[], float c[], float r[], float u[],
	unsigned long n);
float trncst(float x[], float y[], int iorder[], int ncity, int n[]);
void trnspt(int iorder[], int ncity, int n[]);
void ttest(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *t, float *prob);
void tutest(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *t, float *prob);
void twofft(float data1[], float data2[], float fft1[], float fft2[],
	unsigned long n);
void vander(double x[], double w[], double q[], int n);
void vegas(float regn[], int ndim, float (*fxn)(float [], float), int init,
	unsigned long ncall, int itmx, int nprn, float *tgral, float *sd,
	float *chi2a);
void voltra(int n, int m, float t0, float h, float *t, float **f,
	float (*g)(int, float), float (*ak)(int, int, float, float));
void wt1(float a[], unsigned long n, int isign,
	void (*wtstep)(float [], unsigned long, int));
void wtn(float a[], unsigned long nn[], int ndim, int isign,
	void (*wtstep)(float [], unsigned long, int));
void wwghts(float wghts[], int n, float h,
	void (*kermom)(double [], double ,int));
int zbrac(float (*func)(float), float *x1, float *x2);
void zbrak(float (*fx)(float), float x1, float x2, int n, float xb1[],
	float xb2[], int *nb);
float zbrent(float (*func)(float), float x1, float x2, float tol);
void zrhqr(float a[], int m, float rtr[], float rti[]);
float zriddr(float (*func)(float), float x1, float x2, float xacc);
void zroots(fcomplex a[], int m, fcomplex roots[], int polish);

#else /* ANSI */
/* traditional - K&R */

void addint();
void airy();
void amebsa();
void amoeba();
float amotry();
float amotsa();
void anneal();
double anorm2();
void arcmak();
void arcode();
void arcsum();
void asolve();
void atimes();
void avevar();
void balanc();
void banbks();
void bandec();
void banmul();
void bcucof();
void bcuint();
void beschb();
float bessi();
float bessi0();
float bessi1();
void bessik();
float bessj();
float bessj0();
float bessj1();
void bessjy();
float bessk();
float bessk0();
float bessk1();
float bessy();
float bessy0();
float bessy1();
float beta();
float betacf();
float betai();
float bico();
void bksub();
float bnldev();
float brent();
void broydn();
void bsstep();
void caldat();
void chder();
float chebev();
void chebft();
void chebpc();
void chint();
float chixy();
void choldc();
void cholsl();
void chsone();
void chstwo();
void cisi();
void cntab1();
void cntab2();
void convlv();
void copy();
void correl();
void cosft();
void cosft1();
void cosft2();
void covsrt();
void crank();
void cyclic();
void daub4();
float dawson();
float dbrent();
void ddpoly();
int decchk();
void derivs();
float df1dim();
void dfour1();
void dfpmin();
float dfridr();
void dftcor();
void dftint();
void difeq();
void dlinmin();
double dpythag();
void drealft();
void dsprsax();
void dsprstx();
void dsvbksb();
void dsvdcmp();
void eclass();
void eclazz();
float ei();
void eigsrt();
float elle();
float ellf();
float ellpi();
void elmhes();
float erfcc();
float erff();
float erffc();
void eulsum();
float evlmem();
float expdev();
float expint();
float f1();
float f1dim();
float f2();
float f3();
float factln();
float factrl();
void fasper();
void fdjac();
void fgauss();
void fill0();
void fit();
void fitexy();
void fixrts();
void fleg();
void flmoon();
float fmin();
void four1();
void fourew();
void fourfs();
void fourn();
void fpoly();
void fred2();
float fredin();
void frenel();
void frprmn();
void ftest();
float gamdev();
float gammln();
float gammp();
float gammq();
float gasdev();
void gaucof();
void gauher();
void gaujac();
void gaulag();
void gauleg();
void gaussj();
void gcf();
float golden();
void gser();
void hpsel();
void hpsort();
void hqr();
void hufapp();
void hufdec();
void hufenc();
void hufmak();
void hunt();
void hypdrv();
fcomplex hypgeo();
void hypser();
unsigned short icrc();
unsigned short icrc1();
int igray();
void iindexx();
void indexx();
void interp();
int irbit1();
int irbit2();
void jacobi();
void jacobn();
long julday();
void kendl1();
void kendl2();
void kermom();
void ks2d1s();
void ks2d2s();
void ksone();
void kstwo();
void laguer();
void lfit();
void linbcg();
void linmin();
void lnsrch();
void load();
void load1();
void load2();
void locate();
void lop();
void lubksb();
void ludcmp();
void machar();
void matadd();
void matsub();
void medfit();
void memcof();
int metrop();
void mgfas();
void mglin();
float midexp();
float midinf();
float midpnt();
float midsql();
float midsqu();
void miser();
void mmid();
void mnbrak();
void mnewt();
void moment();
void mp2dfr();
void mpadd();
void mpdiv();
void mpinv();
void mplsh();
void mpmov();
void mpmul();
void mpneg();
void mppi();
void mprove();
void mpsad();
void mpsdv();
void mpsmu();
void mpsqrt();
void mpsub();
void mrqcof();
void mrqmin();
void newt();
void odeint();
void orthog();
void pade();
void pccheb();
void pcshft();
void pearsn();
void period();
void piksr2();
void piksrt();
void pinvs();
float plgndr();
float poidev();
void polcoe();
void polcof();
void poldiv();
void polin2();
void polint();
void powell();
void predic();
float probks();
void psdes();
void pwt();
void pwtset();
float pythag();
void pzextr();
float qgaus();
void qrdcmp();
float qromb();
float qromo();
void qroot();
void qrsolv();
void qrupdt();
float qsimp();
float qtrap();
float quad3d();
void quadct();
void quadmx();
void quadvl();
float ran0();
float ran1();
float ran2();
float ran3();
float ran4();
void rank();
void ranpt();
void ratint();
void ratlsq();
double ratval();
float rc();
float rd();
void realft();
void rebin();
void red();
void relax();
void relax2();
void resid();
float revcst();
void reverse();
float rf();
float rj();
void rk4();
void rkck();
void rkdumb();
void rkqs();
void rlft3();
float rofunc();
void rotate();
void rsolv();
void rstrct();
float rtbis();
float rtflsp();
float rtnewt();
float rtsafe();
float rtsec();
void rzextr();
void savgol();
void score();
void scrsho();
float select();
float selip();
void shell();
void shoot();
void shootf();
void simp1();
void simp2();
void simp3();
void simplx();
void simpr();
void sinft();
void slvsm2();
void slvsml();
void sncndn();
double snrm();
void sobseq();
void solvde();
void sor();
void sort();
void sort2();
void sort3();
void spctrm();
void spear();
void sphbes();
void splie2();
void splin2();
void spline();
void splint();
void spread();
void sprsax();
void sprsin();
void sprspm();
void sprstm();
void sprstp();
void sprstx();
void stifbs();
void stiff();
void stoerm();
void svbksb();
void svdcmp();
void svdfit();
void svdvar();
void toeplz();
void tptest();
void tqli();
float trapzd();
void tred2();
void tridag();
float trncst();
void trnspt();
void ttest();
void tutest();
void twofft();
void vander();
void vegas();
void voltra();
void wt1();
void wtn();
void wwghts();
int zbrac();
void zbrak();
float zbrent();
void zrhqr();
float zriddr();
void zroots();

#endif /* ANSI */

#endif /* _NR_H_ */
                                                                                                                                                                                                                                                                                 nrutil.h                                                                                            0100600 0000764 0000144 00000006402 07675067760 011616  0                                                                                                    ustar   chris                           users                                                                                                                                                                                                                  #ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

void nrerror(char error_text[]);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);

#else /* ANSI */
/* traditional - K&R */

void nrerror();
float *vector();
float **matrix();
float **submatrix();
float **convert_matrix();
float ***f3tensor();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
unsigned char *cvector();
unsigned long *lvector();
void free_vector();
void free_dvector();
void free_ivector();
void free_cvector();
void free_lvector();
void free_matrix();
void free_submatrix();
void free_convert_matrix();
void free_dmatrix();
void free_imatrix();
void free_f3tensor();

#endif /* ANSI */

#endif /* _NR_UTILS_H_ */
                                                                                                                                                                                                                                                              tdoa.h                                                                                              0100640 0000764 0000144 00000002563 07675734145 011236  0                                                                                                    ustar   chris                           users                                                                                                                                                                                                                  /******************************************************************/
/* file name: tdoa.h                                              */
/******************************************************************/

/***********************************************/
/* defines                                     */
/***********************************************/
#define TDOA_SUCCESS               (0)
#define TDOA_NEGATIVE_DISCRIMINANT (-1)
#define TDOA_IJ_IS_ZERO            (1)
#define TDOA_IK_IS_ZERO            (2)
#define TDOA_KJ_IS_ZERO            (4)
#define TDOA_KL_IS_ZERO            (8)

#define TDOA_FT_PER_NS  (0.98357105643L)
#define TDOA_M_PER_NS   (0.299792458L)

/***********************************************/
/* typedefs                                    */
/***********************************************/
struct vector3
{
  int x;
  int y;
  int z;
};

struct vector4
{
  double x;
  double y;
  double z;
  double t;
};

/***********************************************/
/* function prototypes                         */
/***********************************************/
extern int computePosition(struct vector4 *viPtr,
                           struct vector4 *vjPtr,
                           struct vector4 *vkPtr,
                           struct vector4 *vlPtr,
                           struct vector3 *v1Ptr,
                           struct vector3 *v2Ptr);


                                                                                                                                             tdoa2d.h                                                                                            0100644 0000764 0000144 00000002054 07676077734 011470  0                                                                                                    ustar   chris                           users                                                                                                                                                                                                                  /******************************************************************/
/* file name: tdoa.h                                              */
/******************************************************************/

/***********************************************/
/* defines                                     */
/***********************************************/
#define TDOA_SUCCESS               (0)
#define TDOA_FAIL                  (-1)

#define TDOA_FT_PER_NS  (0.98357105643L)
#define TDOA_M_PER_NS   (0.299792458L)

/***********************************************/
/* typedefs                                    */
/***********************************************/
struct vector2
{
  double x;
  double y;
};


/***********************************************/
/* function prototypes                         */
/***********************************************/
extern int computePosition(struct vector2 *v1Ptr,
                           struct vector2 *v2Ptr,
                           struct vector2 *v3Ptr,
		           double t12,
		           double t13);


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    wls.h                                                                                               0100600 0000764 0000144 00000003012 07676241460 011070  0                                                                                                    ustar   chris                           users                                                                                                                                                                                                                  /*************************************************************************
 file name: wls.h     
*************************************************************************/

#ifndef __WLS_H__
#define __WLS_H__

/*************************************************************************
 defines    
*************************************************************************/
#define WLS_MAX_RECEIVERS   (50 + 1)
#define WLS_LAMDA_TOLERANCE (0.1)
#define WLS_MAX_ITERATIONS  (100)

#define TDOA_FT_PER_NS      (0.98357105643L)
#define TDOA_M_PER_NS       (0.299792458L)

/*************************************************************************
 structure definitions     
*************************************************************************/     
struct vector3
{
  float x;   /* x coordinate */
  float y;   /* y coordinate */
  float dt;  /* time difference with respect to reference */  
};

/*************************************************************************
 function prototypes     
*************************************************************************/
void wlsInit(void);

void wlsComputePosition(int m,float xn,float yn,struct vector3 *v,
                        float *x,float *y);

float wlsComputeNomRange(float xn,float xi,float x0,
                         float yn,float yi,float y0);
			 
void wlsComputeDirection(float xn,float xi,float x0,
                         float yn,float yi,float y0,
			 float *cx,float *cy);
			 
void wlsComputeDeltaXy(int m,float **c,float **dR,float **dXy);

#endif /* __WLS_H__ */
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      gendata.c                                                                                           0100644 0000764 0000144 00000002117 07676164116 011700  0                                                                                                    ustar   chris                           users                                                                                                                                                                                                                  #include <stdio.h>
#include <math.h>
#include "wls.h"

int main(int argc, char **argv)
{
  char buffer[80];
  int i, n;
  double x, y;
  double xo, yo;
  double deltax, deltay;
  double t[10];

  /* read number of receivers */
  gets(buffer);
  sscanf(buffer,"%d",&n);
       
  /* get position of source */
  gets(buffer);
  sscanf(buffer,"%lf %lf %lf %lf",&x,&y);

  /* indicate number of time differences */
  printf("%d\n",(n-1));
  
  /* indicate initial position estimate */
  printf("0 0\n");
  
  for (i = 0; i < n; i++)
  {
    /* get location of observer i */
    gets(buffer);
    sscanf(buffer,"%lf %lf",&xo,&yo);
  
    /* compute difference vector */
    deltax = x - xo;
    deltay = y - yo;
  
    /* compute distance */
    t[i] = (deltax * deltax) + (deltay * deltay);
    t[i] = sqrt(t[i]);
    
    /* convert to time */
    t[i] /= TDOA_FT_PER_NS;
        
    /* output results position */
    printf("%lf %lf\n",xo,yo);
  } /* for */

  /* compute time differences */
  for (i = 1; i < n; i++)
  {
    printf("%20.14lf\n",(t[i]-t[0]));
  } /* for */
    
  exit(0);

} /* main */
                                                                                                                                                                                                                                                                                                                                                                                                                                                 lubksb.c                                                                                            0100600 0000764 0000144 00000000622 07675104457 011547  0                                                                                                    ustar   chris                           users                                                                                                                                                                                                                  void lubksb(float **a, int n, int *indx, float b[])
{
	int i,ii=0,ip,j;
	float sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 89::9. */
                                                                                                              ludcmp.c                                                                                            0100600 0000764 0000144 00000002113 07675104457 011546  0                                                                                                    ustar   chris                           users                                                                                                                                                                                                                  #include <math.h>
#define NRANSI
#include "nrutil.h"
#define TINY 1.0e-20;

void ludcmp(float **a, int n, int *indx, float *d)
{
	int i,imax,j,k;
	float big,dum,sum,temp;
	float *vv;

	vv=vector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_vector(vv,1,n);
}
#undef TINY
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 89::9. */
                                                                                                                                                                                                                                                                                                                                                                                                                                                     matutil.c                                                                                           0100644 0000764 0000144 00000012673 07676241666 011772  0                                                                                                    ustar   chris                           users                                                                                                                                                                                                                  /*************************************************************************
 file name: matutil.c
*************************************************************************/

#include "nrutil.h"
#include "matutil.h"

static float amCopy[MAT_MAX_SIZE][MAT_MAX_SIZE];
static float *aCopy[MAT_MAX_SIZE]; 

/*************************************************************************

  Name: matInit
  
  Purpose: The purpose of this function is to initialize all static
  storage associated with the matrix system.
  
  Calling Sequence: matInit()
  
  Inputs:

      None.
      
  Outputs:
  
      None.
  
*************************************************************************/
void matInit(void)
{
  int i;

  for (i = 1; i < MAT_MAX_SIZE; i++)
  {
    aCopy[i] = amCopy[i];
  } /* for */
  
  return;
  
} /* matInit */

/*************************************************************************

  Name: matCopy
  
  Purpose: The purpose of this function is to copy the matrix a into b.
  
  Calling Sequence: matCopy(m,n,a,b)
  
  Inputs:

      m - The number of rows in matrix a and b.
      
      n - The number of columns in matrix a, and b.
      
      a - An array of pointers to the rows of a[1..m][1..n].
      
      b - An array of pointers to the rows of b[1..n][1..n], which
      represents the copy of a.
      
  Outputs:
  
      None.
  
*************************************************************************/
void matCopy(int m,int n,float **a,float **b)
{
  int i, j;
  
  for (i = 1; i <= m; i++)
  {
    for (j = 1; j <= n; j++)
    {
      b[i][j] = a[i][j];
    } /* for */
  } /* for */

  return;

} /* matAdd */

/*************************************************************************

  Name: matTranspose
  
  Purpose: The purpose of this function is to compute the transpose of a
  matrix.  
  
  Calling Sequence: matTranspose(m,n,a,aT)
  
  Inputs:

      m - The number of rows in a.
      
      n - The number of columns in a.
      
      a - An array of pointers to the rows of a[1..m][1..n].
        
      aT - An array of pointers to the rows of aT[1..n][1..m], which
      represents the transpose of a.
       
  Outputs:

      None.
         
*************************************************************************/
void matTranspose(int m,int n,float **a,float **aT)
{
  int i, j;
  
  for (i = 1; i <= m; i++)
  {
    for (j = 1; j <= n; j++)
    {
      aT[j][i] = a[i][j];
    } /* for */
  } /* for */
  
  return;

} /* matTranspose */

/*************************************************************************

  Name: matAdd
  
  Purpose: The purpose of this function is to add two matrices,
  a and b.
  
  Calling Sequence: matAdd(m,n,a,b,c)
  
  Inputs:

      m - The number of rows in matrix a and b.
      
      n - The number of columns in matrix a, and b.
      
      a - An array of pointers to the rows of a[1..m][1..n].
      
      b - An array of pointers to the rows of b[1..n][1..n].
      
      c - An array of pointers to the rows of c[1..n][1..n], which
      represents the result a + b.
       
  Outputs:
  
      None.
  
*************************************************************************/
void matAdd(int m,int n,float **a,float **b,float **c)
{
  int i, j;
  
  for (i = 1; i <= m; i++)
  {
    for (j = 1; j <= n; j++)
    {
      c[i][j] = a[i][j] + b[i][j];
    } /* for */
  } /* for */

  return;

} /* matAdd */

/*************************************************************************

  Name: matMultiply
  
  Purpose: The purpose of this function is to multiply two matrices,
  a and b.
  
  Calling Sequence: matMultiply(m,n,p,a,b,c)
  
  Inputs:

      m - The number of rows in matrix a.
      
      n - The number of columns in matrix a, and number of rows in
      matrix b.
      
      p - The number of columns in matrix b.

      a - An array of pointers to the rows of a[1..m][1..n].
      
      b - An array of pointers to the rows of b[1..n][1..p].
      
      c - An array of pointers to the rows of c[1..m][1..p], which
      represents the result a * b.
       
  Outputs:
  
      None.
  
*************************************************************************/
void matMultiply(int m,int n,int p,float **a,float **b,float **c)
{
  int i, j, k;
  
  for (i = 1; i <= m; i++)
  {
    for (j = 1; j <= p; j++)
    {
      c[i][j] = 0.0;
      
      for (k = 1; k <= n; k++)
      {
        c[i][j] += a[i][k] * b[k][j]; 
      } /* for */
    } /* for */ 
  } /* for */
  
  return;
  
} /* matMultiply */

/*************************************************************************

  Name: matInverse
  
  Purpose: The purpose of this function is to invert a matrix a.
  
  Calling Sequence: matInverse(n,a,aInv)
  
  Inputs:

      n - The number of rows and columns in matrix a.
      
      a - An array of pointers to the rows of a[1..n][1..n].
      
      aInv - An array of pointers to the rows of aInv[1..n][1..n], which
      represents the result of the inverse of a.
       
  Outputs:
  
      None.
  
*************************************************************************/
void matInverse(int n,float **a,float **aInv)
{
  float col[10], d;
  int i, j, indx[10];


  matCopy(n,n,a,aCopy);
  
  ludcmp(aCopy,n,&indx,&d);
  
  for(j = 1; j <= n; j++)
  {
    for (i = 1; i <= n; i++)
    {
      col[i] = 0.0;
    } /* for */
    
    col[j] = 1.0;
    lubksb(aCopy,n,&indx,col);
    
    for (i = 1; i <= n; i++)
    {
      aInv[i][j] = col[i];
    } /* for */
  } /* for */
  
  return;
  
} /* matInverse */


                                                                     nrutil.c                                                                                            0100600 0000764 0000144 00000037515 07675071476 011621  0                                                                                                    ustar   chris                           users                                                                                                                                                                                                                  #if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	float **m;

	/* allocate array of pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

#else /* ANSI */
/* traditional - K&R */

#include <stdio.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(nl,nh)
long nh,nl;
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(nl,nh)
long nh,nl;
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(nl,nh)
long nh,nl;
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((unsigned int)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((unsigned int)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((unsigned int)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
float **a;
long newcl,newrl,oldch,oldcl,oldrh,oldrl;
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	float **m;

	/* allocate array of pointers to rows */
	m=(float **) malloc((unsigned int) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(a,nrl,nrh,ncl,nch)
float *a;
long nch,ncl,nrh,nrl;
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((unsigned int) ((nrow+NR_END)*sizeof(float*)));
	if (!m)	nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

float ***f3tensor(nrl,nrh,ncl,nch,ndl,ndh)
long nch,ncl,ndh,ndl,nrh,nrl;
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((unsigned int)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((unsigned int)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(v,nl,nh)
float *v;
long nh,nl;
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(v,nl,nh)
int *v;
long nh,nl;
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(v,nl,nh)
long nh,nl;
unsigned char *v;
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(v,nl,nh)
long nh,nl;
unsigned long *v;
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(v,nl,nh)
double *v;
long nh,nl;
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
long nch,ncl,nrh,nrl;
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
long nch,ncl,nrh,nrl;
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
long nch,ncl,nrh,nrl;
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(b,nrl,nrh,ncl,nch)
float **b;
long nch,ncl,nrh,nrl;
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(b,nrl,nrh,ncl,nch)
float **b;
long nch,ncl,nrh,nrl;
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(t,nrl,nrh,ncl,nch,ndl,ndh)
float ***t;
long nch,ncl,ndh,ndl,nrh,nrl;
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

#endif /* ANSI */
                                                                                                                                                                                   wls.c                                                                                               0100644 0000764 0000144 00000024563 07676244023 011107  0                                                                                                    ustar   chris                           users                                                                                                                                                                                                                  #include <stdio.h>
#include <math.h>
#include "matutil.h"
#include "wls.h"

float dXym[3][2], dRm[WLS_MAX_RECEIVERS][2];
float *dXy[3], *dR[WLS_MAX_RECEIVERS];

float cm[WLS_MAX_RECEIVERS][3];
float *c[WLS_MAX_RECEIVERS];

float cTm[3][WLS_MAX_RECEIVERS], cTxcm[3][3], cTxcInvm[3][3], cTxdRm[3][2];
float *cT[3], *cTxc[3], *cTxcInv[3], *cTxdR[3];

/*************************************************************************

  Name: main program
  
  Purpose: The purpose of this program is to compute the position of
  an RF emitter, in 2 space, using the time difference of arrivial (TDOA)
  technique.  The position parameters have units of feet, and
  the time parameters have units of nanoseconds.
  
  Calling Sequence: wls < inputFile
  
  Inputs:
    
    A text file with the following format:
      m
      xe ye  
      x0 y0
      x1 y1
      x2 y2
        .
	.
      xm-1 ym-1
      t01
      t02
       .
      t0m-1
      	     
      The tuples, [xi yi]:i=0,1,2,...m-1 correspond to the locations
      of receivers, and the quantities t0i:i=0,1,2,...m-1, correspond
      to the differences of time of arrival between the reference
      receiver, located at (x0,y0), and a receiver at location (xi,yi).
      Note that the receiver, located at (x0,y0) is designated as the
      reference receiver for which all time differences are computed.
      The tuple [xe ye] is the initial estimate of the position of the
      emitter.      
  
  Outputs:
  
     Printed output of two solutions [x1 y1 z1], [x2 y2 z2], which locate
     the RF emitter. 

*************************************************************************/
int main(void)
{
  char buffer[80];
  struct vector3 v[WLS_MAX_RECEIVERS];  
  int i, m;
  float xe, ye;
  float x, y;
  double distance, theta;
  
  /* initialize the system */
  wlsInit();  

  /* read number of receivers */  
  gets(buffer);
  sscanf(buffer,"%d",&m);
  
  /* read initial estimate of receiver position */
  gets(buffer);
  sscanf(buffer,"%f %f\n",&xe,&ye);
  
  /* read receiver locations */
  for (i = 0; i <= m; i++)
  {
    gets(buffer);
    sscanf(buffer,"%f %f\n",&v[i].x,&v[i].y);
  } /* for */
      
  /* read receiver time differences with respect to reference receiver */
  for (i = 1; i <= m; i++)
  {
    gets(buffer);
    sscanf(buffer,"%f %f\n",&v[i].dt);
  } /* for */

  /* determine location of emitter */
  wlsComputePosition(m,xe,ye,v,&x,&y);

  /* output in cartesian coordinates */
  printf("\n[x y] = [%f %f]\n",x,y);

  /* compute distance */
  distance = sqrt((x*x) + (y*y));

  /* fix to avoid division by zero */
  if (x == 0)
  {
    x = 0.001;
  } /* if */
    
  /* compute bearing */
  theta = atan(y / x) * 180 / M_PI;
    
  /* output in polar coordinates */
  printf("[r theta] = [%f %f]\n",distance,theta);

  exit(0);
  
} /* main */

/*************************************************************************

  Name: wlsInit
  
  Purpose: The purpose of this function is to initialize all data
  structures.
  
  Calling Sequence: wlsInit()
  
  Inputs:

      None.
      
  Outputs:
  
      None.
  
*************************************************************************/
void wlsInit(void)
{
  int i, j;

  /* initialize the matrix system */
  matInit();
  
  for (i = 0; i < 3; i++)
  {
    dXy[i] = dXym[i];
  } /* for */
  
  for (i = 0; i < WLS_MAX_RECEIVERS; i++)
  {
    dR[i] = dRm[i];
    c[i] = cm[i];
  } /* for */

  for (i = 0; i < 3; i++)
  {
    cT[i] = cTm[i];
    cTxc[i] = cTxcm[i];
    cTxcInv[i] = cTxcInvm[i];
    cTxdR[i] = cTxdRm[i];    
  } /* for */

  return;
  
} /* wlsInit */

/*************************************************************************

  Name: wlsComputePosition
  
  Purpose: The purpose of this function is to compute the position of
  an emitter.  A least squares method is used so that RF multipath,
  timing jitter, and receiver position errors can be dealt with.
    
  Calling Sequence: wlsComputePosition(m,xn,yn,v,x,y)
  
  Inputs:

      m - The number of time differences.
      
      xn - The initial x position estimate of the emitter.
      
      yn - The nominal y position best estimate of the emitter.
      
      v - A 3 vector that represents the position of a receiver and the
      time difference of arrival between the receiver and the reference
      receiver.  Note that v[0] represents the coordinates of the
      reference receiver.
      
      x - A pointer to storage of the returned x coordinate of the
      emitter.
      
      y - A pointer to storage of the returned y coordinate of the
      emitter.
                 
  Outputs:
  
      None.
  
*************************************************************************/
void wlsComputePosition(int m,float xn,float yn,struct vector3 *v,
                        float *x,float *y)
{
  int i, done, loopCount;
  float r0i[10], rn0i[10];
  float dx, dy, lamda;
    
  /* used for while loop */
  done = 0;
  loopCount = 0;
  
  /* compute r0i */
  for (i = 1; i <=m; i++)
  {
    /* convert to ft */
    r0i[i] = (100000 * v[i].dt) / 101670;
  } /* for */
  
  while (!done)
  {
    for (i = 1; i <= m; i++)
    {
      /* compute differential direction cosine matrix */
      wlsComputeDirection(xn,v[i].x,v[0].x,yn,
                          v[i].y,v[0].y,&c[i][1],&c[i][2]);    
 
      /* compute nominal differential pseudo range */
      rn0i[i] = wlsComputeNomRange(xn,v[i].x,v[0].x,yn,v[i].y,v[0].y); 

      /* compute incremental differential pseudo range */
      dR[i][1] = r0i[i] - rn0i[i];
    } /* for */
   
    /* solve the least squares equation for dx and dy */
    wlsComputeDeltaXy(m,c,dR,dXy);  
    
    /* compute new nominal position location */
    xn += dXy[1][1];
    yn += dXy[2][1];

    /* compute distance increment */
    dx = dXy[1][1];
    dy = dXy[2][1];
    lamda = sqrt((dx * dx) + (dy * dy));    

    /* break out of loop when desired resolution is met */
    if (lamda < WLS_LAMDA_TOLERANCE)
    {
      /* set return values */
      *x = xn;
      *y = yn;      
      /* bail out */
      break;
    } /* if */ 
     
    /* indicate that one more iteration occurred */
    loopCount++;
    
    /* break out of loop if too many iterations have occurred */
    if (loopCount >= WLS_MAX_ITERATIONS)
    {
      /* indicate failure */
      printf("too many iterations\n");      
      /* bail out */
      break;
    } /* if */ 
  } /* while */
  
  return;
  
} /* wlsComputePosition */

/*************************************************************************

  Name: wlsComputeNomRange
  
  Purpose: The purpose of this function is to compute the nominal
  differential pseudo range.
  
  Calling Sequence: rn0i = wlsComputeNomRange(xn,xi,x0,yn,yi.y0)
  
  Inputs:

     xn - The nominal x position best estimate of the emitter.
      
      xi - The x position of receiver i.
      
      x0 - The x position of the reference receiver.
      
      yn - The nominal y position best estimate of the emitter.
      
      yi - The y position of receiver i.
      
      y0 - The y position of the reference receiver.
           
  Outputs:
  
      rn0i - The nominal differential pseudo range.
  
*************************************************************************/
float wlsComputeNomRange(float xn,float xi,float x0,
                         float yn,float yi,float y0)
{
  float xni, xn0;
  float yni, yn0;
  float dni, dn0;
  float rn0i;

  /* compute distance differences */  
  xni = xn - xi;
  xn0 = xn - x0;
  yni = yn - yi;
  yn0 = yn - y0;
  
  /* compute distance from receiver i */
  dni = sqrt((xni * xni) + (yni * yni));
  
  /* compute distance from reference receiver */
  dn0 = sqrt((xn0 * xn0) + (yn0 * yn0));
  
  /* compute pseudo range */
  rn0i = dni - dn0;
    
  return (rn0i);
  
} /* wlsComputeNomRange */
				  
/*************************************************************************

  Name: wlsComputeDirection
  
  Purpose: The purpose of this function is to compute the differential
  cosine and differential sine associated with a particular receiver
  position with respect to the reference receiver and the emitter.
  
  Calling Sequence: wlsComputeDirection(xn,xi,x0,yn,yi,y0,cx,cy)
  
  Inputs:

      xn - The nominal x position best estimate of the emitter.
      
      xi - The x position of receiver i.
      
      x0 - The x position of the reference receiver.
      
      yn - The nominal y position best estimate of the emitter.
      
      yi - The y position of receiver i.
      
      y0 - The y position of the reference receiver.

      cx - A pointer to the direction cosine associated with the receiver.
      
      cy - A pointer to the direction sine associated with the reciever.
  
                    
  Outputs:
  
      None.
  
*************************************************************************/
void wlsComputeDirection(float xn,float xi,float x0,
                         float yn,float yi,float y0,
			 float *cx,float *cy)
{
  float xni, xn0;
  float yni, yn0;
  float dni, dn0;

  /* compute distance differences */  
  xni = xn - xi;
  xn0 = xn - x0;
  yni = yn - yi;
  yn0 = yn - y0;
  
  /* compute distance from receiver i */
  dni = sqrt((xni * xni) + (yni * yni));
  
  /* compute distance from reference receiver */
  dn0 = sqrt((xn0 * xn0) + (yn0 * yn0));
  
  /* compute direction cosine */
  *cx = (xni / dni) - (xn0 / dn0);
  
  /* compute direction sine */
  *cy = (yni / dni) - (yn0 / dn0);
  
  return;
  
} /* wlsComputeDirection */

/*************************************************************************

  Name: wlsComputeDeltaXy
  
  Purpose: The purpose of this function is to compute the incremental
  position correction to a position estimate.
  
  Calling Sequence: wlsComputeDeltaXy(m,c,dR,dXy)
  
  Inputs:

      m - The number of time differences.
      
      c - An array of pointers to the rows of c[1..m][1..2].
      
      dR - An array of pointers to the rows of dR[1..m][1..1].
      
      dXy - An array of pointers to the rows of dXy[1..2][1..1].
     
      
  Outputs:
  
      None.
  
*************************************************************************/
void wlsComputeDeltaXy(int m,float **c,float **dR,float **dXy)
{

  /* compute C' */
  matTranspose(m,2,c,cT);
  
  /* compute cTxc = C'.C */
  matMultiply(2,m,2,cT,c,cTxc);
  
  /* compute cInv = Inv(C'.C) */
  matInverse(2,cTxc,cTxcInv);
  
  /* compute [c'][dr] = C'.dR */
  matMultiply(2,m,1,cT,dR,cTxdR);
  
  /* compute dXy = Inv(C'.C).C'.dR */
  matMultiply(2,2,1,cTxcInv,cTxdR,dXy);
  
  return;

} /* wlsComputeDeltaXy */
                                                                                                                                             d1.inp                                                                                              0100640 0000764 0000144 00000000070 07676243750 011137  0                                                                                                    ustar   chris                           users                                                                                                                                                                                                                  4
-200 -300
1000 1000
-1000 -1000
-1000 1000
1000 -1000
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        d2.inp                                                                                              0100644 0000764 0000144 00000000126 07676146036 011145  0                                                                                                    ustar   chris                           users                                                                                                                                                                                                                  10
80 85
100 100
-100 -100
100 -100
30 30
50 50
-253 -125
-400 155
20 20
400 0
-400 0
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          