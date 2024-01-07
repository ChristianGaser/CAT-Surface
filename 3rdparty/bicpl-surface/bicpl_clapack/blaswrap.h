/* CLAPACK 3.0 BLAS wrapper macros
 * Feb 5, 2000
 */

//#define NO_BLAS_WRAP

#ifndef __BLASWRAP_H
#define __BLASWRAP_H

#ifndef NO_BLAS_WRAP

#define dbdsqr_ bicpl_dbdsqr_
#define dgebd2_ bicpl_dgebd2_
#define dgebrd_ bicpl_dgebrd_
#define dgelq2_ bicpl_dgelq2_
#define dgelqf_ bicpl_dgelqf_
#define dgeqr2_ bicpl_dgeqr2_
#define dgeqrf_ bicpl_dgeqrf_
#define dlabrd_ bicpl_dlabrd_
#define dlacpy_ bicpl_dlacpy_
#define dlamch_ bicpl_dlamch_
#define dlange_ bicpl_dlange_
#define dlapy2_ bicpl_dlapy2_
#define dlarfb_ bicpl_dlarfb_
#define dlarf_  bicpl_dlarf_
#define dlarfg_ bicpl_dlarfg_
#define dlarft_ bicpl_dlarft_
#define dlartg_ bicpl_dlartg_
#define dlas2_  bicpl_dlas2_
#define dlascl_ bicpl_dlascl_
#define dlaset_ bicpl_dlaset_
#define dlasq1_ bicpl_dlasq1_
#define dlasq2_ bicpl_dlasq2_
#define dlasq3_ bicpl_dlasq3_
#define dlasq4_ bicpl_dlasq4_
#define dlasq5_ bicpl_dlasq5_
#define dlasq6_ bicpl_dlasq6_
#define dlasr_  bicpl_dlasr_
#define dlasrt_ bicpl_dlasrt_
#define dlassq_ bicpl_dlassq_
#define dlasv2_ bicpl_dlasv2_
#define dorg2r_ bicpl_dorg2r_
#define dorgbr_ bicpl_dorgbr_
#define dorgl2_ bicpl_dorgl2_
#define dorglq_ bicpl_dorglq_
#define dorgqr_ bicpl_dorgqr_
#define dorm2r_ bicpl_dorm2r_
#define dormbr_ bicpl_dormbr_
#define dorml2_ bicpl_dorml2_
#define dormlq_ bicpl_dormlq_
#define dormqr_ bicpl_dormqr_
#define d_sign  bicpl_d_sign
#define exit_   bicpl_exit_
#define F77_aloc bicpl_F77_aloc
#define ieeeck_ bicpl_ieeeck_
#define ilaenv_ bicpl_ilaenv_
#define lsame_  bicpl_lsame_
#define pow_dd  bicpl_pow_dd
#define pow_di  bicpl_pow_di
#define s_cat   bicpl_s_cat
#define s_cmp   bicpl_s_cmp
#define s_copy  bicpl_s_copy
#define xerbla_ bicpl_xerbla_
 
/* BLAS1 routines */
#define drotg_  bicpl_drotg
#define drotmg_ bicpl_drotmg
#define drot_   bicpl_drot
#define drotm_  bicpl_drotm
#define dswap_  bicpl_dswap
#define dscal_  bicpl_dscal
#define zdscal_ bicpl_zdscal
#define dcopy_  bicpl_dcopy
#define dnrm2_  bicpl_dnrm2
 
/* BLAS2 routines */
#define dgemv_ bicpl_dgemv
#define dtrmv_ bicpl_dtrmv
#define dger_  bicpl_dger
 
/* BLAS3 routines */
#define dgemm_ bicpl_dgemm
#define dtrmm_ bicpl_dtrmm

#endif /* NO_BLAS_WRAP */

#endif /* __BLASWRAP_H */
