/* glibc_compat.h
 *
 * Keep the Linux binaries of the "Build Binaries" workflow runnable on an older
 * glibc than the one of the build container.
 *
 * An ELF binary records, for every symbol it imports, the glibc symbol version
 * it was linked against, and the dynamic loader treats that as a hard minimum.
 * The binaries are built in manylinux_2_28 (glibc 2.28), and calling expf()
 * there binds expf@GLIBC_2.27, so CAT_VolMarchingCubes and CAT_VolThicknessPbt
 * failed on every older system with
 *
 *   /lib64/libm.so.6: version `GLIBC_2.27' not found
 *
 * which covers RHEL/CentOS 7, Ubuntu 16.04 and SLES 12 -- still ordinary
 * cluster installations, on which the other CAT tools and the CAT mex-files
 * (see cat_glibc_compat.h in CAT) run. Nothing in CAT-Surface needs anything
 * new; the whole incompatibility comes from changes inside glibc itself:
 *
 *   expf, exp2f, logf, log2f, powf
 *       glibc 2.27 added new implementations under a new version.
 *
 *   exp, exp2, log, log2, pow
 *       glibc 2.29 did the same for the double versions (not reached in
 *       manylinux_2_28, but pinned so that a newer build image does not raise
 *       the floor again).
 *
 *   pthread_create, pthread_join, pthread_exit
 *       glibc 2.34 merged libpthread into libc and gave these a new version
 *       (again only relevant for a newer build image).
 *
 * glibc keeps the previous versions as compat symbols indefinitely, so asking
 * for them by name is no downgrade: on a current glibc they call the same
 * implementation and only add the old errno handling.
 *
 * The workflow injects this header into every translation unit with
 * "-include", because .symver rewrites the *reference* and has to be seen by
 * each source that makes one of these calls. It is only meant for the
 * dynamically linked CI build: a static link (as used for local
 * cross-compilation with -static) cannot resolve versioned references.
 *
 * What remains is memcpy@GLIBC_2.14 (x86-64), and the workflow checks this
 * floor for every binary.
 *
 * Define CAT_NO_GLIBC_COMPAT to switch this off.
 */

#ifndef CAT_GLIBC_COMPAT_H
#define CAT_GLIBC_COMPAT_H

#if defined(__linux__) && !defined(CAT_NO_GLIBC_COMPAT)

/* __GLIBC__ is only defined once a libc header has been seen. */
#include <limits.h>

#if defined(__GLIBC__) && defined(__ELF__)

/* The baseline version node is the one that architecture's ABI started at, so
   it exists in every glibc able to run the binary at all: x86-64 entered at
   2.2.5, aarch64 only at 2.17. On anything else the header does nothing. */
#if defined(__x86_64__) && !defined(__ILP32__)
#define CAT_GLIBC_BASE "GLIBC_2.2.5"
#elif defined(__aarch64__)
#define CAT_GLIBC_BASE "GLIBC_2.17"
#endif

#ifdef CAT_GLIBC_BASE

/* A .symver for a symbol the translation unit never calls emits nothing. */
__asm__(".symver expf,expf@" CAT_GLIBC_BASE);
__asm__(".symver exp2f,exp2f@" CAT_GLIBC_BASE);
__asm__(".symver logf,logf@" CAT_GLIBC_BASE);
__asm__(".symver log2f,log2f@" CAT_GLIBC_BASE);
__asm__(".symver powf,powf@" CAT_GLIBC_BASE);
__asm__(".symver exp,exp@" CAT_GLIBC_BASE);
__asm__(".symver exp2,exp2@" CAT_GLIBC_BASE);
__asm__(".symver log,log@" CAT_GLIBC_BASE);
__asm__(".symver log2,log2@" CAT_GLIBC_BASE);
__asm__(".symver pow,pow@" CAT_GLIBC_BASE);
__asm__(".symver pthread_create,pthread_create@" CAT_GLIBC_BASE);
__asm__(".symver pthread_join,pthread_join@" CAT_GLIBC_BASE);
__asm__(".symver pthread_exit,pthread_exit@" CAT_GLIBC_BASE);

#endif /* CAT_GLIBC_BASE */

#endif /* __GLIBC__ && __ELF__ */

#endif /* __linux__ && !CAT_NO_GLIBC_COMPAT */

#endif /* CAT_GLIBC_COMPAT_H */
