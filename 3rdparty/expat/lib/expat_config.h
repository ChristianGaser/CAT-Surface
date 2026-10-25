/* expat_config.h for CAT-Surface.
 *
 * Expat is compiled straight into libCAT, without its own configure script,
 * so this file replaces the one configure would generate.  It has to work on
 * every platform CAT-Surface is built for (macOS, Linux, MinGW, Emscripten),
 * which is why it decides from compiler macros only:
 *
 * - The byte order comes from the compiler.
 * - The hash salt comes from /dev/urandom on every POSIX system (Emscripten
 *   emulates it); Windows uses rand_s(), which xmlparse.c selects itself for
 *   _WIN32.  Expat falls back to time and process id if a read fails.
 * - General entities, DTDs and namespaces are on, with the usual context of
 *   1024 bytes, as in a default upstream build.
 */
#ifndef EXPAT_CONFIG_H
#define EXPAT_CONFIG_H 1

#if defined(__BYTE_ORDER__) && defined(__ORDER_BIG_ENDIAN__) \
    && (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
#  define BYTEORDER 4321
#  define WORDS_BIGENDIAN 1
#else
#  define BYTEORDER 1234
#endif

#define XML_CONTEXT_BYTES 1024
#define XML_DTD 1
#define XML_GE 1
#define XML_NS 1

#if ! defined(_WIN32)
#  define XML_DEV_URANDOM 1
#endif

#endif /* EXPAT_CONFIG_H */
