/* src/include/config.h.  Generated from config.h.in by configure.  */
/* src/include/config.h.in.  Generated from configure.ac by autoheader.  */

#ifndef R_CONFIG_H
#define R_CONFIG_H

/* Version of C Compiler */
#define CC_VER "gcc (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0"

/* Define to 1 if using 'alloca.c'. */
/* #undef C_ALLOCA */

/* C stack direction: 1 (down) or -1 (up) */
#define C_STACK_DIRECTION 1

/* Define to 1 if translation of program messages to the user's native
   language is requested. */
#define ENABLE_NLS 1

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* C type used for Fortran character lengths */
#define FC_LEN_T size_t

/* Version of Fortran Compiler */
#define FC_VER "GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0"

/* Define to 1 if you have the `access' function. */
#define HAVE_ACCESS 1

/* Define to 1 if you have 'alloca', as a function or macro. */
#define HAVE_ALLOCA 1

/* Define to 1 if <alloca.h> works. */
#define HAVE_ALLOCA_H 1

/* Define if you have the Aqua headers and libraries, and want to include
   support for R.app and for the quartz() device to be built. */
/* #undef HAVE_AQUA */

/* Define to 1 if you have the `argz_count' function. */
#define HAVE_ARGZ_COUNT 1

/* Define to 1 if you have the <argz.h> header file. */
#define HAVE_ARGZ_H 1

/* Define to 1 if you have the `argz_next' function. */
#define HAVE_ARGZ_NEXT 1

/* Define to 1 if you have the `argz_stringify' function. */
#define HAVE_ARGZ_STRINGIFY 1

/* Define to 1 if you have the <arpa/inet.h> header file. */
#define HAVE_ARPA_INET_H 1

/* Define to 1 if you have the `asprintf' function. */
#define HAVE_ASPRINTF 1

/* Define to 1 if you have the `atan2pi' function. */
/* #undef HAVE_ATAN2PI */

/* Define to 1 if you have the `atanpi' function. */
/* #undef HAVE_ATANPI */

/* Define to 1 if the compiler understands __builtin_expect. (For intl) */
#define HAVE_BUILTIN_EXPECT 1

/* Define to 1 if you have the <bzlib.h> header file. */
#define HAVE_BZLIB_H 1

/* Define to 1 if you have the `cabs' function. */
#define HAVE_CABS 1

/* Define to 1 if you have the `cacos' function. */
#define HAVE_CACOS 1

/* Define to 1 if you have cairo-ps. */
/* #undef HAVE_CAIRO_PDF */

/* Define to 1 if you have cairo-pdf. */
/* #undef HAVE_CAIRO_PS */

/* Define to 1 if you have cairo-svg. */
/* #undef HAVE_CAIRO_SVG */

/* Define to 1 if you have the `carg' function. */
#define HAVE_CARG 1

/* Define to 1 if you have the `casin' function. */
#define HAVE_CASIN 1

/* Define to 1 if you have the `catan' function. */
#define HAVE_CATAN 1

/* Define to 1 if you have the `ccos' function. */
#define HAVE_CCOS 1

/* Define to 1 if you have the `ccosh' function. */
#define HAVE_CCOSH 1

/* Define to 1 if you have the `cexp' function. */
#define HAVE_CEXP 1

/* Define to 1 if you have the OS X function CFLocaleCopyCurrent in the
   CoreFoundation framework. (For intl) */
/* #undef HAVE_CFLOCALECOPYCURRENT */

/* Define to 1 if you have the OS X function CFPreferencesCopyAppValue in the
   CoreFoundation framework. (For intl) */
/* #undef HAVE_CFPREFERENCESCOPYAPPVALUE */

/* Define to 1 if you have the `chdir' function. */
#define HAVE_CHDIR 1

/* Define to 1 if you have the `chmod' function. */
#define HAVE_CHMOD 1

/* Define to 1 if you have the `clock_gettime' function. */
#define HAVE_CLOCK_GETTIME 1

/* Define to 1 if you have the `clog' function. */
#define HAVE_CLOG 1

/* Defined if framework CoreFoundation is present */
/* #undef HAVE_COREFOUNDATION_FW */

/* Define to 1 if you have the `cospi' function. */
/* #undef HAVE_COSPI */

/* Define to 1 if you have the `cpow' function. */
#define HAVE_CPOW 1

/* Define to 1 if you have the `csin' function. */
#define HAVE_CSIN 1

/* Define to 1 if you have the `csinh' function. */
#define HAVE_CSINH 1

/* Define to 1 if you have the `csqrt' function. */
#define HAVE_CSQRT 1

/* Define to 1 if you have the `ctan' function. */
#define HAVE_CTAN 1

/* Define to 1 if you have the `ctanh' function. */
#define HAVE_CTANH 1

/* Define to 1 if you have the <curl/curl.h> header file. */
#define HAVE_CURL_CURL_H 1

/* Define if the GNU dcgettext() function is already present or preinstalled.
   */
#define HAVE_DCGETTEXT 1

/* Define to 1 if you have the declaration of `alloca', and to 0 if you don't.
   */
#define HAVE_DECL_ALLOCA 1

/* Define to 1 if you have the declaration of `dladdr', and to 0 if you don't.
   */
#define HAVE_DECL_DLADDR 1

/* Define to 1 if you have the declaration of `dlsym', and to 0 if you don't.
   */
#define HAVE_DECL_DLSYM 1

/* Define to 1 if you have the declaration of `feof_unlocked', and to 0 if you
   don't. (For intl) */
#define HAVE_DECL_FEOF_UNLOCKED 1

/* Define to 1 if you have the declaration of `fgets_unlocked', and to 0 if
   you don't. (For intl) */
#define HAVE_DECL_FGETS_UNLOCKED 1

/* Define to 1 if you have the declaration of `getc_unlocked', and to 0 if you
   don't. (For intl) */
#define HAVE_DECL_GETC_UNLOCKED 1

/* Define to 1 if you have the declaration of `isfinite', and to 0 if you
   don't. */
#define HAVE_DECL_ISFINITE 1

/* Define to 1 if you have the declaration of `isnan', and to 0 if you don't.
   */
#define HAVE_DECL_ISNAN 1

/* Define to 1 if you have the declaration of `mkdtemp', and to 0 if you
   don't. */
#define HAVE_DECL_MKDTEMP 1

/* Define to 1 if you have the declaration of `putenv', and to 0 if you don't.
   */
#define HAVE_DECL_PUTENV 1

/* Define to 1 if you have the declaration of `realpath', and to 0 if you
   don't. */
#define HAVE_DECL_REALPATH 1

/* Define to 1 if you have the declaration of `RTLD_DEFAULT', and to 0 if you
   don't. */
#define HAVE_DECL_RTLD_DEFAULT 1

/* Define to 1 if you have the declaration of `RTLD_NEXT', and to 0 if you
   don't. */
#define HAVE_DECL_RTLD_NEXT 1

/* Define to 1 if you have the declaration of `siglongjmp', and to 0 if you
   don't. */
#define HAVE_DECL_SIGLONGJMP 1

/* Define to 1 if you have the declaration of `sigsetjmp', and to 0 if you
   don't. */
#define HAVE_DECL_SIGSETJMP 1

/* Define to 1 if you have the declaration of `SIZE_MAX', and to 0 if you
   don't. */
#define HAVE_DECL_SIZE_MAX 1

/* Define to 1 if you have the declaration of `strdup', and to 0 if you don't.
   */
#define HAVE_DECL_STRDUP 1

/* Define to 1 if you have the declaration of `strncasecmp', and to 0 if you
   don't. */
#define HAVE_DECL_STRNCASECMP 1

/* Define to 1 if you have the declaration of `vasprintf', and to 0 if you
   don't. */
#define HAVE_DECL_VASPRINTF 1

/* Define to 1 if you have the declaration of `_snprintf', and to 0 if you
   don't. (For intl) */
#define HAVE_DECL__SNPRINTF 0

/* Define to 1 if you have the declaration of `_snwprintf', and to 0 if you
   don't. (For intl) */
#define HAVE_DECL__SNWPRINTF 0

/* Define to 1 if you have the <dirent.h> header file, and it defines `DIR'.
   */
#define HAVE_DIRENT_H 1

/* Define to 1 if you have the `dladdr' function. */
#define HAVE_DLADDR 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the `dlsym' function. */
#define HAVE_DLSYM 1

/* Define to 1 if you have the <elf.h> header file. */
#define HAVE_ELF_H 1

/* Define to 1 if you have the <errno.h> header file. */
#define HAVE_ERRNO_H 1

/* Define to 1 if you have the `execv' function. */
#define HAVE_EXECV 1

/* Define to 1 if you have the `exp10' function. */
#define HAVE_EXP10 1

/* Define to 1 if you have the `expm1' function. */
#define HAVE_EXPM1 1

/* Define if your Fortran compiler appends an extra_underscore to external
   names containing an underscore. */
/* #undef HAVE_F77_EXTRA_UNDERSCORE */

/* Define if your Fortran compiler appends an underscore to external names. */
#define HAVE_F77_UNDERSCORE 1

/* Define to 1 if you have the `fcntl' function. */
#define HAVE_FCNTL 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define to 1 if you have the `fdopen' function. */
#define HAVE_FDOPEN 1

/* Define to 1 if you have the <features.h> header file. */
#define HAVE_FEATURES_H 1

/* Define to 1 if you have the <floatingpoint.h> header file. */
/* #undef HAVE_FLOATINGPOINT_H */

/* Define if C's Rcomplex and Fortran's COMPLEX*16 can be interchanged, and
   can do arithmetic on the latter. */
#define HAVE_FORTRAN_DOUBLE_COMPLEX 1

/* Define to 1 if fseeko (and presumably ftello) exists and is declared. */
#define HAVE_FSEEKO 1

/* Define to 1 if you have the `ftello' function. */
#define HAVE_FTELLO 1

/* Define to 1 if you have the `ftruncate' function. */
#define HAVE_FTRUNCATE 1

/* Define to 1 if you have the `fwprintf' function. */
#define HAVE_FWPRINTF 1

/* Define to 1 if you have the `getcwd' function. */
#define HAVE_GETCWD 1

/* Define to 1 if you have the `getegid' function. */
#define HAVE_GETEGID 1

/* Define to 1 if you have the `geteuid' function. */
#define HAVE_GETEUID 1

/* Define to 1 if you have the `getgid' function. */
#define HAVE_GETGID 1

/* Define to 1 if you have the `getgrgid' function. */
#define HAVE_GETGRGID 1

/* Define to 1 if you have the `getline' function. */
#define HAVE_GETLINE 1

/* Define to 1 if you have the `getpagesize' function. */
#define HAVE_GETPAGESIZE 1

/* Define to 1 if you have the `getpriority' function. */
#define HAVE_GETPRIORITY 1

/* Define to 1 if you have the `getpwnam' function. */
#define HAVE_GETPWNAM 1

/* Define to 1 if you have the `getpwuid' function. */
#define HAVE_GETPWUID 1

/* Define to 1 if you have the `getrlimit' function. */
#define HAVE_GETRLIMIT 1

/* Define to 1 if you have the `getrusage' function. */
#define HAVE_GETRUSAGE 1

/* Define if the GNU gettext() function is already present or preinstalled. */
#define HAVE_GETTEXT 1

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the `getuid' function. */
#define HAVE_GETUID 1

/* Define if you have the GNU C library version >= 2. This is needed to fix a
   problem with getting the prototype of strptime(). */
#define HAVE_GLIBC2 1

/* Define to 1 if you have the `glob' function. */
#define HAVE_GLOB 1

/* Define to 1 if you have the <glob.h> header file. */
#define HAVE_GLOB_H 1

/* Define to 1 if you have the `gmtime_r' function. */
#define HAVE_GMTIME_R 1

/* Define to 1 if you have the <grp.h> header file. */
#define HAVE_GRP_H 1

/* Define to 1 if you have the `history_truncate_file' function. */
#define HAVE_HISTORY_TRUNCATE_FILE 1

/* Define to 1 if you have the `hypot' function. */
#define HAVE_HYPOT 1

/* Define if you have the iconv() function. */
#define HAVE_ICONV 1

/* Define if you have the `iconvlist' function. */
/* #undef HAVE_ICONVLIST */

/* Define if `iconv' accepts "CP1252". */
#define HAVE_ICONV_CP1252 1

/* Define to 1 if you have the <iconv.h> header file. */
#define HAVE_ICONV_H 1

/* Define to 1 if the system has the type `int64_t'. */
#define HAVE_INT64_T 1

/* Define if you have the 'intmax_t' type in <stdint.h> or <inttypes.h>. (For
   intl) */
#define HAVE_INTMAX_T 1

/* Define to 1 if the system has the type `intptr_t'. */
#define HAVE_INTPTR_T 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if <inttypes.h> exists, doesn't clash with <sys/types.h>, and
   declares uintmax_t. (For intl) */
#define HAVE_INTTYPES_H_WITH_UINTMAX 1

/* Define to 1 if the system has the type `int_fast64_t'. */
#define HAVE_INT_FAST64_T 1

/* Define to 1 if you have the `isblank' function. */
#define HAVE_ISBLANK 1

/* Define to 1 if you have the `isnan' function. */
#define HAVE_ISNAN 1

/* Define to 1 if you have the `iswctype' function. */
#define HAVE_ISWCTYPE 1

/* Define if you have the JPEG headers and libraries. */
/* #undef HAVE_JPEG */

/* Define if KERN_USRSTACK sysctl is supported. */
/* #undef HAVE_KERN_USRSTACK */

/* Define if you have KeySym defined in X11. */
#define HAVE_KEYSYM 1

/* Define to 1 if you have the `kill' function. */
#define HAVE_KILL 1

/* Define if you have <langinfo.h> and nl_langinfo(CODESET). */
#define HAVE_LANGINFO_CODESET 1

/* Define to 1 if you have the <langinfo.h> header file. */
#define HAVE_LANGINFO_H 1

/* Define if your <locale.h> file defines LC_MESSAGES. */
#define HAVE_LC_MESSAGES 1

/* Define if your system has libcurl >= 7.28.0 with support for https. */
#define HAVE_LIBCURL 1

/* Define if __libc_stack_end is visible. */
#define HAVE_LIBC_STACK_END 1

/* Define to 1 if you have the `dl' library (-ldl). */
#define HAVE_LIBDL 1

/* Define to 1 if you have the `icucore' library (-licucore). */
/* #undef HAVE_LIBICUCORE */

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the `ncurses' library (-lncurses). */
/* #undef HAVE_LIBNCURSES */

/* Define to 1 if you have the `readline' library (-lreadline). */
#define HAVE_LIBREADLINE 1

/* Define to 1 if you have the `rt' library (-lrt). */
#define HAVE_LIBRT 1

/* Define to 1 if you have the `sunmath' library (-lsunmath). */
/* #undef HAVE_LIBSUNMATH */

/* Define to 1 if you have the `termcap' library (-ltermcap). */
/* #undef HAVE_LIBTERMCAP */

/* Define to 1 if you have the `termlib' library (-ltermlib). */
/* #undef HAVE_LIBTERMLIB */

/* Define to 1 if you have the `tinfo' library (-ltinfo). */
/* #undef HAVE_LIBTINFO */

/* Define to 1 if you have the `tk' library (-ltk). */
/* #undef HAVE_LIBTK */

/* Define to 1 if you have the <limits.h> header file. */
#define HAVE_LIMITS_H 1

/* Define to 1 if you have the `link' function. */
#define HAVE_LINK 1

/* Define to 1 if you have the <locale.h> header file. */
#define HAVE_LOCALE_H 1

/* Define to 1 if you have the `localtime_r' function. */
#define HAVE_LOCALTIME_R 1

/* Define to 1 if you have the `log10' function. */
#define HAVE_LOG10 1

/* Define to 1 if you have the `log1p' function. */
#define HAVE_LOG1P 1

/* Define to 1 if you have the `log1pl' function. */
#define HAVE_LOG1PL 1

/* Define to 1 if you have the `log2' function. */
#define HAVE_LOG2 1

/* Define if you wish to use the 'long double' type. */
#define HAVE_LONG_DOUBLE 1

/* Define to 1 if the system has the type `long long int'. (For intl) */
#define HAVE_LONG_LONG_INT 1

/* Define if your system has lzma >= 5.0.3. */
#define HAVE_LZMA 1

/* Define to 1 if you have the <lzma.h> header file. */
#define HAVE_LZMA_H 1

/* Define to 1 if you have the `matherr' function. */
/* #undef HAVE_MATHERR */

/* Define to 1 if you have the `mbrtowc' function. */
#define HAVE_MBRTOWC 1

/* Define to 1 if the system has the type `mbstate_t'. */
#define HAVE_MBSTATE_T 1

/* Define to 1 if you have the `mbstowcs' function. */
#define HAVE_MBSTOWCS 1

/* Define to 1 if you have the `mempcpy' function. */
#define HAVE_MEMPCPY 1

/* Define to 1 if you have the <minix/config.h> header file. */
/* #undef HAVE_MINIX_CONFIG_H */

/* Define to 1 if you have the `mkdtemp' function. */
#define HAVE_MKDTEMP 1

/* Define to 1 if you have the `mkfifo' function. */
#define HAVE_MKFIFO 1

/* Define to 1 if you have a working `mmap' system call. */
#define HAVE_MMAP 1

/* Define to 1 if you have the `munmap' function. */
#define HAVE_MUNMAP 1

/* Define to 1 if you have the <ndir.h> header file, and it defines `DIR'. */
/* #undef HAVE_NDIR_H */

/* Define to 1 if you have the `nearbyint' function. */
#define HAVE_NEARBYINT 1

/* Define to 1 if you have the `nearbyintl' function. */
#define HAVE_NEARBYINTL 1

/* Define to 1 if you have the <netdb.h> header file. */
#define HAVE_NETDB_H 1

/* Define to 1 if you have the <netinet/in.h> header file. */
#define HAVE_NETINET_IN_H 1

/* Define to 1 if you have the `nl_langinfo' function. */
#define HAVE_NL_LANGINFO 1

/* Define if you have <langinfo.h> and it defines the NL_LOCALE_NAME macro if
   _GNU_SOURCE is defined. */
#define HAVE_NL_LOCALE_NAME 1

/* Define if module-loading does not need an underscore to be prepended to
   external names. */
/* #undef HAVE_NO_SYMBOL_UNDERSCORE */

/* Define if you have off_t, fseeko and ftello. */
#define HAVE_OFF_T 1

/* Define if you have C OpenMP support. */
#define HAVE_OPENMP 1

/* Define if your OpenMP 4 implementation fully supports SIMD reduction */
#define HAVE_OPENMP_SIMDRED 1

/* Define to 1 if you have pangocairo. */
/* #undef HAVE_PANGOCAIRO */

/* Define if your system has pcre2. */
#define HAVE_PCRE2 1

/* Define to 1 if you have the <pcre.h> header file. */
/* #undef HAVE_PCRE_H */

/* Define to 1 if you have the <pcre/pcre.h> header file. */
/* #undef HAVE_PCRE_PCRE_H */

/* Define if you have the PNG headers and libraries. */
#define HAVE_PNG 1

/* Define to 1 if you have the `popen' function. */
#define HAVE_POPEN 1

/* Define if your system time functions do not count leap seconds, as required
   by POSIX. */
#define HAVE_POSIX_LEAPSECONDS 1

/* Define if your printf() function supports format strings with positions.
   (For intl) */
#define HAVE_POSIX_PRINTF 1

/* Define if you have POSIX.1 compatible sigsetjmp/siglongjmp. */
#define HAVE_POSIX_SETJMP 1

/* Define to 1 if you have the `powl' function. */
#define HAVE_POWL 1

/* Define to 1 if you have the `pown' function. */
/* #undef HAVE_POWN */

/* Define if have support for POSIX threads. */
#define HAVE_PTHREAD 1

/* Define if the <pthread.h> defines PTHREAD_MUTEX_RECURSIVE. (For intl) */
/* #undef HAVE_PTHREAD_MUTEX_RECURSIVE */

/* Define if the POSIX multithreading library has read/write locks. (For intl)
   */
/* #undef HAVE_PTHREAD_RWLOCK */

/* Define to 1 if you have the `putenv' function. */
#define HAVE_PUTENV 1

/* Define if putenv("FOO") can unset an environment variable */
#define HAVE_PUTENV_UNSET 1

/* Define if putenv("FOO=") can unset an environment variable */
/* #undef HAVE_PUTENV_UNSET2 */

/* Define to 1 if you have the <pwd.h> header file. */
#define HAVE_PWD_H 1

/* Define to 1 if you have the <readline/history.h> header file. */
#define HAVE_READLINE_HISTORY_H 1

/* Define to 1 if you have the <readline/readline.h> header file. */
#define HAVE_READLINE_READLINE_H 1

/* Define to 1 if you have the `readlink' function. */
#define HAVE_READLINK 1

/* Define to 1 if you have the `realpath' function. */
#define HAVE_REALPATH 1

/* Define to 1 if you have the `rintl' function. */
#define HAVE_RINTL 1

/* Define to 1 if you have the `rl_completion_matches' function. */
#define HAVE_RL_COMPLETION_MATCHES 1

/* Define to 1 if you have the `rl_resize_terminal' function. */
#define HAVE_RL_RESIZE_TERMINAL 1

/* Define to 1 if you have the `rl_sort_completion_matches' function. */
#define HAVE_RL_SORT_COMPLETION_MATCHES 1

/* Define to 1 if you have the `sched_getaffinity' function. */
#define HAVE_SCHED_GETAFFINITY 1

/* Define to 1 if you have the <sched.h> header file. */
#define HAVE_SCHED_H 1

/* Define to 1 if you have the `sched_setaffinity' function. */
#define HAVE_SCHED_SETAFFINITY 1

/* Define to 1 if you have the `select' function. */
#define HAVE_SELECT 1

/* Define to 1 if you have the `setenv' function. */
#define HAVE_SETENV 1

/* Define to 1 if you have the `setitimer' function. */
#define HAVE_SETITIMER 1

/* Define to 1 if you have the `setlocale' function. */
#define HAVE_SETLOCALE 1

/* Define to 1 if you have the `setrlimit' function. */
#define HAVE_SETRLIMIT 1

/* Define to 1 if you have the `sigaction' function. */
#define HAVE_SIGACTION 1

/* Define to 1 if you have the `sigaltstack' function. */
#define HAVE_SIGALTSTACK 1

/* Define to 1 if you have the `sigemptyset' function. */
#define HAVE_SIGEMPTYSET 1

/* Define to 1 if you have the `sinpi' function. */
/* #undef HAVE_SINPI */

/* Define to 1 if you have the `snprintf' function. */
#define HAVE_SNPRINTF 1

/* Define to 1 if the system has the type `stack_t'. */
#define HAVE_STACK_T 1

/* Define to 1 if you have the `stat' function. */
#define HAVE_STAT 1

/* Define to 1 if you have the <stdalign.h> header file. */
#define HAVE_STDALIGN_H 1

/* Define to 1 if you have the <stdarg.h> header file. */
#define HAVE_STDARG_H 1

/* Define to 1 if you have the <stdbool.h> header file. */
#define HAVE_STDBOOL_H 1

/* Define to 1 if you have the <stddef.h> header file. */
#define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define if <stdint.h> exists, doesn't clash with <sys/types.h>, and declares
   uintmax_t. (For intl) */
#define HAVE_STDINT_H_WITH_UINTMAX 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `stpcpy' function. */
#define HAVE_STPCPY 1

/* Define to 1 if you have the `strcasecmp' function. */
#define HAVE_STRCASECMP 1

/* Define to 1 if you have the `strdup' function. */
#define HAVE_STRDUP 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strncasecmp' function. */
#define HAVE_STRNCASECMP 1

/* Define to 1 if you have the `strtoul' function. */
#define HAVE_STRTOUL 1

/* Define to 1 if `st_atimensec' is a member of `struct stat'. */
/* #undef HAVE_STRUCT_STAT_ST_ATIMENSEC */

/* Define to 1 if `st_atimespec.tv_nsec' is a member of `struct stat'. */
/* #undef HAVE_STRUCT_STAT_ST_ATIMESPEC_TV_NSEC */

/* Define to 1 if `st_atim.st__tim.tv_nsec' is a member of `struct stat'. */
/* #undef HAVE_STRUCT_STAT_ST_ATIM_ST__TIM_TV_NSEC */

/* Define to 1 if `st_atim.tv_nsec' is a member of `struct stat'. */
#define HAVE_STRUCT_STAT_ST_ATIM_TV_NSEC 1

/* Define to 1 if you have the <sunmath.h> header file. */
/* #undef HAVE_SUNMATH_H */

/* Define to 1 if you have the `symlink' function. */
#define HAVE_SYMLINK 1

/* Define to 1 if you have the `sysconf' function. */
#define HAVE_SYSCONF 1

/* Define to 1 if you have the <sys/dir.h> header file, and it defines `DIR'.
   */
/* #undef HAVE_SYS_DIR_H */

/* Define to 1 if you have the <sys/ndir.h> header file, and it defines `DIR'.
   */
/* #undef HAVE_SYS_NDIR_H */

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/resource.h> header file. */
#define HAVE_SYS_RESOURCE_H 1

/* Define to 1 if you have the <sys/select.h> header file. */
#define HAVE_SYS_SELECT_H 1

/* Define to 1 if you have the <sys/socket.h> header file. */
#define HAVE_SYS_SOCKET_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/times.h> header file. */
#define HAVE_SYS_TIMES_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <sys/utsname.h> header file. */
#define HAVE_SYS_UTSNAME_H 1

/* Define to 1 if you have <sys/wait.h> that is POSIX.1 compatible. */
#define HAVE_SYS_WAIT_H 1

/* Define to 1 if you have the `tanpi' function. */
/* #undef HAVE_TANPI */

/* Define if you have the Tcl/Tk headers and libraries and want Tcl/Tk support
   to be built. */
/* #undef HAVE_TCLTK */

/* Define to 1 if you have the <thread.h> header file. */
/* #undef HAVE_THREAD_H */

/* Define to 1 if you have the `thr_stksegment' function. */
/* #undef HAVE_THR_STKSEGMENT */

/* Define this if libtiff is available. */
/* #undef HAVE_TIFF */

/* Define to 1 if you have the <tiffio.h> header file. */
/* #undef HAVE_TIFFIO_H */

/* Define to 1 if you have the `tilde_expand_word' function. */
#define HAVE_TILDE_EXPAND_WORD 1

/* Define to 1 if you have the `times' function. */
#define HAVE_TIMES 1

/* Define to 1 if you have the `timespec_get' function. */
#define HAVE_TIMESPEC_GET 1

/* Define to 1 if your 'struct tm' has tm_gmtoff. */
#define HAVE_TM_GMTOFF 1

/* Define to 1 if your 'struct tm' has tm_zone. */
#define HAVE_TM_ZONE 1

/* Define if your system has tre. */
/* #undef HAVE_TRE */

/* Define to 1 if you have the <tre/tre.h> header file. */
/* #undef HAVE_TRE_TRE_H */

/* Define to 1 if you have the `tsearch' function. */
#define HAVE_TSEARCH 1

/* Define if you have the 'uintmax_t' type in <stdint.h> or <inttypes.h>. (For
   intl) */
#define HAVE_UINTMAX_T 1

/* Define to 1 if the system has the type `uintptr_t'. */
#define HAVE_UINTPTR_T 1

/* Define to 1 if you have the `umask' function. */
#define HAVE_UMASK 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the `unsetenv' function. */
#define HAVE_UNSETENV 1

/* Define if you have the 'unsigned long long' type. (For intl) */
#define HAVE_UNSIGNED_LONG_LONG 1

/* Define to 1 if the system has the type `unsigned long long int'. (For intl)
   */
#define HAVE_UNSIGNED_LONG_LONG_INT 1

/* Define to 1 if you have the `utime' function. */
#define HAVE_UTIME 1

/* Define to 1 if you have the `utimensat' function. */
#define HAVE_UTIMENSAT 1

/* Define to 1 if you have the `utimes' function. */
#define HAVE_UTIMES 1

/* Define to 1 if you have the <utime.h> header file. */
#define HAVE_UTIME_H 1

/* Define to 1 if you have the <valgrind/memcheck.h> header file. */
/* #undef HAVE_VALGRIND_MEMCHECK_H */

/* Define to 1 if you have the `vasprintf' function. */
#define HAVE_VASPRINTF 1

/* Define to 1 if you have the `va_copy' function. */
#define HAVE_VA_COPY 1

/* Define to 1 or 0, depending whether the compiler supports simple visibility
   declarations. (For intl) */
#define HAVE_VISIBILITY 1

/* Define to 1 if __attribute__((visibility())) is supported */
#define HAVE_VISIBILITY_ATTRIBUTE 1

/* Define to 1 if you have the <wchar.h> header file. */
#define HAVE_WCHAR_H 1

/* Define if you have the 'wchar_t' type. (For intl) */
#define HAVE_WCHAR_T 1

/* Define to 1 if you have the `wcrtomb' function. */
#define HAVE_WCRTOMB 1

/* Define to 1 if you have the `wcscoll' function. */
#define HAVE_WCSCOLL 1

/* Define to 1 if you have the `wcsftime' function. */
#define HAVE_WCSFTIME 1

/* Define to 1 if you have the `wcslen' function. */
#define HAVE_WCSLEN 1

/* Define to 1 if you have the `wcstod' function. */
#define HAVE_WCSTOD 1

/* Define to 1 if you have the `wcstombs' function. */
#define HAVE_WCSTOMBS 1

/* Define to 1 if you have the `wcswidth' function. */
#define HAVE_WCSWIDTH 1

/* Define to 1 if you have the `wctrans' function. */
#define HAVE_WCTRANS 1

/* Define to 1 if the system has the type `wctrans_t'. */
#define HAVE_WCTRANS_T 1

/* Define to 1 if you have the `wctype' function. */
#define HAVE_WCTYPE 1

/* Define to 1 if you have the <wctype.h> header file. */
#define HAVE_WCTYPE_H 1

/* Define to 1 if you have the `wcwidth' function. */
#define HAVE_WCWIDTH 1

/* Define if you have the 'wint_t' type. (For intl) */
#define HAVE_WINT_T 1

/* Define to 1 if you have cairo. */
/* #undef HAVE_WORKING_CAIRO */

/* Define if calloc(0) returns a null pointer. */
#define HAVE_WORKING_CALLOC 1

/* Define if ctanh() exists and is working correctly. */
#define HAVE_WORKING_CTANH 1

/* Define if your ftell works correctly on files opened for append. */
#define HAVE_WORKING_FTELL 1

/* Define if isfinite() is correct for -Inf/NaN/Inf. */
#define HAVE_WORKING_ISFINITE 1

/* Define if log1p() exists and is accurate enough. */
#define HAVE_WORKING_LOG1P 1

/* Define if your mktime works correctly after 2037. */
#define HAVE_WORKING_MKTIME_AFTER_2037 1

/* Define if your mktime works correctly before 1902. */
#define HAVE_WORKING_MKTIME_BEFORE_1902 1

/* Define if your mktime works correctly before 1970. */
#define HAVE_WORKING_MKTIME_BEFORE_1970 1

/* Define if sigaction() is complete enough for R's usage */
#define HAVE_WORKING_SIGACTION 1

/* Define to 1 if you have cairo including Xlib. */
/* #undef HAVE_WORKING_X11_CAIRO */

/* Define if you have the X11 headers and libraries, and want the X11 GUI to
   be built. */
#define HAVE_X11 1

/* Define if you have the X11/Xmu headers and libraries. */
#define HAVE_X11_Xmu 1

/* Define to 1 if you have the `__cospi' function. */
/* #undef HAVE___COSPI */

/* Define to 1 if you have the `__fsetlocking' function. */
#define HAVE___FSETLOCKING 1

/* Define to 1 if you have the `__sinpi' function. */
/* #undef HAVE___SINPI */

/* Define to 1 if you have the `__tanpi' function. */
/* #undef HAVE___TANPI */

/* Define as const if the declaration of iconv() needs const. */
#define ICONV_CONST 

/* Define if you have IEEE 754 floating point arithmetic. */
#define IEEE_754 1

/* Define if integer division by zero raises signal SIGFPE. (For intl) */
#define INTDIV0_RAISES_SIGFPE 0

/* Define if you have 32 bit ints. */
#define INT_32_BITS 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Define if mktime sets errno. */
#define MKTIME_SETS_ERRNO 1

/* Define to disable Valgrind instrumentation */
#define NVALGRIND 1

/* Define if using GNU-style Objective C runtime. */
/* #undef OBJC_GNU_RUNTIME */

/* Define if using NeXT/Apple-style Objective C runtime. */
/* #undef OBJC_NEXT_RUNTIME */

/* Name of package */
#define PACKAGE "R"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "https://bugs.r-project.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "R"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "R 4.3.3"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "R"

/* Define to the home page for this package. */
#define PACKAGE_URL "https://www.r-project.org"

/* Define to the version of this package. */
#define PACKAGE_VERSION "4.3.3"

/* PCRE2 code unit width wanted. */
#define PCRE2_CODE_UNIT_WIDTH 8

/* Define if <inttypes.h> exists and defines unusable PRI* macros. (For intl)
   */
/* #undef PRI_MACROS_BROKEN */

/* Define if the pthread_in_use() detection is hard. (For intl) */
/* #undef PTHREAD_IN_USE_DETECTION_HARD */

/* Define this to use architecture-dependent subdirectories of this name. */
#define R_ARCH ""

/* Define this to be the name of the CPU of your system. */
#define R_CPU "x86_64"

/* Define as `inline', or `__inline__' or `__inline' if that's what the C
   compiler calls it, or to nothing if it is not supported. */
#define R_INLINE inline

/* Define this to enable memory profiling. */
#define R_MEMORY_PROFILING 1

/* Define this to be the name of the OS of your system. */
#define R_OS "linux-gnu"

/* Define this to be the canonical name (cpu-vendor-os) of your system. */
#define R_PLATFORM "x86_64-pc-linux-gnu"

/* Define this to be printing command on your system. */
#define R_PRINTCMD "lpr"

/* Define this to enable R-level profiling. */
#define R_PROFILING 1

/* Type for socket lengths: socklen_t, sock_t, int? */
#define R_SOCKLEN_T socklen_t

/* Define this to be the name of the vendor of your system. */
#define R_VENDOR "pc"

/* Define this to be the extension used for shared objects on your system. */
#define SHLIB_EXT ".so"

/* The size of `double', as computed by sizeof. */
#define SIZEOF_DOUBLE 8

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG 8

/* The size of `long double', as computed by sizeof. */
#define SIZEOF_LONG_DOUBLE 16

/* The size of `long long', as computed by sizeof. */
#define SIZEOF_LONG_LONG 8

/* The size of `size_t', as computed by sizeof. */
#define SIZEOF_SIZE_T 8

/* Define as the maximum value of type 'size_t', if the system doesn't define
   it. (For intl) */
/* #undef SIZE_MAX */

/* If using the C implementation of alloca, define if you know the
   direction of stack growth for your system; otherwise it will be
   automatically deduced at runtime.
	STACK_DIRECTION > 0 => grows toward higher addresses
	STACK_DIRECTION < 0 => grows toward lower addresses
	STACK_DIRECTION = 0 => direction of growth unknown */
/* #undef STACK_DIRECTION */

/* Define to 1 if all of the C90 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it. */
#define STDC_HEADERS 1

/* Define if you have C/C++/Fortran OpenMP support for package code. */
#define SUPPORT_OPENMP 1

/* Define to enable provoking compile errors on write barrier violation. */
/* #undef TESTING_WRITE_BARRIER */

/* Define to 1 if the type of the st_atim member of a struct stat is struct
   timespec. */
#define TYPEOF_STRUCT_STAT_ST_ATIM_IS_STRUCT_TIMESPEC 1

/* Define to use ICU for collation. */
/* #undef USE_ICU */

/* Define to use Apple's ICU. */
/* #undef USE_ICU_APPLE */

/* Define to use internal time-zone code */
/* #undef USE_INTERNAL_MKTIME */

/* Define if the POSIX multithreading library can be used. (For intl) */
/* #undef USE_POSIX_THREADS */

/* Define if references to the POSIX multithreading library should be made
   weak. (For intl) */
/* #undef USE_POSIX_THREADS_WEAK */

/* Define if the GNU Pth multithreading library can be used. (For intl) */
/* #undef USE_PTH_THREADS */

/* Define if references to the GNU Pth multithreading library should be made
   weak. (For intl) */
/* #undef USE_PTH_THREADS_WEAK */

/* Define to 1 to use internal code for `towlower' and `towupper'. */
/* #undef USE_RI18N_CASE */

/* Define to 1 to use internal `iswprint' etc. */
/* #undef USE_RI18N_FNS */

/* Define to 1 to use internal `wcwidth' */
#define USE_RI18N_WIDTH 1

/* Define if the old Solaris multithreading library can be used. (For intl) */
/* #undef USE_SOLARIS_THREADS */

/* Define if references to the old Solaris multithreading library should be
   made weak. (For intl) */
/* #undef USE_SOLARIS_THREADS_WEAK */

/* Enable extensions on AIX 3, Interix.  */
#ifndef _ALL_SOURCE
# define _ALL_SOURCE 1
#endif
/* Enable general extensions on macOS.  */
#ifndef _DARWIN_C_SOURCE
# define _DARWIN_C_SOURCE 1
#endif
/* Enable general extensions on Solaris.  */
#ifndef __EXTENSIONS__
# define __EXTENSIONS__ 1
#endif
/* Enable GNU extensions on systems that have them.  */
#ifndef _GNU_SOURCE
# define _GNU_SOURCE 1
#endif
/* Enable X/Open compliant socket functions that do not require linking
   with -lxnet on HP-UX 11.11.  */
#ifndef _HPUX_ALT_XOPEN_SOCKET_API
# define _HPUX_ALT_XOPEN_SOCKET_API 1
#endif
/* Identify the host operating system as Minix.
   This macro does not affect the system headers' behavior.
   A future release of Autoconf may stop defining this macro.  */
#ifndef _MINIX
/* # undef _MINIX */
#endif
/* Enable general extensions on NetBSD.
   Enable NetBSD compatibility extensions on Minix.  */
#ifndef _NETBSD_SOURCE
# define _NETBSD_SOURCE 1
#endif
/* Enable OpenBSD compatibility extensions on NetBSD.
   Oddly enough, this does nothing on OpenBSD.  */
#ifndef _OPENBSD_SOURCE
# define _OPENBSD_SOURCE 1
#endif
/* Define to 1 if needed for POSIX-compatible behavior.  */
#ifndef _POSIX_SOURCE
/* # undef _POSIX_SOURCE */
#endif
/* Define to 2 if needed for POSIX-compatible behavior.  */
#ifndef _POSIX_1_SOURCE
/* # undef _POSIX_1_SOURCE */
#endif
/* Enable POSIX-compatible threading on Solaris.  */
#ifndef _POSIX_PTHREAD_SEMANTICS
# define _POSIX_PTHREAD_SEMANTICS 1
#endif
/* Enable extensions specified by ISO/IEC TS 18661-5:2014.  */
#ifndef __STDC_WANT_IEC_60559_ATTRIBS_EXT__
# define __STDC_WANT_IEC_60559_ATTRIBS_EXT__ 1
#endif
/* Enable extensions specified by ISO/IEC TS 18661-1:2014.  */
#ifndef __STDC_WANT_IEC_60559_BFP_EXT__
# define __STDC_WANT_IEC_60559_BFP_EXT__ 1
#endif
/* Enable extensions specified by ISO/IEC TS 18661-2:2015.  */
#ifndef __STDC_WANT_IEC_60559_DFP_EXT__
# define __STDC_WANT_IEC_60559_DFP_EXT__ 1
#endif
/* Enable extensions specified by ISO/IEC TS 18661-4:2015.  */
#ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
# define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1
#endif
/* Enable extensions specified by ISO/IEC TS 18661-3:2015.  */
#ifndef __STDC_WANT_IEC_60559_TYPES_EXT__
# define __STDC_WANT_IEC_60559_TYPES_EXT__ 1
#endif
/* Enable extensions specified by ISO/IEC TR 24731-2:2010.  */
#ifndef __STDC_WANT_LIB_EXT2__
# define __STDC_WANT_LIB_EXT2__ 1
#endif
/* Enable extensions specified by ISO/IEC 24747:2009.  */
#ifndef __STDC_WANT_MATH_SPEC_FUNCS__
# define __STDC_WANT_MATH_SPEC_FUNCS__ 1
#endif
/* Enable extensions on HP NonStop.  */
#ifndef _TANDEM_SOURCE
# define _TANDEM_SOURCE 1
#endif
/* Enable X/Open extensions.  Define to 500 only if necessary
   to make mbstate_t available.  */
#ifndef _XOPEN_SOURCE
/* # undef _XOPEN_SOURCE */
#endif


/* Define if the Win32 multithreading API can be used. (For intl) */
/* #undef USE_WIN32_THREADS */

/* Define according to your operating system type. */
#define Unix 1

/* Define as 1 or 2 to specify levels of Valgrind instrumentation */
#define VALGRIND_LEVEL 0

/* Version number of package */
#define VERSION "4.3.3"

/* Define to 1 if your processor stores words with the most significant byte
   first (like Motorola and SPARC, unlike Intel and VAX). */
/* #undef WORDS_BIGENDIAN */

/* Define according to your operating system type. */
/* #undef Win32 */

/* Define to 1 if the X Window System is missing or not being used. */
/* #undef X_DISPLAY_MISSING */

/* Number of bits in a file offset, on hosts where this is settable. */
/* #undef _FILE_OFFSET_BITS */

/* Define to 1 to make fseeko visible on some hosts (e.g. glibc 2.2). */
/* #undef _LARGEFILE_SOURCE */

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */

/* Define for Solaris 2.5.1 so the uint64_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT64_T */

/* Define to 'long' if <sys/types.h> does not define. Apparently necessary to
   fix a GCC bug on AIX? */
/* #undef blkcnt_t */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define as a signed integer type capable of holding a process identifier. */
/* #undef pid_t */

/* Define as the type of the result of subtracting two pointers, if the system
   doesn't define it. (For intl) */
/* #undef ptrdiff_t */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to the type of an unsigned integer type of width exactly 64 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint64_t */

/* Define to unsigned long or unsigned long long if <stdint.h> and
   <inttypes.h> don't define. (For intl) */
/* #undef uintmax_t */


#endif /* not R_CONFIG_H */


#define __libc_lock_t                   gl_lock_t
#define __libc_lock_define              gl_lock_define
#define __libc_lock_define_initialized  gl_lock_define_initialized
#define __libc_lock_init                gl_lock_init
#define __libc_lock_lock                gl_lock_lock
#define __libc_lock_unlock              gl_lock_unlock
#define __libc_lock_recursive_t                   gl_recursive_lock_t
#define __libc_lock_define_recursive              gl_recursive_lock_define
#define __libc_lock_define_initialized_recursive  gl_recursive_lock_define_initialized
#define __libc_lock_init_recursive                gl_recursive_lock_init
#define __libc_lock_lock_recursive                gl_recursive_lock_lock
#define __libc_lock_unlock_recursive              gl_recursive_lock_unlock
#define glthread_in_use  libintl_thread_in_use
#define glthread_lock_init     libintl_lock_init
#define glthread_lock_lock     libintl_lock_lock
#define glthread_lock_unlock   libintl_lock_unlock
#define glthread_lock_destroy  libintl_lock_destroy
#define glthread_rwlock_init     libintl_rwlock_init
#define glthread_rwlock_rdlock   libintl_rwlock_rdlock
#define glthread_rwlock_wrlock   libintl_rwlock_wrlock
#define glthread_rwlock_unlock   libintl_rwlock_unlock
#define glthread_rwlock_destroy  libintl_rwlock_destroy
#define glthread_recursive_lock_init     libintl_recursive_lock_init
#define glthread_recursive_lock_lock     libintl_recursive_lock_lock
#define glthread_recursive_lock_unlock   libintl_recursive_lock_unlock
#define glthread_recursive_lock_destroy  libintl_recursive_lock_destroy
#define glthread_once                 libintl_once
#define glthread_once_call            libintl_once_call
#define glthread_once_singlethreaded  libintl_once_singlethreaded

