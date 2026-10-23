#define _GNU_SOURCE
#include <dlfcn.h>
#include <fcntl.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

/* Redirect stale build-time ROOT/GenFit header paths to CVMFS installs. */
static const char kOldGenFit[] = "/gpfs/mnt/gpfs02/sphenix/user/phnxbld/workarea/sPHENIX_SL7.3/almalinux-9.2/gcc-14.2.0/need_root_version/root-6.32.06/GenFit-02-00-05/";
static const char kGenFitInclude[] = "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/genfit2-02-00-05_root-6.32.06/include/";
static const char kOldRootBuild[] = "/gpfs/mnt/gpfs02/sphenix/user/phnxbld/workarea/sPHENIX_SL7.3/almalinux-9.2/gcc-14.2.0/root-6.32.06/root-6.32.06/build/include/";
static const char kGxxInclude[] = "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/gcc/14.2.0-2f0a0/x86_64-el9/include/c++/14.2.0/";

static const char* basename_of(const char* path) { const char* p = strrchr(path, '/'); return p ? p + 1 : path; }
static const char* redirect_path(const char* path, char* out, size_t size)
{
  if (!path || path[0] != '/') return path;
  if (strncmp(path, kOldGenFit, sizeof(kOldGenFit) - 1) == 0) { snprintf(out, size, "%s%s", kGenFitInclude, basename_of(path)); return out; }
  if (strncmp(path, kOldRootBuild, sizeof(kOldRootBuild) - 1) == 0) { snprintf(out, size, "%s%s", kGxxInclude, basename_of(path)); return out; }
  return path;
}

#define RESOLVE(symbol) static __typeof__(&symbol) real_##symbol = NULL; if (!real_##symbol) real_##symbol = (__typeof__(&symbol)) dlsym(RTLD_NEXT, #symbol)

int open(const char* p, int f, ...) { RESOLVE(open); char b[PATH_MAX]; const char* q=redirect_path(p,b,sizeof(b)); if(f&O_CREAT){va_list a;va_start(a,f);mode_t m=(mode_t)va_arg(a,int);va_end(a);return real_open(q,f,m);}return real_open(q,f); }
int open64(const char* p, int f, ...) { RESOLVE(open64); char b[PATH_MAX]; const char* q=redirect_path(p,b,sizeof(b)); if(f&O_CREAT){va_list a;va_start(a,f);mode_t m=(mode_t)va_arg(a,int);va_end(a);return real_open64(q,f,m);}return real_open64(q,f); }
int openat(int d,const char* p,int f,...) { RESOLVE(openat); char b[PATH_MAX];const char*q=redirect_path(p,b,sizeof(b));if(f&O_CREAT){va_list a;va_start(a,f);mode_t m=(mode_t)va_arg(a,int);va_end(a);return real_openat(d,q,f,m);}return real_openat(d,q,f); }
int openat64(int d,const char* p,int f,...) { RESOLVE(openat64); char b[PATH_MAX];const char*q=redirect_path(p,b,sizeof(b));if(f&O_CREAT){va_list a;va_start(a,f);mode_t m=(mode_t)va_arg(a,int);va_end(a);return real_openat64(d,q,f,m);}return real_openat64(d,q,f); }
int stat(const char*p,struct stat*s){RESOLVE(stat);char b[PATH_MAX];return real_stat(redirect_path(p,b,sizeof(b)),s);}
int lstat(const char*p,struct stat*s){RESOLVE(lstat);char b[PATH_MAX];return real_lstat(redirect_path(p,b,sizeof(b)),s);}
int access(const char*p,int m){RESOLVE(access);char b[PATH_MAX];return real_access(redirect_path(p,b,sizeof(b)),m);}
int faccessat(int d,const char*p,int m,int f){RESOLVE(faccessat);char b[PATH_MAX];return real_faccessat(d,redirect_path(p,b,sizeof(b)),m,f);}
FILE* fopen(const char*p,const char*m){RESOLVE(fopen);char b[PATH_MAX];return real_fopen(redirect_path(p,b,sizeof(b)),m);}
FILE* fopen64(const char*p,const char*m){RESOLVE(fopen64);char b[PATH_MAX];return real_fopen64(redirect_path(p,b,sizeof(b)),m);}
